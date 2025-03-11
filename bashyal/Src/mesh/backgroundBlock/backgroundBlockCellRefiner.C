#include "backgroundBlock.H"
#include "foamCGALConverter.H"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h>
#include <map>
#include <vector>

// Ensure labelList is available
// #include <labelList.H>

namespace Bashyal
{

    // Define CGAL types
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
    typedef CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> Nef_polyhedron;
    typedef CGAL::Plane_3<Kernel> Plane_3;
    typedef FoamCGALConverter<Kernel> Converter;

    void backgroundBlock::decomposeToConvex()
    {
        // Ensure faces are triangulated for CGAL compatibility
        this->triangulateFaces();

        // Convert OpenFOAM geometry to CGAL Polyhedron_3 using FoamCGALConverter
        Polyhedron originalPolyhedron;
        try
        {
            originalPolyhedron = Converter::toCGALPolyhedron(points_, triFaces_);
        }
        catch (const std::exception &e)
        {
            FatalErrorInFunction << "Failed to convert to CGAL Polyhedron: " << e.what() << exit(Foam::FatalError);
        }

        // Validate the constructed polyhedron
        if (!originalPolyhedron.is_valid())
        {
            FatalErrorInFunction << "Invalid polyhedron construction" << exit(Foam::FatalError);
        }

        // Convert to Nef_polyhedron_3 for convex decomposition
        Nef_polyhedron N(originalPolyhedron);
        CGAL::convex_decomposition_3(N);

        // Extract convex polyhedra from the decomposition
        std::vector<Polyhedron> convexPolyhedra;
        Nef_polyhedron::Volume_const_iterator ci = ++N.volumes_begin();
        for (; ci != N.volumes_end(); ++ci)
        {
            if (ci->mark())
            {
                Polyhedron cp;
                N.convert_inner_shell_to_polyhedron(ci->shells_begin(), cp);
                convexPolyhedra.push_back(cp);
            }
        }

        Foam::label nCells = static_cast<Foam::label>(convexPolyhedra.size());
        if (nCells == 0)
        {
            FatalErrorInFunction << "No convex parts found" << exit(Foam::FatalError);
        }

        // Collect unique points from all convex polyhedra
        Foam::pointField allPoints;
        Foam::HashTable<Foam::label, Foam::point> pointMap;
        for (const Polyhedron &cp : convexPolyhedra)
        {
            for (Polyhedron::Vertex_const_iterator vi = cp.vertices_begin(); vi != cp.vertices_end(); ++vi)
            {
                Foam::point pt = Converter::toFoamPoint(vi->point());
                if (!pointMap.found(pt))
                {
                    pointMap.insert(pt, allPoints.size());
                    allPoints.append(pt);
                }
            }
        }

        // Build face ownership map to distinguish internal and boundary faces
        struct FaceKey : public std::vector<Foam::label>
        {
            FaceKey(const std::vector<Foam::label> &v) : std::vector<Foam::label>(v) { std::sort(begin(), end()); }
            bool operator<(const FaceKey &other) const { return std::lexicographical_compare(begin(), end(), other.begin(), other.end()); }
        };
        std::map<FaceKey, std::vector<std::pair<Foam::label, Foam::face>>> faceMap;

        for (Foam::label cellI = 0; cellI < nCells; ++cellI)
        {
            const Polyhedron &cp = convexPolyhedra[cellI];
            for (Polyhedron::Facet_const_iterator fi = cp.facets_begin(); fi != cp.facets_end(); ++fi)
            {
                std::vector<Foam::label> pointIndices;
                auto circ = fi->facet_begin();
                do
                {
                    Foam::point pt = Converter::toFoamPoint(circ->vertex()->point());
                    pointIndices.push_back(pointMap[pt]);
                    ++circ;
                } while (circ != fi->facet_begin());

                Foam::face f(pointIndices.size());
                for (Foam::label i = 0; i < pointIndices.size(); ++i)
                {
                    f[i] = pointIndices[i];
                }
                FaceKey key(pointIndices);
                faceMap[key].push_back(std::make_pair(cellI, f));
            }
        }

        // Separate boundary faces and construct topology
        Foam::faceList newFaces;
        Foam::labelList newOwners;
        Foam::labelList newNeighbours;
        Foam::faceList boundaryFaces;
        Foam::label newNboundaries = 0;

        for (const auto &entry : faceMap)
        {
            const FaceKey &key = entry.first;
            const auto &instances = entry.second;
            if (instances.size() == 1)
            {
                // Boundary face
                Foam::label faceI = newFaces.size();
                newFaces.append(instances[0].second);
                newOwners.append(instances[0].first);
                newNeighbours.append(-1);
                boundaryFaces.append(instances[0].second);
                newNboundaries++;
            }
            else if (instances.size() == 2)
            {
                // Internal face
                Foam::label owner = instances[0].first;
                Foam::label neighbour = instances[1].first;
                Foam::label faceI = newFaces.size();
                newFaces.append(instances[0].second);
                newOwners.append(std::min(owner, neighbour));
                newNeighbours.append(std::max(owner, neighbour));
            }
            else
            {
                FatalErrorInFunction << "Face shared by >2 cells" << exit(Foam::FatalError);
            }
        }

        // Generate new patches
        Foam::List<int> newPatches;
        mapNewFacesToBoundaries(
            allPoints, newFaces, newNeighbours, newNboundaries, // New polyhedron data
            newPatches                                          // Output patches
        );

        // Update class attributes with new topology
        points_ = allPoints;
        faces_ = newFaces;
        owners_ = newOwners;
        neighbours_ = newNeighbours;
        patches_ = newPatches;
        ncells_ = nCells;
        edited_ = true;
        if (nCells > 1)
        {
            multiple_ = true;
        }
    }

    void backgroundBlock::mapNewFacesToBoundaries(
        const Foam::pointField &newPoints,
        const Foam::faceList &newFaces,
        const Foam::labelList &newNeighbours,
        const Foam::label newNboundaries,
        Foam::List<int> &newPatches)
    {
        // Define CGAL kernel
        typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

        // Define a tolerance for floating-point comparisons
        const double tolerance = 1e-10;

        // Helper function to compute the plane of a face from its points
        auto computePlane = [](const Foam::face &face, const Foam::pointField &points) -> CGAL::Plane_3<Kernel>
        {
            if (face.size() < 3)
            {
                throw std::invalid_argument("Face must have at least 3 points to define a plane.");
            }
            CGAL::Point_3<Kernel> p0 = Converter::toCGALPoint(points[face[0]]);
            CGAL::Point_3<Kernel> p1 = Converter::toCGALPoint(points[face[1]]);
            CGAL::Point_3<Kernel> p2 = Converter::toCGALPoint(points[face[2]]);
            return CGAL::Plane_3<Kernel>(p0, p1, p2); // Plane defined by three points
        };

        // Step 1: Compute planes for all faces of polyhedron A
        std::vector<CGAL::Plane_3<Kernel>> planes(nboundaries_);
        Foam::List<int> boundaryPatches(nboundaries_);
        Foam::label loopCount;
        for (Foam::label i = 0; i < faces_.size(); ++i)
        {
            if (neighbours_[i] == -1)
            {
                planes[loopCount] = computePlane(faces_[i], points_);
                boundaryPatches[loopCount] = patches_[i];
                loopCount++;
            }
        }

        // Step 3: Resize newPatches to match the number of new faces and initialize with a default value
        newPatches.setSize(newFaces.size());
        forAll(newPatches, i)
        {
            newPatches[i] = 0; // Default value for unmatched faces
        }

        // Step 4: Map each new face to a boundary and assign patch names
        for (Foam::label i = 0; i < newFaces.size(); ++i)
        {
            if (newNeighbours[i] == -1)
            {
                const Foam::face &newFace = newFaces[i];
                CGAL::Plane_3<Kernel> newPlane = computePlane(newFace, newPoints);

                // Helper lambda to find a matching plane and return its index
                auto findMatchingPlane = [&](const std::vector<CGAL::Plane_3<Kernel>> &planes) -> Foam::label
                {
                    for (Foam::label j = 0; j < planes.size(); ++j)
                    {
                        const CGAL::Plane_3<Kernel> &plane = planes[j];
                        // Check if planes are parallel (cross product of normal vectors near zero)
                        if (CGAL::cross_product(plane.orthogonal_vector(), newPlane.orthogonal_vector()).squared_length() < tolerance)
                        {
                            // Check if a point from the new face lies on the plane
                            CGAL::Point_3<Kernel> p = Converter::toCGALPoint(newPoints[newFace[0]]);
                            // Compute the plane equation value at point p
                            auto plane_eq_value = plane.a() * p.x() + plane.b() * p.y() + plane.c() * p.z() + plane.d();
                            if (std::abs(CGAL::to_double(plane_eq_value)) < tolerance)
                            {
                                return j; // Return index of matching plane
                            }
                        }
                    }
                    return -1; // No match found
                };

                // Check for a match in polyhedron A
                Foam::label match = findMatchingPlane(planes);
                if (match != -1)
                {
                    // Assign the patch name from the matching face in A
                    newPatches[i] = boundaryPatches[match];
                }
            }
        }
    }

} // namespace Bashyal
