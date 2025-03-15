#include "backgroundBlock.H"
#include "foamCGALConverter.H"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h>
#include <map>
#include <vector>

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
        Foam::pointField triPoints;
        Foam::faceList triFaces;
        // Ensure faces are triangulated for CGAL compatibility
        this->triangulateFaces(points_, faces_, triPoints, triFaces);

        // Convert OpenFOAM geometry to CGAL Polyhedron_3 using FoamCGALConverter
        Polyhedron originalPolyhedron;
        try
        {
            originalPolyhedron = Converter::toCGALPolyhedron(triPoints, triFaces);
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
        double volumeTolerance = 1e-8; // Minimum volume threshold
        std::vector<Polyhedron> convexPolyhedra;
        Nef_polyhedron::Volume_const_iterator ci = ++N.volumes_begin();
        for (; ci != N.volumes_end(); ++ci)
        {
            if (ci->mark())
            {
                Polyhedron cp;
                N.convert_inner_shell_to_polyhedron(ci->shells_begin(), cp);
                // Compute volume using a copy
                Polyhedron tempPoly = cp;
                if (!CGAL::is_triangle_mesh(tempPoly))
                {
                    CGAL::Polygon_mesh_processing::triangulate_faces(tempPoly);
                }
                auto volume = CGAL::Polygon_mesh_processing::volume(tempPoly);

                // Convert the exact volume to a double for comparison
                double volumeDouble = CGAL::to_double(volume);

                // Compare the volume with the tolerance
                if (std::abs(volumeDouble) > volumeTolerance)
                {
                    convexPolyhedra.push_back(cp); // cp is still original
                }
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

        Foam::faceList newFaces;
        Foam::labelList newOwners;
        Foam::labelList newNeighbours;
        Foam::List<int> newPatches;
        Foam::label newNboundaries = 0;
        // Foam::faceList boundaryFaces;

        // Foam::HashTable<label, face> faceOwnerMap_; // For tracking boundary faces
        Foam::HashTable<Foam::label, Foam::face> facePositionMap; // For tracking face Indices

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
                for (unsigned int i = 0; i < pointIndices.size(); ++i)
                {
                    f[i] = pointIndices[i];
                }

                Foam::face faceCopy = Foam::face(f);

                // Sort the face points to ensure unique representation (important for shared faces)
                std::sort(faceCopy.begin(), faceCopy.end());

                if (facePositionMap.found(faceCopy))
                {
                    Foam::label faceLocation = facePositionMap[faceCopy];
                    newNeighbours[faceLocation] = cellI;
                    newPatches[faceLocation] = 0;
                    newNboundaries--;
                    // Face with same orientation exists; add current cell and face to its list
                    // faceMap[faceCopy].push_back(std::make_pair(cellI, f));
                }
                else
                {
                    // Neither f nor its reverse exists; insert new entry with f as key
                    newFaces.append(f);
                    newOwners.append(cellI);
                    newNeighbours.append(-1);
                    facePositionMap.insert(faceCopy, newFaces.size() - 1);
                    newPatches.append(-99);
                    newNboundaries++;
                    // faceMap.insert(f, {std::make_pair(cellI, f)});
                }
                // Note: The original line `FaceKey key(pointIndices); faceMap[key].push_back(...)`
                // is replaced by the else clause above due to key type mismatch
            }
        }

        mapNewFacesToBoundaries(
            allPoints, newFaces, newNeighbours,
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
        Foam::label loopCount = 0;
        for (Foam::label i = 0; i < faces_.size(); ++i)
        {
            if (neighbours_[i] == -1)
            {
                planes[loopCount] = computePlane(faces_[i], points_);
                boundaryPatches[loopCount] = patches_[i];
                loopCount++;
            }
        }

        // Step 4: Map each new face to a boundary and assign patch names
        for (Foam::label i = 0; i < newFaces.size(); ++i)
        {
            if (newPatches[i] == -99)
            {
                const Foam::face &newFace = newFaces[i];
                CGAL::Plane_3<Kernel> newPlane = computePlane(newFace, newPoints);

                // Helper lambda to find a matching plane and return its index
                auto findMatchingPlane = [&](const std::vector<CGAL::Plane_3<Kernel>> &planes) -> Foam::label
                {
                    for (unsigned int j = 0; j < planes.size(); ++j)
                    {
                        const CGAL::Plane_3<Kernel> &plane = planes[j];
                        // Check if planes are parallel (cross product of normal vectors near zero)
                        if (CGAL::cross_product(plane.orthogonal_vector(), newPlane.orthogonal_vector()).squared_length() < tolerance)
                        {
                            bool allPointsOnPlane = true;
                            // Loop through all points in newFace
                            for (Foam::label pointIndex : newFace)
                            {
                                CGAL::Point_3<Kernel> p = Converter::toCGALPoint(newPoints[pointIndex]);
                                // Compute the plane equation value at point p
                                auto plane_eq_value = plane.a() * p.x() + plane.b() * p.y() + plane.c() * p.z() + plane.d();
                                // Check if the point lies on the plane (within tolerance)
                                if (std::abs(CGAL::to_double(plane_eq_value)) >= tolerance)
                                {
                                    allPointsOnPlane = false;
                                    break; // Exit early if any point does not lie on the plane
                                }
                            }
                            // If all points lie on the plane, return the plane index
                            if (allPointsOnPlane)
                            {
                                return j;
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
