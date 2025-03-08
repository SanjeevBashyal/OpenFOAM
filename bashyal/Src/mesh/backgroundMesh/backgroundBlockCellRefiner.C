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
    typedef FoamCGALConverter<Kernel> Converter;

    // Custom modifier class for building the polyhedron
    template <class HDS>
    class BuildPolyhedron : public CGAL::Modifier_base<HDS>
    {
    private:
        const std::vector<Kernel::Point_3> &points;
        const Foam::faceList &faces;

    public:
        BuildPolyhedron(const std::vector<Kernel::Point_3> &pts, const Foam::faceList &fcs)
            : points(pts), faces(fcs) {}

        void operator()(HDS &hds) override
        {
            CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);
            builder.begin_surface(points.size(), faces.size());

            // Add vertices
            for (const auto &point : points)
            {
                builder.add_vertex(point);
            }

            // Add facets
            for (const auto &face : faces)
            {
                builder.begin_facet();
                for (Foam::label idx : face)
                {
                    builder.add_vertex_to_facet(idx);
                }
                builder.end_facet();
            }

            builder.end_surface();
        }
    };

    void backgroundBlock::decomposeToConvex()
    {
        // Step 1: Map OpenFOAM points to their indices
        std::map<Kernel::Point_3, Foam::label> pointToIndex;
        forAll(points_, i)
        {
            pointToIndex[Converter::toCGALPoint(points_[i])] = i;
        }

        // Step 2: Convert points_ to CGAL points
        std::vector<Kernel::Point_3> cgalPoints;
        for (const Foam::point &p : points_)
        {
            cgalPoints.push_back(Converter::toCGALPoint(p));
        }

        this->triangulateFaces();

        // Step 3: Build Polyhedron_3 using the custom modifier
        Polyhedron P;
        BuildPolyhedron<Polyhedron::HalfedgeDS> builder(cgalPoints, triFaces_);
        P.delegate(builder);

        if (!P.is_valid())
        {
            FatalErrorInFunction << "Polyhedron construction failed" << exit(Foam::FatalError);
        }

        // Step 4: Convert to Nef_polyhedron_3 and decompose
        Nef_polyhedron N(P);
        CGAL::convex_decomposition_3(N);

        // Step 5: Extract convex parts from N
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
        Foam::HashTable<Foam::label, Foam::point> pointMap;
        Foam::pointField points; // Stores the vertices of the cube

        for (Foam::label cellI = 0; cellI < nCells; ++cellI)
        {
            const Polyhedron &cp = convexPolyhedra[cellI];
            for (Polyhedron::Vertex_const_iterator vertex = cp.vertices_begin(); vertex != cp.vertices_end(); ++vertex)
            {
                Foam::point pt = Converter::toFoamPoint(vertex->point()); // Fixed: use vertex->point()
                if (!pointMap.found(pt))
                {
                    pointMap.insert(pt, points.size());
                    points.append(pt); // No rounding here to match original intent
                }
            }
        }

        // Foam::HashTable<Foam::label, Foam::face> faceOwnerMap_;    // For tracking boundary faces
        // Foam::HashTable<Foam::label, Foam::face> facePositionMap_; // For tracking face Indices

        // Step 6: Collect faces from convex polyhedra
        struct FaceKey : public std::vector<Foam::label>
        {
            FaceKey(const std::vector<Foam::label> &v) : std::vector<Foam::label>(v) { std::sort(begin(), end()); }
        };
        std::map<FaceKey, std::vector<std::pair<Foam::label, Foam::face>>> faceMap;

        for (Foam::label cellI = 0; cellI < nCells; ++cellI)
        {
            const Polyhedron &cp = convexPolyhedra[cellI];
            for (Polyhedron::Facet_const_iterator facet = cp.facets_begin(); facet != cp.facets_end(); ++facet)
            {
                std::vector<Foam::label> pointIndices;
                auto circ = facet->facet_begin();
                do
                {
                    Kernel::Point_3 p = circ->vertex()->point();
                    pointIndices.push_back(pointMap[Converter::toFoamPoint(p)]);
                    ++circ;
                } while (circ != facet->facet_begin());

                // Convert std::vector<Foam::label> to Foam::labelList
                Foam::labelList faceLabels(pointIndices.size());
                for (Foam::label i = 0; i < pointIndices.size(); ++i)
                {
                    faceLabels[i] = pointIndices[i];
                }

                // Create the face using Foam::labelList
                Foam::face f(faceLabels);
                FaceKey key(pointIndices);
                faceMap[key].push_back(std::make_pair(cellI, f));
            }
        }

        // Step 7: Map original faces to their patch and stringPtr
        std::map<FaceKey, Foam::label> originalFaceMap;
        forAll(faces_, faceI)
        {
            const Foam::face &f = faces_[faceI];
            std::vector<Foam::label> pointIndices(f.begin(), f.end());
            FaceKey key(pointIndices);
            originalFaceMap[key] = faceI;
        }

        // Step 8: Build new face lists
        Foam::faceList newFaces;
        Foam::labelList newOwners;
        Foam::labelList newNeighbours;
        Foam::wordList newPatches;
        Foam::wordList newStringPtrs;

        for (const auto &entry : faceMap)
        {
            const FaceKey &key = entry.first;
            const std::vector<std::pair<Foam::label, Foam::face>> &instances = entry.second;
            if (instances.size() == 1)
            {
                // Boundary face
                Foam::label cellI = instances[0].first;
                const Foam::face &f = instances[0].second;
                Foam::label faceI = newFaces.size();
                newFaces.append(f);
                newOwners.append(cellI);
                newNeighbours.append(-1);
                auto it = originalFaceMap.find(key);
                newPatches.append(Foam::word::null);
                newStringPtrs.append(Foam::word::null);
                // if (it != originalFaceMap.end())
                // {
                //     Foam::label origFaceI = it->second;
                //     newPatches.append(patches_[origFaceI]);
                //     newStringPtrs.append(stringPtrs_[origFaceI]);
                // }
                // else
                // {
                //     FatalErrorInFunction << "Boundary face not in original faces" << exit(Foam::FatalError);
                // }
            }
            else if (instances.size() == 2)
            {
                // Internal face
                Foam::label cellI = instances[0].first;
                Foam::label cellB = instances[1].first;
                Foam::label owner = std::min(cellI, cellB);
                Foam::label neighbour = std::max(cellI, cellB);
                const Foam::face &fOwner = (cellI == owner) ? instances[0].second : instances[1].second;
                Foam::label faceI = newFaces.size();
                newFaces.append(fOwner);
                newOwners.append(owner);
                newNeighbours.append(neighbour);
                newPatches.append(Foam::word::null);
                newStringPtrs.append(Foam::word::null);
            }
            else
            {
                FatalErrorInFunction << "Face shared by >2 cells" << exit(Foam::FatalError);
            }
        }

        // Step 9: Update class attributes
        points_ = points;
        faces_ = newFaces;
        owners_ = newOwners;
        neighbours_ = newNeighbours;
        patches_ = newPatches;
        stringPtrs_ = newStringPtrs;
        ncells_ = convexPolyhedra.size();
        edited_ = true;
    }

} // namespace Bashyal