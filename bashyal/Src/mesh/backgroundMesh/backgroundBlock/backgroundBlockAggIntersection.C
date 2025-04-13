#include "backgroundBlock.H"
#include "foamCGALConverter.H" // Assumed converter header
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h>

namespace Bashyal
{
    // Define CGAL types with an exact kernel for robustness
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
    typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
    typedef FoamCGALConverter<Kernel> Converter; // Adjust namespace as needed

    // Subtracts a single aggregate from the current nef_
    void backgroundBlock::subtractAggregate(
        const aggregate &agg, // Pass aggregate object by const reference
        int aggregateIdentifier)
    {
        if (dead_)
            return;

        // Access the Nef polyhedron from the aggregate object
        const CGAL::Nef_polyhedron_3<Kernel_t> &aggNef = agg.nef_;

        if (aggNef.is_empty())
        {
            // Optional Warning:
            // Foam::WarningInFunc() << "Skipping subtraction of aggregate with empty Nef from block " << blockID_ << Foam::endl;
            return;
        }

        try
        {
            // Perform subtraction: nef_ = nef_ - aggNef
            nef_ = nef_.difference(aggNef);

            // Check result
            if (nef_.is_empty())
            {
                // Foam::Info << "Block " << blockID_ << " became empty after subtracting aggregate ID " << aggregateIdentifier << "." << Foam::endl;
                dead_ = true;
                points_.clear();
                faces_.clear();
                patches_.clear();
                owners_.clear();
                neighbours_.clear();
                ncells_ = 0;
                edited_ = true; // Mark as edited even if now dead
                return;
            }
            // Don't check for non-simple here, difference can create them

            // Convert the result back to Foam geometry
            // Pass the aggregateIdentifier for patch mapping
            bool success = updateFoamGeometryFromNef("Aggregate Subtraction", aggregateIdentifier);
            if (!success)
            {
                Foam::WarningInFunc() << "Failed to update Foam geometry after subtracting aggregate "
                                      << aggregateIdentifier << " from block " << blockID_ << ". Marking as dead." << Foam::endl;
                dead_ = true;
                points_.clear();
                faces_.clear();
                patches_.clear();
                owners_.clear();
                neighbours_.clear();
                ncells_ = 0;
            }
            edited_ = true; // Mark as edited because an aggregate interaction occurred
        }
        catch (const std::exception &e)
        {
            Foam::WarningInFunc() << "Error subtracting aggregate " << aggregateIdentifier
                                  << " for block " << blockID_ << ": " << e.what() << ". Keeping previous block state." << Foam::endl;
            // Don't change nef_ or Foam geometry on error
        }
        catch (...)
        {
            Foam::WarningInFunc() << "Unknown error subtracting aggregate " << aggregateIdentifier
                                  << " for block " << blockID_ << ". Keeping previous block state." << Foam::endl;
            // Don't change nef_ or Foam geometry on error
        }
    }

    void backgroundBlock::intersectClosedSurfaceCGAL(
        const Foam::faceList &faces,
        const Foam::pointField &points,
        unsigned int identifier)
    {
        Foam::pointField triPointsA;
        Foam::faceList triFacesA;
        // Ensure faces are triangulated for CGAL compatibility
        this->triangulateFaces(points_, faces_, triPointsA, triFacesA);

        Foam::pointField triPointsB;
        Foam::faceList triFacesB;
        // Ensure faces are triangulated for CGAL compatibility
        this->triangulateFaces(points, faces, triPointsB, triFacesB);

        // Step 1: Convert block's geometry (A) to CGAL Polyhedron
        Polyhedron A;
        try
        {
            A = Converter::toCGALPolyhedron(triPointsA, triFacesA);
            if (!A.is_valid())
            {
                Foam::Warning << "Block polyhedron A is invalid." << Foam::endl;
                dead_ = true;
                ncells_ = 0;
                return;
            }
        }
        catch (const std::exception &e)
        {
            Foam::Warning << "Failed to convert block to Polyhedron: " << e.what() << Foam::endl;
            dead_ = true;
            ncells_ = 0;
            return;
        }

        // Step 2: Convert input geometry (B) to CGAL Polyhedron
        Polyhedron B;
        try
        {
            B = Converter::toCGALPolyhedron(triPointsB, triFacesB);
            if (!B.is_valid())
            {
                Foam::Warning << "Input polyhedron B is invalid." << Foam::endl;
                dead_ = true;
                ncells_ = 0;
                return;
            }
        }
        catch (const std::exception &e)
        {
            Foam::Warning << "Failed to convert input to Polyhedron: " << e.what() << Foam::endl;
            dead_ = true;
            ncells_ = 0;
            return;
        }

        // Step 3: Convert to Nef_polyhedron_3 for subtraction
        Nef_polyhedron nefA(A);
        Nef_polyhedron nefB(B);

        // Step 4: Compute difference A - B
        Nef_polyhedron N = nefA - nefB;

        // Step 5: Handle the result
        if (N.number_of_volumes() <= 1)
        { // Only the exterior volume, i.e., result is empty
            points_.clear();
            faces_.clear();
            patches_.clear();
            dead_ = true;
            ncells_ = 0;
            edited_ = true;
            Foam::Info << "Subtraction resulted in an empty polyhedron. Block marked as dead." << Foam::endl;
            return;
        }

        CGAL::convex_decomposition_3(N);

        // Extract convex polyhedra from the decomposition
        double volumeTolerance = 1e-12; // Minimum volume threshold
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

        Foam::List<int> patches;
        patches.setSize(faces.size());
        for (int i = 0; i < faces.size(); ++i)
        {
            patches[i] = identifier;
        }

        mapNewFacesToBoundaries(
            allPoints, newFaces, newNeighbours, // New polyhedron data
            points_, faces_, patches_,          // Original polyhedron data
            points, faces, patches,             // Input polyhedron data
            newPatches                          // Output patches
        );

        // Foam::faceList outFaces;
        // Foam::labelList outOwners;
        // Foam::labelList outNeighbours;
        // Foam::List<int> outPatches;

        // reorderToUpperTriangularNeighboursOnly(
        //     newFaces, newOwners, newNeighbours, newPatches,
        //     outFaces, outOwners, outNeighbours, outPatches);

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
        const Foam::pointField &pointsA,
        const Foam::faceList &facesA,
        const Foam::List<int> &patchesA,
        const Foam::pointField &pointsB,
        const Foam::faceList &facesB,
        const Foam::List<int> &patchesB,
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

        // Step 1: Compute planes for boundary faces of polyhedron A (using neighbours_)
        std::vector<CGAL::Plane_3<Kernel>> planesA;
        std::vector<int> boundaryPatchesA;
        for (Foam::label i = 0; i < facesA.size(); ++i)
        {
            if (neighbours_[i] == -1) // Only boundary faces
            {
                planesA.push_back(computePlane(facesA[i], pointsA));
                boundaryPatchesA.push_back(patchesA[i]);
            }
        }

        // Step 2: Compute planes for all faces of polyhedron B (all are boundary faces)
        std::vector<CGAL::Plane_3<Kernel>> planesB(facesB.size());
        for (Foam::label i = 0; i < facesB.size(); ++i)
        {
            planesB[i] = computePlane(facesB[i], pointsB);
        }

        // Step 3: Resize newPatches and initialize with 0
        newPatches.setSize(newFaces.size(), 0);

        // Step 4: Map each new boundary face to a boundary and assign patches
        for (Foam::label i = 0; i < newFaces.size(); ++i)
        {
            if (newPatches[i] == -99) // Only process boundary faces
            {
                const Foam::face &newFace = newFaces[i];
                CGAL::Plane_3<Kernel> newPlane = computePlane(newFace, newPoints);

                // Helper lambda to find a matching plane
                auto findMatchingPlane = [&](const std::vector<CGAL::Plane_3<Kernel>> &planes) -> Foam::label
                {
                    for (unsigned int j = 0; j < planes.size(); ++j)
                    {
                        const CGAL::Plane_3<Kernel> &plane = planes[j];
                        if (CGAL::cross_product(plane.orthogonal_vector(), newPlane.orthogonal_vector()).squared_length() < tolerance)
                        {
                            bool allPointsOnPlane = true;
                            for (Foam::label pointIndex : newFace)
                            {
                                CGAL::Point_3<Kernel> p = Converter::toCGALPoint(newPoints[pointIndex]);
                                auto plane_eq_value = plane.a() * p.x() + plane.b() * p.y() + plane.c() * p.z() + plane.d();
                                if (std::abs(CGAL::to_double(plane_eq_value)) >= tolerance)
                                {
                                    allPointsOnPlane = false;
                                    break;
                                }
                            }
                            if (allPointsOnPlane)
                            {
                                return j;
                            }
                        }
                    }
                    return -1;
                };

                // Check for a match in polyhedron A
                Foam::label matchA = findMatchingPlane(planesA);
                if (matchA != -1)
                {
                    newPatches[i] = boundaryPatchesA[matchA];
                }
                else
                {
                    // Check for a match in polyhedron B
                    Foam::label matchB = findMatchingPlane(planesB);
                    if (matchB != -1)
                    {
                        newPatches[i] = patchesB[matchB];
                    }
                    // Unmatched boundary faces remain 0
                }
            }
            // Internal faces remain 0
        }
    }

    void backgroundBlock::reorderToUpperTriangularNeighboursOnly(
        const Foam::faceList &faces,
        const Foam::labelList &owners,
        const Foam::labelList &neighbours,
        const Foam::List<int> &patches,
        Foam::faceList &outFaces,
        Foam::labelList &outOwners,
        Foam::labelList &outNeighbours,
        Foam::List<int> &outPatches)
    {
        // Set the size of output lists to match input
        outFaces.setSize(faces.size());
        outOwners.setSize(owners.size());
        outNeighbours.setSize(neighbours.size());
        outPatches.setSize(patches.size());

        // Initialize position tracker for output lists
        // Foam::label outPos = 0;

        // Variable to store the current owner handle
        Foam::label currentOwner = -1;

        Foam::label newOwnerStartIndex = 0;

        // Start looping through all faces
        for (Foam::label faceI = 0; faceI < faces.size(); ++faceI)
        {
            Foam::label faceOwner = owners[faceI];

            // Check if we’ve encountered a new owner
            if (faceOwner != currentOwner)
            {
                // Update the current owner
                currentOwner = faceOwner;

                // Place the first face of this owner directly
                outFaces[faceI] = faces[faceI];
                outOwners[faceI] = owners[faceI];
                outNeighbours[faceI] = neighbours[faceI];
                outPatches[faceI] = patches[faceI];
                newOwnerStartIndex = faceI;
            }
            else
            {
                // Compare with previous neighbours for the same owner
                Foam::label currentNeighbour = neighbours[faceI];
                Foam::label insertPos = faceI;

                // Find the correct insertion position to maintain ascending order
                for (Foam::label j = newOwnerStartIndex; j <= faceI; j++)
                {
                    if (currentNeighbour < outNeighbours[j])
                    {
                        insertPos = j;
                        break;
                    }
                    else
                    {
                        continue; // Stop when we find a smaller or equal neighbour
                    }
                }

                // Shift elements to the right to make space for insertion
                for (Foam::label k = faceI; k > insertPos; --k)
                {
                    outFaces[k] = outFaces[k - 1];
                    outOwners[k] = outOwners[k - 1];
                    outNeighbours[k] = outNeighbours[k - 1];
                    outPatches[k] = outPatches[k - 1];
                }

                // Insert the current face at the correct position
                outFaces[insertPos] = faces[faceI];
                outOwners[insertPos] = owners[faceI];
                outNeighbours[insertPos] = currentNeighbour;
                outPatches[insertPos] = patches[faceI];
            }
        }
    }

    //     void backgroundBlock::reorderToUpperTriangular(
    //         const Foam::faceList& faces,
    //         const Foam::labelList& owners,
    //         const Foam::labelList& neighbours,
    //         const Foam::List<int>& patches,
    //         Foam::faceList& outFaces,
    //         Foam::labelList& outOwners,
    //         Foam::labelList& outNeighbours,
    //         Foam::List<int>& outPatches
    //     )
    //     {
    //         // Step 1: Ensure owner < neighbour for internal faces and reverse face if needed
    //         Foam::faceList modFaces = faces;
    //         Foam::labelList modOwners = owners;
    //         Foam::labelList modNeighbours = neighbours;
    //         Foam::List<int> modPatches = patches;

    //         for (Foam::label faceI = 0; faceI < faces.size(); ++faceI)
    //         {
    //             if (modNeighbours[faceI] != -1)  // Internal face
    //             {
    //                 if (modOwners[faceI] > modNeighbours[faceI])
    //                 {
    //                     // Swap owner and neighbour
    //                     std::swap(modOwners[faceI], modNeighbours[faceI]);
    //                     // Reverse face orientation
    //                     std::reverse(modFaces[faceI].begin(), modFaces[faceI].end());
    //                 }
    //             }
    //         }

    //         // Step 2: Create a list of indices to sort
    //         std::vector<Foam::label> indices(faces.size());
    //         for (Foam::label i = 0; i < faces.size(); ++i)
    //         {
    //             indices[i] = i;
    //         }

    //         // Step 3: Sort indices based on neighbour for each owner group
    //         Foam::label currentOwner = -1;
    //         Foam::label groupStart = 0;

    //         for (Foam::label i = 0; i <= indices.size(); ++i)
    //         {
    //             if (i == indices.size() || modOwners[indices[i]] != currentOwner)
    //             {
    //                 // Sort the previous group by neighbour
    //                 if (i > groupStart)
    //                 {
    //                     std::sort(indices.begin() + groupStart, indices.begin() + i,
    //                         [&](Foam::label a, Foam::label b)
    //                         {
    //                             return modNeighbours[a] < modNeighbours[b];
    //                         });
    //                 }
    //                 // Start a new group
    //                 if (i < indices.size())
    //                 {
    //                     currentOwner = modOwners[indices[i]];
    //                     groupStart = i;
    //                 }
    //             }
    //         }

    //         // Step 4: Reorder all lists based on sorted indices
    //         outFaces.setSize(faces.size());
    //         outOwners.setSize(faces.size());
    //         outNeighbours.setSize(faces.size());
    //         outPatches.setSize(faces.size());

    //         for (Foam::label i = 0; i < indices.size(); ++i)
    //         {
    //             Foam::label originalIndex = indices[i];
    //             outFaces[i] = modFaces[originalIndex];
    //             outOwners[i] = modOwners[originalIndex];
    //             outNeighbours[i] = modNeighbours[originalIndex];
    //             outPatches[i] = modPatches[originalIndex];
    //         }
    //     }
    // }
}
