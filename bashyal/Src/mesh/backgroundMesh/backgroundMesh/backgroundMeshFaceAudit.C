#include "backgroundMesh.H"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <list>
#include <map>

using namespace Foam;

namespace Bashyal
{
    // Define CGAL types
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef Kernel::Point_2 Point_2;
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;

    void backgroundMesh::getAuditedBlockData(
        Foam::label i, Foam::label j, Foam::label k,
        Foam::pointField &updatedPoints,
        Foam::faceList &updatedFaces,
        Foam::labelList &updatedOwners,
        Foam::labelList &updatedNeighbours,
        Foam::List<int> &updatedPatches)
    {
        // Access the block without modifying it
        backgroundBlock &block = *backgroundBlocks_[i][j][k];
        const Foam::Vector<int> &identity = block.identity_;

        // Initialize updatedPoints with the block's current points
        updatedPoints = block.points_;

        // Temporary lists for updated data
        Foam::faceList newFaces;
        Foam::labelList newOwners;
        Foam::labelList newNeighbours;
        Foam::List<int> newPatches;

        // Process interface faces (patchType 1, 2, 3)
        for (int patchType = 3; patchType >= 1; --patchType)
        {
            if (isFaceMeshBoundary(identity, patchType))
                continue; // Skip if this is a mesh boundary

            // Collect faces with the current patchType from the block
            Foam::faceList blockFaces;
            Foam::labelList blockFaceOwners;
            forAll(block.faces_, faceI)
            {
                if (block.patches_[faceI] == patchType)
                {
                    blockFaces.append(block.faces_[faceI]);
                    blockFaceOwners.append(block.owners_[faceI]); // Local owner index
                }
            }

            if (blockFaces.empty())
                continue;

            // Identify the neighboring block
            backgroundBlock *neighborBlock = nullptr;
            if (patchType == 1 && i < dim_[0] - 1)
                neighborBlock = backgroundBlocks_[i + 1][j][k].get();
            else if (patchType == 2 && j < dim_[1] - 1)
                neighborBlock = backgroundBlocks_[i][j + 1][k].get();
            else if (patchType == 3 && k < dim_[2] - 1)
                neighborBlock = backgroundBlocks_[i][j][k + 1].get();

            if (neighborBlock)
            {
                // Collect neighbor faces with the corresponding negative patchType
                Foam::faceList neighborFaces;
                Foam::labelList neighborFaceOwners;
                forAll(neighborBlock->faces_, nFaceI)
                {
                    if (neighborBlock->patches_[nFaceI] == -patchType)
                    {
                        neighborFaces.append(neighborBlock->faces_[nFaceI].reverseFace());
                        neighborFaceOwners.append(neighborBlock->globalNCells_ + neighborBlock->owners_[nFaceI]); // Global neighbor index
                    }
                }

                // Perform face auditing (intersections)
                Foam::faceList intersectedFaces;
                Foam::labelList intersectedOwners;     // Local owners
                Foam::labelList intersectedNeighbours; // Global neighbours
                Foam::List<int> intersectedPatches;

                auditSinglePatchFaces(
                    updatedPoints, // Passed by reference to add new points
                    blockFaces,
                    blockFaceOwners,
                    neighborBlock->points_,
                    neighborFaces,
                    neighborFaceOwners,
                    patchType,
                    intersectedFaces,
                    intersectedOwners,
                    intersectedNeighbours,
                    intersectedPatches);

                // Append intersected faces to the new lists
                forAll(intersectedFaces, idx)
                {
                    newFaces.append(intersectedFaces[idx]);
                    newOwners.append(block.globalNCells_ + intersectedOwners[idx]); // Local owner
                    newNeighbours.append(intersectedNeighbours[idx]);               // Global neighbour
                    newPatches.append(intersectedPatches[idx]);
                }
            }
        }

        // Copy faces that are not patchType 1, 2, or 3
        forAll(block.faces_, faceI)
        {
            int patchType = block.patches_[faceI];
            if (patchType == 0)
            {
                newFaces.append(block.faces_[faceI]);
                newOwners.append(block.globalNCells_ + block.owners_[faceI]);
                newNeighbours.append(block.globalNCells_ + block.neighbours_[faceI]);
                newPatches.append(patchType);
            }
            else if (isFaceMeshBoundary(identity, patchType))
            {
                newFaces.append(block.faces_[faceI]);
                newOwners.append(block.globalNCells_ + block.owners_[faceI]);
                newNeighbours.append(block.neighbours_[faceI]);
                newPatches.append(patchType);
            }
        }

        // Assign the updated data to output parameters
        updatedFaces = newFaces;
        updatedOwners = newOwners;
        updatedNeighbours = newNeighbours;
        updatedPatches = newPatches;
    }

    void backgroundMesh::auditSinglePatchFaces(
        Foam::pointField &blockPoints,
        const Foam::faceList &blockFaces,
        const Foam::labelList &blockFaceOwners,
        const Foam::pointField &neighborPoints,
        const Foam::faceList &neighborFaces,
        const Foam::labelList &neighborFaceOwners,
        int patchType,
        Foam::faceList &intersectedFaces,
        Foam::labelList &intersectedOwners,
        Foam::labelList &intersectedNeighbours,
        Foam::List<int> &intersectedPatches)
    {
        // Define area tolerance
        const double areaTolerance = 1e-12;

        // Lambda to project 3D points to 2D based on patchType
        auto projectTo2D = [&](const Foam::point &p) -> Point_2
        {
            if (patchType == 1)
                return Point_2(p.y(), p.z()); // YZ-plane
            if (patchType == 2)
                return Point_2(p.x(), p.z()); // XZ-plane
            return Point_2(p.x(), p.y());     // XY-plane (default)
        };

        // Convert block faces to Polygon_2
        std::vector<Polygon_2> blockPwhs;
        std::vector<std::vector<Foam::label>> blockPointIndices;
        for (Foam::label i = 0; i < blockFaces.size(); ++i)
        {
            const auto &face = blockFaces[i];
            Polygon_2 poly;
            std::vector<Foam::label> indices;
            for (const auto &idx : face)
            {
                poly.push_back(projectTo2D(blockPoints[idx]));
                indices.push_back(idx);
            }
            blockPwhs.push_back(poly);
            blockPointIndices.push_back(indices);
        }

        // Convert neighbor faces to Polygon_2
        std::vector<Polygon_2> neighborPwhs;
        std::vector<std::vector<Foam::label>> neighborPointIndices;
        for (Foam::label i = 0; i < neighborFaces.size(); ++i)
        {
            const auto &face = neighborFaces[i];
            Polygon_2 polyN;
            std::vector<Foam::label> indices;
            for (const auto &idx : face)
            {
                polyN.push_back(projectTo2D(neighborPoints[idx]));
                indices.push_back(idx);
            }
            neighborPwhs.push_back(polyN);
            neighborPointIndices.push_back(indices);
        }

        // Compute intersections
        for (size_t i = 0; i < blockPwhs.size(); ++i)
        {
            auto &blockPwh = blockPwhs[i];
            Foam::label blockOwner = blockFaceOwners[i];
            const auto &blockIndices = blockPointIndices[i];

            bool originalReverseFlag = false;

            // Ensure polygons are counter-clockwise for intersection
            if (blockPwh.orientation() != CGAL::COUNTERCLOCKWISE)
            {
                blockPwh.reverse_orientation();
                originalReverseFlag = true;
            }

            for (size_t j = 0; j < neighborPwhs.size(); ++j)
            {
                auto &neighborPwh = neighborPwhs[j];
                Foam::label neighborOwner = neighborFaceOwners[j];

                if (neighborPwh.orientation() != CGAL::COUNTERCLOCKWISE)
                {
                    neighborPwh.reverse_orientation();
                }

                // Perform intersection
                std::list<Polygon_with_holes_2> intersection_result;
                CGAL::intersection(blockPwh, neighborPwh, std::back_inserter(intersection_result));

                // Process each intersection result
                for (auto &pwh : intersection_result)
                {
                    auto &outer = pwh.outer_boundary();
                    if (originalReverseFlag)
                    {
                        outer.reverse_orientation();
                    }

                    // Compute the area of the outer boundary
                    double area = std::abs(CGAL::to_double(outer.area()));

                    // Proceed only if area exceeds tolerance
                    if (area > areaTolerance)
                    {
                        Foam::face intersectedFace;
                        for (const auto &pt : outer)
                        {
                            // Convert exact coordinates to double
                            double x_val = CGAL::to_double(pt.x());
                            double y_val = CGAL::to_double(pt.y());
                            double x, y, z;
                            if (patchType == 1)
                            {                                         // YZ-plane
                                x = blockPoints[blockIndices[0]].x(); // Fixed x
                                y = x_val;
                                z = y_val;
                            }
                            else if (patchType == 2)
                            {                                         // XZ-plane
                                y = blockPoints[blockIndices[0]].y(); // Fixed y
                                x = x_val;
                                z = y_val;
                            }
                            else
                            {                                         // XY-plane
                                z = blockPoints[blockIndices[0]].z(); // Fixed z
                                x = x_val;
                                y = y_val;
                            }
                            Foam::point p3d(x, y, z);
                            Foam::label ptIdx = findOrAddPoint(blockPoints, p3d);
                            intersectedFace.push_back(ptIdx);
                        }
                        if (intersectedFace.size() >= 3)
                        {
                            intersectedFaces.push_back(intersectedFace);
                            intersectedOwners.push_back(blockOwner);
                            intersectedNeighbours.push_back(neighborOwner);
                            intersectedPatches.push_back(patchType);
                        }
                    }
                    // Note: Holes in pwh are not processed, matching original behavior
                }
            }
        }
    }

    Foam::label backgroundMesh::findOrAddPoint(Foam::pointField &blockPoints, const Foam::point &p)
    {
        const scalar tol = 1e-6; // Tolerance for point comparison
        for (Foam::label i = 0; i < blockPoints.size(); ++i)
        {
            if (mag(blockPoints[i] - p) < tol)
            {
                return i; // Return index if point exists
            }
        }
        // Add new point and return its index
        blockPoints.append(p);
        return blockPoints.size() - 1;
    }
}
