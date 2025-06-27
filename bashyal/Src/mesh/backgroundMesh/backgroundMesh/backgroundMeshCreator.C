#include "backgroundMesh.H"
#include "quickMesh.H"
#include <algorithm>
#include <vector>
#include <set>

using namespace Foam;

namespace Bashyal
{

    void backgroundMesh::developBlocks()
    {
        // Loop through each level of the List structure
        for (auto &blockListLevel1 : backgroundBlocks_)
        {
            for (auto &blockListLevel2 : blockListLevel1)
            {
                for (auto &blockPtr : blockListLevel2)
                {
                    // Dereference the autoPtr to get the actual backgroundBlock
                    backgroundBlock &block = *blockPtr;
                    block.develop();
                }
            }
        }
    }

    void backgroundMesh::auditCells()
    {
        // this->reset();
        Foam::label cellsCount = 0;

        // Loop through each level of the List structure
        for (auto &blockListLevel1 : backgroundBlocks_)
        {
            for (auto &blockListLevel2 : blockListLevel1)
            {
                for (auto &blockPtr : blockListLevel2)
                {
                    // Dereference the autoPtr to get the actual backgroundBlock
                    backgroundBlock &block = *blockPtr;
                    block.globalNCells_ = cellsCount;

                    cellsCount = cellsCount + block.ncells_;
                }
            }
        }
    }

    void backgroundMesh::developMesh()
    {
        this->reset();
        this->auditCells();

        // Loop through each block
        for (Foam::label i = 0; i < dim_[0]; i++)
        {
            for (Foam::label j = 0; j < dim_[1]; j++)
            {
                for (Foam::label k = 0; k < dim_[2]; k++)
                {
                    backgroundBlock &block = *backgroundBlocks_[i][j][k];

                    // Get audited data for this block
                    Foam::pointField updatedPoints;
                    Foam::faceList updatedFaces;
                    Foam::labelList updatedOwners;
                    Foam::labelList updatedNeighbours;
                    Foam::List<int> updatedPatches;

                    getAuditedBlockData(i, j, k, updatedPoints, updatedFaces, updatedOwners, updatedNeighbours, updatedPatches);

                    // Add points and update pointMap_ for unique point indices
                    addPoints(updatedPoints);

                    // Add faces with updated point indices
                    addFaces(block.identity_, updatedPoints, updatedFaces, updatedOwners, updatedNeighbours, updatedPatches);

                    cellCount_ = cellCount_ + block.ncells_;
                }
            }
        }

        Foam::faceList outFaces;
        Foam::labelList outOwners;
        Foam::labelList outNeighbours;
        Foam::List<int> outPatches;

        reorderToUpperTriangularInternal(
            globalFaces_, globalOwners_, globalNeighbours_,
            outFaces, outOwners, outNeighbours);

        globalFaces_ = outFaces;
        globalOwners_ = outOwners;
        globalNeighbours_ = outNeighbours;
    }

    void backgroundMesh::addPoints(const pointField &blockPoints)
    {
        for (const auto &pt : blockPoints)
        {
            Foam::point roundPt = this->roundPoint(pt);
            if (!pointMap_.found(roundPt))
            {
                pointMap_.insert(roundPt, globalPoints_.size());
                globalPoints_.append(roundPt);
            }
        }
    }

    void backgroundMesh::addFaces(const Foam::Vector<int> &identity, const pointField &blockPoints, const faceList &blockFaces, const labelList &owners, const labelList &neighbours, const List<int> &blockPatches)
    {
        int faceCount = 0;
        for (const auto &faceI : blockFaces)
        {
            // Adjust face points according to the globalPoints_
            Foam::face globalFace;
            for (const auto &pt : faceI)
            {
                // Retrieve the point coordinates using pt as an index into blockPoints
                Foam::point pointCoords = blockPoints[pt]; // Using blockPoints to get the actual coordinates

                const Foam::point roundPt = this->roundPoint(pointCoords);

                // Add to pointMap_ if not already present, and retrieve the global index
                label globalPointIdx;
                if (pointMap_.found(roundPt))
                {
                    globalPointIdx = pointMap_[roundPt];
                }
                else
                {
                    // If not found, add the point and assign it a new global index
                    globalPointIdx = globalPoints_.size();
                    globalPoints_.append(roundPt);
                    pointMap_.insert(roundPt, globalPointIdx);
                }

                // Add global point index to the globalFace
                globalFace.append(globalPointIdx);
            }

            int patchType = blockPatches[faceCount];

            if (isFaceMeshBoundary(identity, patchType))
            {
                // Internal face: add to globalFaces_ with owner and neighbor
                boundaryFaces_.append(globalFace);
                boundaryOwners_.append(owners[faceCount]);
                boolBoundaryFaces_.append(true);
                boundaryPatches_.append(patchType); // Store patch type for boundary face
            }
            else
            {
                // Internal face: add to globalFaces_ with owner and neighbor
                globalFaces_.append(globalFace);
                globalOwners_.append(owners[faceCount]);
                boolBoundaryFaces_.append(false);
                globalNeighbours_.append(neighbours[faceCount]);
            }
            faceCount++;
        }
    }

    // Find the common edge between two faces
    std::pair<Foam::label, Foam::label> findCommonEdge(const Foam::face &face1, const Foam::face &face2)
    {
        for (Foam::label i = 0; i < face1.size(); ++i)
        {
            Foam::label p1 = face1[i];
            Foam::label p2 = face1[(i + 1) % face1.size()];
            for (Foam::label j = 0; j < face2.size(); ++j)
            {
                Foam::label q1 = face2[j];
                Foam::label q2 = face2[(j + 1) % face2.size()];
                if ((p1 == q2 && p2 == q1) || (p1 == q1 && p2 == q2))
                {
                    return {p1, p2};
                }
            }
        }
        return {-1, -1}; // No common edge
    }

    // Merge two faces around their common edge
    Foam::face mergeFacesAroundEdge(const Foam::face &face1, const Foam::face &face2,
                                    const std::pair<Foam::label, Foam::label> &commonEdge)
    {
        Foam::label startIdx1 = face1.find(commonEdge.first);
        Foam::label startIdx2 = face2.find(commonEdge.second);
        Foam::face mergedFace;

        // Traverse face1 after the common edge
        for (Foam::label i = (startIdx1 + 1) % face1.size(); i != startIdx1; i = (i + 1) % face1.size())
        {
            mergedFace.append(face1[i]);
        }

        // Traverse face2 after the common edge
        for (Foam::label i = (startIdx2 + 1) % face2.size(); i != startIdx2; i = (i + 1) % face2.size())
        {
            mergedFace.append(face2[i]);
        }

        // Close the face with common edge points
        // mergedFace.append(commonEdge.first);
        // mergedFace.append(commonEdge.second);

        return mergedFace;
    }

    void backgroundMesh::reorderToUpperTriangularInternal(
        const Foam::faceList &faces,
        const Foam::labelList &owners,
        const Foam::labelList &neighbours,
        Foam::faceList &outFaces,
        Foam::labelList &outOwners,
        Foam::labelList &outNeighbours)
    {
        // Create and sort indices by owner, then neighbor
        Foam::labelList indices(faces.size());
        forAll(indices, i) { indices[i] = i; }
        auto compare = [&owners, &neighbours](Foam::label a, Foam::label b) {
            return (owners[a] < owners[b]) || 
                   (owners[a] == owners[b] && neighbours[a] < neighbours[b]);
        };
        std::sort(indices.begin(), indices.end(), compare);
    
        // Temporary storage for merged data
        std::vector<Foam::face> mergedFaces;
        std::vector<Foam::label> mergedOwners;
        std::vector<Foam::label> mergedNeighbours;
    
        // Process faces and merge duplicates
        for (Foam::label idx = 0; idx < indices.size(); ++idx) {
            Foam::label i = indices[idx];
            Foam::label owner = owners[i];
            Foam::label neighbour = neighbours[i];
            Foam::face currentFace = faces[i];
    
            if (idx == 0 || 
                owners[indices[idx - 1]] != owner || 
                neighbours[indices[idx - 1]] != neighbour) {
                // New owner-neighbor pair
                mergedFaces.push_back(currentFace);
                mergedOwners.push_back(owner);
                mergedNeighbours.push_back(neighbour);
            } else {
                // Duplicate pair, merge with the last face
                Foam::face &lastFace = mergedFaces.back();
    
                // Find the common edge
                std::pair<Foam::label, Foam::label> commonEdge = findCommonEdge(lastFace, currentFace);
                if (commonEdge.first != -1 && commonEdge.second != -1) {
                    // Merge faces around the common edge
                    lastFace = mergeFacesAroundEdge(lastFace, currentFace, commonEdge);
                } else {
                    Foam::Warning << "No common edge found between faces with same owner and neighbor." << Foam::endl;
                }
            }
        }
    
        // Assign to output
        outFaces.setSize(mergedFaces.size());
        outOwners.setSize(mergedOwners.size());
        outNeighbours.setSize(mergedNeighbours.size());
        for (unsigned int i = 0; i < mergedFaces.size(); ++i) {
            outFaces[i] = mergedFaces[i];
            outOwners[i] = mergedOwners[i];
            outNeighbours[i] = mergedNeighbours[i];
        }
    }

    void backgroundMesh::reset()
    {
        cellCount_ = 0;
        globalPoints_.clear();
        globalFaces_.clear();
        globalOwners_.clear();
        globalNeighbours_.clear();

        boundaryFaces_.clear();
        boundaryOwners_.clear();
        boolBoundaryFaces_.clear();
        boundaryPatches_.clear();

        pointMap_.clear();
        // faceMap_.clear();
        // faceOwnerMap_.clear();
    }

    point backgroundMesh::roundPoint(point value)
    {
        return point(this->roundToSeven(value[0]), this->roundToSeven(value[1]), this->roundToSeven(value[2]));
    }

    scalar backgroundMesh::roundToSeven(scalar value)
    {
        return std::round(value * std::pow(10.0, 9)) / std::pow(10.0, 9);
    }

}
