#include "backgroundMesh.H"
#include "quickMesh.H"

using namespace Foam;

namespace Bashyal
{

    void backgroundMesh::auditMesh()
    {
        this->reset();
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
                    block.decomposeToConvex();

                    cellsCount = cellsCount + block.ncells_;

                    // auditBlockInternalBoundaryFaces(block);
                }
            }
        }

        // Write polyMesh here using globalPoints_, globalFaces_, globalOwners_, globalNeighbours_, and allFaces_
        // writePolyMesh();
    }

    void backgroundMesh::developMesh()
    {
        this->reset();

        // Loop through each level of the List structure
        for (auto &blockListLevel1 : backgroundBlocks_)
        {
            for (auto &blockListLevel2 : blockListLevel1)
            {
                for (auto &blockPtr : blockListLevel2)
                {
                    // Dereference the autoPtr to get the actual backgroundBlock
                    backgroundBlock &block = *blockPtr;
                    block.decomposeToConvex();

                    // Add points and update pointMap_ for unique point indices
                    addPoints(block.getPoints());

                    // Add faces with updated point indices
                    addFaces(block.globalNCells_, block.identity_, block.getPoints(), block.getFaces(), block.getOwners(), block.getNeighbours(), block.getPatches());

                    cellCount_ = cellCount_ + block.ncells_;
                }
            }
        }

        // Write polyMesh here using globalPoints_, globalFaces_, globalOwners_, globalNeighbours_, and allFaces_
        // writePolyMesh();
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

    void backgroundMesh::addFaces(Foam::label globalNCells, const Foam::Vector<int> &identity, const pointField &blockPoints, const faceList &blockFaces, const labelList &owners, const labelList &neighbours, const List<int> &blockPatches)
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
                boundaryOwners_.append(globalNCells + owners[faceCount]);
                boolBoundaryFaces_.append(true);
                boundaryPatches_.append(patchType); // Store patch type for boundary face
            }
            else
            {
                // Internal face: add to globalFaces_ with owner and neighbor
                globalFaces_.append(globalFace);
                globalOwners_.append(globalNCells + owners[faceCount]);
                globalNeighbours_.append(globalNCells + neighbours[faceCount]);
                boolBoundaryFaces_.append(false);
            }
            faceCount++;
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
        return std::round(value * std::pow(10.0, 7)) / std::pow(10.0, 7);
    }

}
