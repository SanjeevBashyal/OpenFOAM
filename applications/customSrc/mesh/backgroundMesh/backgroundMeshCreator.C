#include "backgroundMesh.H"
#include "quickMesh.H"

using namespace Foam;

namespace Bashyal
{
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

                    // Add points and update pointMap_ for unique point indices
                    addPoints(block.getPoints());

                    // Add faces with updated point indices
                    addFaces(block.getPoints(), block.getFaces(), block.getPatches(), block.getStringPtrs());

                    cellCount_ = cellCount_ + block.ncells_;
                }
            }
        }

        // Write polyMesh here using globalPoints_, globalFaces_, globalOwners_, globalNeighbours_, and boundaryFaces_
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

    void backgroundMesh::addFaces(const pointField &blockPoints, const faceList &blockFaces, const wordList &blockPatches, const wordList &stringPtrs)
    {
        label currentCell = cellCount_;
        int count = 0;

        for (const auto &face : blockFaces)
        {
            // Adjust face points according to the globalPoints_
            Foam::face globalFace;
            for (const auto &pt : face)
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

            Foam::face globalFaceCopy = Foam::face(globalFace);

            // Sort the face points to ensure unique representation (important for shared faces)
            std::sort(globalFaceCopy.begin(), globalFaceCopy.end());

            // Check if the face exists in the boundaryFaceMap_ (first encounter)
            if (boundaryFaceMap_.found(globalFaceCopy))
            {
                // The face is a boundary face but now we encounter it again
                // Pop the face from boundaryFaces_ and add it to globalFaces_ (becomes internal face)
                label faceOwnerIndex = boundaryFaceMap_[globalFaceCopy];

                label ownerFaceLocation = boundaryFaces_.find(globalFace);

                globalFaces_.append(boundaryFaces_[ownerFaceLocation]); // Add to global faces
                globalOwners_.append(faceOwnerIndex);                   // Keep the previous owner
                globalNeighbours_.append(currentCell);                  // Set neighbor for internal face

                // Remove from boundaryFaces_
                boundaryFaces_.remove(ownerFaceLocation);
                boundaryPatches_.remove(ownerFaceLocation);

                boundaryFaceMap_.erase(globalFaceCopy); // Remove from boundaryFaceMap_
            }
            else if (faceMap_.found(globalFace))
            {
                // If face is found in the global faces, it is a fully shared face (already internal)
                // label faceIndex = faceMap_[globalFace];
                // globalNeighbours_[faceIndex] = currentCell;
            }
            else
            {
                // New boundary face, add to boundaryFaces_ and boundaryFaceMap_
                boundaryFaceMap_.insert(globalFaceCopy, currentCell);
                boundaryFaces_.append(globalFace);

                boundaryPatches_.append(blockPatches[count]);
                // globalOwners_.append(currentCell); // Set owner for this boundary face

                if (!(stringPtrs[count] == word::null))
                {
                    if (stringPtrMap_.found(stringPtrs[count]))
                    {
                        stringPtrMap_[stringPtrs[count]].append(globalFace);
                    }
                    else
                    {
                        faceList newFaceList;
                        newFaceList.append(globalFace);
                        stringPtrMap_.insert(stringPtrs[count], newFaceList);
                    }
                }
            }
            count++;
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
        boundaryPatches_.clear();
        pointMap_.clear();
        faceMap_.clear();
        boundaryFaceMap_.clear();
        stringPtrMap_.clear();
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
