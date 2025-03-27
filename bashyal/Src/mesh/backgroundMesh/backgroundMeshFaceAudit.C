#include "backgroundMesh.H"
#include <algorithm> // For std::sort

namespace Bashyal
{
    void backgroundMesh::auditFaces()
    {
        // Loop through each block in the 3D backgroundBlocks_ list
        for (label i = 0; i < dim_[0]; i++)
        {
            for (label j = 0; j < dim_[1]; j++)
            {
                for (label k = 0; k < dim_[2]; k++)
                {
                    // Dereference the autoPtr to access the backgroundBlock
                    backgroundBlock& block = *(backgroundBlocks_[i][j][k]);
                    const Vector<int>& identity = block.identity_;

                    // Loop through each face in the block
                    forAll(block.faces_, faceI)
                    {
                        int patchType = block.patches_[faceI];

                        // Skip boundary faces (if applicable)
                        if (isFaceMeshBoundary(identity, patchType))
                        {
                            continue;
                        }

                        // Check for internal faces with patchType 1, 2, or 3
                        if (patchType == 1 || patchType == 2 || patchType == 3)
                        {
                            Vector<int> neighborIdentity = identity;
                            int expectedPatchType = -patchType; // 1 -> -1, 2 -> -2, 3 -> -3
                            backgroundBlock* neighborBlock = nullptr;

                            // Identify the neighboring block based on patchType
                            if (patchType == 1 && i < dim_[0] - 1) // Positive X direction
                            {
                                neighborIdentity[0] += 1; // i+1
                                neighborBlock = &*(backgroundBlocks_[i+1][j][k]);
                            }
                            else if (patchType == 2 && j < dim_[1] - 1) // Positive Y direction
                            {
                                neighborIdentity[1] += 1; // j+1
                                neighborBlock = &*(backgroundBlocks_[i][j+1][k]);
                            }
                            else if (patchType == 3 && k < dim_[2] - 1) // Positive Z direction
                            {
                                neighborIdentity[2] += 1; // k+1
                                neighborBlock = &*(backgroundBlocks_[i][j][k+1]);
                            }

                            // If a neighboring block exists, audit the face
                            if (neighborBlock)
                            {
                                // Find the face in the neighbor with the expected patch type
                                label neighborFaceIndex = -1;
                                forAll(neighborBlock->patches_, nFaceI)
                                {
                                    if (neighborBlock->patches_[nFaceI] == expectedPatchType)
                                    {
                                        neighborFaceIndex = nFaceI;
                                        break;
                                    }
                                }

                                if (neighborFaceIndex != -1)
                                {
                                    // Extract and sort points of the current face
                                    const face& currentFace = block.faces_[faceI];
                                    List<point> currentPoints;
                                    forAll(currentFace, ptI)
                                    {
                                        currentPoints.append(roundPoint(block.points_[currentFace[ptI]]));
                                    }
                                    std::sort(currentPoints.begin(), currentPoints.end());

                                    // Extract and sort points of the neighboring face
                                    const face& neighborFace = neighborBlock->faces_[neighborFaceIndex];
                                    List<point> neighborPoints;
                                    forAll(neighborFace, ptI)
                                    {
                                        neighborPoints.append(roundPoint(neighborBlock->points_[neighborFace[ptI]]));
                                    }
                                    std::sort(neighborPoints.begin(), neighborPoints.end());

                                    // Compare the faces
                                    if (currentPoints != neighborPoints)
                                    {
                                        Info << "Mismatch detected between block (" << i << "," << j << "," << k 
                                             << ") face " << faceI << " (patch " << patchType 
                                             << ") and neighbor block (" << neighborIdentity[0] << "," 
                                             << neighborIdentity[1] << "," << neighborIdentity[2] 
                                             << ") face " << neighborFaceIndex << " (patch " << expectedPatchType << ")" << endl;
                                    }
                                }
                                else
                                {
                                    FatalErrorInFunction 
                                        << "No face with patch type " << expectedPatchType 
                                        << " found in neighboring block at " << neighborIdentity 
                                        << exit(FatalError);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}