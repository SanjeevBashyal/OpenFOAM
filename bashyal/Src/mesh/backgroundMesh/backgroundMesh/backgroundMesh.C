#include "backgroundMesh.H"
#include "quickMesh.H"

using namespace Foam;

namespace Bashyal
{
    backgroundMesh::backgroundMesh(Foam::Time *runTime)
        : runTime_(runTime)
    {
        // Read the backgroundMeshDict from the constant directory
        IOdictionary bgMeshDict(
            IOobject(
                "backgroundMeshDict",
                runTime_->system(),
                *runTime_, // Fixed: Dereference Time* to objectRegistry&
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE));

        // Determine the method for setting mesh bounds
        word method = bgMeshDict.getOrDefault<word>("method", "manual");

        if (method == "manual")
        {
            // Read minPoint and maxPoint directly from the dictionary
            meshMin_ = bgMeshDict.get<Foam::point>("minPoint"); // Fixed: Use get<T>
            meshMax_ = bgMeshDict.get<Foam::point>("maxPoint"); // Fixed: Use get<T>
        }
        else
        {
            FatalErrorInFunction
                << "Invalid method '" << method << "' in backgroundMeshDict. "
                << "Valid options are 'manual' or 'fromSTL'."
                << exit(FatalError);
        }

        // Read the resolution (required for both methods)
        resolution_ = bgMeshDict.get<Foam::scalar>("resolution"); // Fixed: Use get<T>

        // Initialize the remaining members
        dim_ = countBlocksPerAxis();
        backgroundBlocks_ = createListPointers();
        vertices_ = createVertices();

        // Populate the backgroundBlocks_ array
        label count = 0;
        for (label i = 0; i < dim_[0]; i++) // Fixed: Use label instead of scalar
        {
            for (label j = 0; j < dim_[1]; j++)
            {
                for (label k = 0; k < dim_[2]; k++)
                {
                    point min = point(
                        meshMin_.x() + i * resolution_,
                        meshMin_.y() + j * resolution_,
                        meshMin_.z() + k * resolution_);
                    point max = point(
                        min.x() + resolution_,
                        min.y() + resolution_,
                        min.z() + resolution_);
                    backgroundBlocks_[i][j][k] = autoPtr<backgroundBlock>(
                        new backgroundBlock(
                            this,
                            Foam::Vector<int>(i, j, k),
                            Foam::boundBox(min, max),
                            count));
                    count++;
                }
            }
        }
    }

    void backgroundMesh::resetBlocks()
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

                    if (block.edited_)
                    {
                        block.reset();
                    }
                }
            }
        }
    }

    void backgroundMesh::setBoundaryPatchType(Foam::dictionary &boundaryDict)
    {
        boundaryDict_ = boundaryDict;
    }

    inline Foam::pointField backgroundMesh::createVertices()
    {
        Foam::pointField Nodes(8);
        Nodes[0] = Foam::point(meshMin_);                                 // Vertex 0: (x0, y0, z0)
        Nodes[1] = Foam::point(meshMax_.x(), meshMin_.y(), meshMin_.z()); // Vertex 1: (x1, y0, z0)
        Nodes[2] = Foam::point(meshMax_.x(), meshMax_.y(), meshMin_.z()); // Vertex 2: (x1, y1, z0)
        Nodes[3] = Foam::point(meshMin_.x(), meshMax_.y(), meshMin_.z()); // Vertex 3: (x0, y1, z0)
        Nodes[4] = Foam::point(meshMin_.x(), meshMin_.y(), meshMax_.z()); // Vertex 4: (x0, y0, z1)
        Nodes[5] = Foam::point(meshMax_.x(), meshMin_.y(), meshMax_.z()); // Vertex 5: (x1, y0, z1)
        Nodes[6] = Foam::point(meshMax_);                                 // Vertex 7: (x1, y1, z1)
        Nodes[7] = Foam::point(meshMin_.x(), meshMax_.y(), meshMax_.z()); // Vertex 6: (x0, y1, z1)

        return Nodes;
    }

    inline Vector<int> backgroundMesh::countBlocksPerAxis() const
    {
        scalar lengthX = meshMax_.x() - meshMin_.x();
        scalar lengthY = meshMax_.y() - meshMin_.y();
        scalar lengthZ = meshMax_.z() - meshMin_.z();

        // Calculate the number of cubes along each axis
        label numCubesX = std::ceil(lengthX / resolution_);
        label numCubesY = std::ceil(lengthY / resolution_);
        label numCubesZ = std::ceil(lengthZ / resolution_);

        return Vector<int>(numCubesX, numCubesY, numCubesZ);
    }

    inline List<List<List<autoPtr<backgroundBlock>>>> backgroundMesh::createListPointers()
    {
        List<List<List<autoPtr<backgroundBlock>>>> nestedList(dim_[0]); // Outer list with size 2

        // Step 2: Initialize inner lists
        for (label i = 0; i < dim_[0]; i++)
        {
            nestedList[i].setSize(dim_[1]); // Middle list with size 2
            for (label j = 0; j < dim_[1]; j++)
            {
                nestedList[i][j].setSize(dim_[2]); // Innermost list with size 2
            }
        }

        return nestedList;
    }

    void backgroundMesh::getBlockIndexRange(const boundBox &bounds, Vector<int> &minIndex, Vector<int> &maxIndex)
    {
        // Calculate raw minimum indices
        int i_min = std::floor((bounds.min().x() - meshMin_.x()) / resolution_);
        int j_min = std::floor((bounds.min().y() - meshMin_.y()) / resolution_);
        int k_min = std::floor((bounds.min().z() - meshMin_.z()) / resolution_);

        // Calculate raw maximum indices
        int i_max = std::floor((bounds.max().x() - meshMin_.x()) / resolution_);
        int j_max = std::floor((bounds.max().y() - meshMin_.y()) / resolution_);
        int k_max = std::floor((bounds.max().z() - meshMin_.z()) / resolution_);

        // Clamp indices to valid range
        minIndex.x() = std::max(0, i_min);
        minIndex.y() = std::max(0, j_min);
        minIndex.z() = std::max(0, k_min);

        maxIndex.x() = std::min(dim_[0] - 1, i_max);
        maxIndex.y() = std::min(dim_[1] - 1, j_max);
        maxIndex.z() = std::min(dim_[2] - 1, k_max);
    }

    bool backgroundMesh::contains(const point &pt)
    {
        return (pt.x() >= meshMin_.x() && pt.x() <= meshMax_.x()) &&
               (pt.y() >= meshMin_.y() && pt.y() <= meshMax_.y()) &&
               (pt.z() >= meshMin_.z() && pt.z() <= meshMax_.z());
    }

}
