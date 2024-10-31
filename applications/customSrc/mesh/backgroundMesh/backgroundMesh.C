#include "backgroundMesh.H"
#include "quickMesh.H"

using namespace Foam;

namespace Bashyal
{
    backgroundMesh::backgroundMesh(Foam::Time *runTime, point meshMin, point meshMax, float resolution)
        : meshMin_(meshMin),
          meshMax_(meshMax),
          resolution_(resolution),
          dim_(countBlocksPerAxis()),
          backgroundBlocks_(createListPointers()),
          runTime_(runTime),
          vertices_(createVertices())
    {
        label count = 0;
        for (scalar i = 0; i < dim_[0]; i++)
        {
            for (scalar j = 0; j < dim_[1]; j++)
            {
                for (scalar k = 0; k < dim_[2]; k++)
                {
                    point min = point(meshMin_.x() + i * resolution_, meshMin_.y() + j * resolution_, meshMin_.z() + k * resolution_);
                    point max = point(min.x() + resolution_, min.y() + resolution_, min.z() + resolution_);
                    backgroundBlocks_[i][j][k] = autoPtr<backgroundBlock>(new backgroundBlock(this, Foam::Vector<int>(i, j, k), Foam::boundBox(min, max), count));
                    count = count + 1;
                }
            }
        }

        Foam::Info << "Here" << Foam::endl;
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

    Vector<int> backgroundMesh::getBlockIndexContainingPoint(point &pt)
    {
        // Ensure the point is within the bounds of the mesh
        if (!this->contains(pt))
        {
            return Vector<int>(dim_[0] - 1, dim_[1] - 1, dim_[2] - 1);
        }

        // Determine the index in each dimension by scaling the point's position
        int i = floor((pt.x() - meshMin_.x()) / resolution_);
        int j = floor((pt.y() - meshMin_.y()) / resolution_);
        int k = floor((pt.z() - meshMin_.z()) / resolution_);

        // Return the index vector (i, j, k)
        return Vector<int>(i, j, k);
    }

    bool backgroundMesh::contains(const point &pt)
    {
        return (pt.x() >= meshMin_.x() && pt.x() <= meshMax_.x()) &&
               (pt.y() >= meshMin_.y() && pt.y() <= meshMax_.y()) &&
               (pt.z() >= meshMin_.z() && pt.z() <= meshMax_.z());
    }

}
