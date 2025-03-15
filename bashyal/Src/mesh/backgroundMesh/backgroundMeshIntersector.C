#include "backgroundMesh.H"

using namespace Foam;
namespace Bashyal
{

    void backgroundMesh::intersectCube(cubeAggregate &cubeAgg)
    {
        boundBox bounds = cubeAgg.getBoundBox();
        Vector<int> minIndex = this->getBlockIndexContainingPoint(bounds.min());
        Vector<int> maxIndex = this->getBlockIndexContainingPoint(bounds.max());
        // const faceList &cubeFaces = cubeAgg.faces_;

        // Iterate only over the relevant blocks within minIndex and maxIndex bounds
        for (int i = minIndex.x(); i <= maxIndex.x(); ++i)
        {
            for (int j = minIndex.y(); j <= maxIndex.y(); ++j)
            {
                for (int k = minIndex.z(); k <= maxIndex.z(); ++k)
                {
                    backgroundBlock &block = *backgroundBlocks_[i][j][k];
                    block.intersectClosedSurface(cubeAgg.faces_, cubeAgg.globalPoints_, cubeAgg.getCentroid(), cubeAgg.identifier_);
                }
            }
        }
    }

    void backgroundMesh::intersectCubeCGAL(cubeAggregate &cubeAgg)
    {
        boundBox bounds = cubeAgg.getBoundBox();
        Vector<int> minIndex = this->getBlockIndexContainingPoint(bounds.min());
        Vector<int> maxIndex = this->getBlockIndexContainingPoint(bounds.max());
        // const faceList &cubeFaces = cubeAgg.faces_;

        // Iterate only over the relevant blocks within minIndex and maxIndex bounds
        for (int i = minIndex.x(); i <= maxIndex.x(); ++i)
        {
            for (int j = minIndex.y(); j <= maxIndex.y(); ++j)
            {
                for (int k = minIndex.z(); k <= maxIndex.z(); ++k)
                {
                    backgroundBlock &block = *backgroundBlocks_[i][j][k];
                    block.intersectClosedSurfaceCGAL(cubeAgg.faces_, cubeAgg.globalPoints_, cubeAgg.identifier_);
                }
            }
        }
    }

    void backgroundMesh::intersectRound(roundAggregate &roundAgg)
    {
        boundBox bounds = roundAgg.getBoundBox();
        Vector<int> minIndex = this->getBlockIndexContainingPoint(bounds.min());
        Vector<int> maxIndex = this->getBlockIndexContainingPoint(bounds.max());
        // const faceList &cubeFaces = cubeAgg.faces_;

        // Iterate only over the relevant blocks within minIndex and maxIndex bounds
        for (int i = minIndex.x(); i <= maxIndex.x(); ++i)
        {
            for (int j = minIndex.y(); j <= maxIndex.y(); ++j)
            {
                for (int k = minIndex.z(); k <= maxIndex.z(); ++k)
                {
                    backgroundBlock &block = *backgroundBlocks_[i][j][k];
                    block.intersectClosedSurface(roundAgg.faces_, roundAgg.globalPoints_, roundAgg.getCentroid(), roundAgg.identifier_);
                }
            }
        }

        // Process intersections (store, visualize, etc.)
        // ...
    }

    void backgroundMesh::intersectCubes(cubeAggregates &cubeAggs)
    {
        int count = 0;
        for (cubeAggregate &cubeAgg : cubeAggs.sediments_)
        {
            this->intersectCube(cubeAgg);
            count++;
        }
    }

}
