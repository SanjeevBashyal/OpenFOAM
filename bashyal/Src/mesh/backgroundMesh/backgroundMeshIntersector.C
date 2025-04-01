#include "backgroundMesh.H"

using namespace Foam;
namespace Bashyal
{

    void backgroundMesh::intersect(aggregate &agg)
    {
        boundBox bounds = agg.getBoundBox();
        Vector<int> minIndex, maxIndex;
        this->getBlockIndexRange(bounds, minIndex, maxIndex);

        for (int i = minIndex.x(); i <= maxIndex.x(); ++i)
        {
            for (int j = minIndex.y(); j <= maxIndex.y(); ++j)
            {
                for (int k = minIndex.z(); k <= maxIndex.z(); ++k)
                {
                    backgroundBlock &block = *backgroundBlocks_[i][j][k];
                    agg.intersectWithBlock(block);
                }
            }
        }
    }

    void backgroundMesh::intersectCubes(cubeAggregates &cubeAggs)
    {
        int count = 0;
        for (cubeAggregate &cubeAgg : cubeAggs.sediments_)
        {
            this->intersect(cubeAgg);
            count++;
        }
    }

}
