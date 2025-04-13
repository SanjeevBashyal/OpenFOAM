#include "backgroundMesh.H"

using namespace Foam;
namespace Bashyal
{

    void backgroundMesh::intersectDomainBoundary(const boundary &domainBoundary, bool keepInside)
    {
        word operation = keepInside ? "intersection" : "difference";
        Info << "Performing domain boundary " << operation << " using boundary '"
             << domainBoundary.name() << "' on all background blocks..." << Foam::endl;

        if (domainBoundary.nef_.is_empty())
        {
            WarningInFunction << "Domain boundary '" << domainBoundary.name() << "' Nef is empty. Skipping operation." << Foam::endl;
            return;
        }

        // Loop through all blocks
        for (label i = 0; i < dim_[0]; i++)
        {
            for (label j = 0; j < dim_[1]; j++)
            {
                for (label k = 0; k < dim_[2]; k++)
                {
                    backgroundBlock &block = *backgroundBlocks_[i][j][k];
                    if (!block.isDead()) // Only intersect non-dead blocks
                    {
                        block.intersectBoundary(domainBoundary, keepInside);
                    }
                }
            }
        }
        Info << "Domain boundary " << operation << " complete." << Foam::endl;
    }

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
