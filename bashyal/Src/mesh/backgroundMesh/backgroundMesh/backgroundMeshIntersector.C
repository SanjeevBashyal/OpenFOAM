#include "backgroundMesh.H"

namespace Bashyal
{

    void backgroundMesh::intersectDomainBoundary(const boundary &domainBoundary, bool keepInside)
    {
        Foam::Info << "Performing domain boundary intersection/difference using boundary '"
             << domainBoundary.name() << "' on all background blocks..." << Foam::endl;

        if (domainBoundary.nef_.is_empty())
        {
            WarningInFunction << "Domain boundary '" << domainBoundary.name() << "' Nef is empty. Skipping operation." << Foam::endl;
            return;
        }

        // Loop through all blocks
        for (Foam::label i = 0; i < dim_[0]; i++)
        {
            for (Foam::label j = 0; j < dim_[1]; j++)
            {
                for (Foam::label k = 0; k < dim_[2]; k++)
                {
                    backgroundBlock &block = *backgroundBlocks_[i][j][k];
                    if (!block.dead_) // Only intersect non-dead blocks
                    {
                        block.intersectBoundary(domainBoundary, keepInside);
                    }
                }
            }
        }
        Foam::Info << "Domain boundary intersection complete." << Foam::endl;
    }

    void backgroundMesh::intersect(aggregate &agg)
    {
        Foam::boundBox bounds = agg.getBoundBox();
        Foam::Vector<int> minIndex, maxIndex;
        this->getBlockIndexRange(bounds, minIndex, maxIndex);

        for (int i = minIndex.x(); i <= maxIndex.x(); ++i)
        {
            for (int j = minIndex.y(); j <= maxIndex.y(); ++j)
            {
                for (int k = minIndex.z(); k <= maxIndex.z(); ++k)
                {
                    backgroundBlock &block = *backgroundBlocks_[i][j][k];
                    block.subtractAggregate(agg);
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
