#include "backgroundBlock.H"
using namespace Foam;
namespace Bashyal
{
    void backgroundBlock::isInsideTriSurface(triSurface &surf)
    {
        triSurfaceSearch querySurf(surf);
        boolList pointInside(querySurf.calcInside(this->points_));
        for (label i = 0; i < bools.size(); ++i)
        {
            if (!bools[i])
            {
                return false; // Return false if any element is false
            }
        }
        return true; // All elements are true
    }

    void backgroundBlock::makeDead()
    {
        this->dead_ = true;
    }
}