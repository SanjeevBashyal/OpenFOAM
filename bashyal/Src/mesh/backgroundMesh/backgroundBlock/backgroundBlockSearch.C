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

    bool backgroundBlock::pointOnList(const Foam::List<Foam::point> &pointList, const Foam::point &checkPoint, Foam::point &outputPoint, const double tolerance)
    {
        for (const Foam::point &pt : pointList)
        {
            // Calculate the distance between the points
            double distance = Foam::mag(pt - checkPoint);

            // Check if the distance is within the specified tolerance
            if (distance <= tolerance)
            {
                outputPoint = pt;
                return true;
            }
        }
        return false;
    }

    bool backgroundBlock::arePointsSame(const Foam::point &point1, const Foam::point &point2, const double tolerance)
    {
        double distance = Foam::mag(point1 - point2);

        // Check if the distance is within the specified tolerance
        if (distance <= tolerance)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}