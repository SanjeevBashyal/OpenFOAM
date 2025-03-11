#include backgroundMesh.H

namespace Bashyal
{
    bool isFacePatchBoundary(const Foam::Vector<int> &identity, int patchType, int maxX, int maxY, int maxZ)
    {
        bool isBoundary = false;

        // Check minimum boundaries
        if (patchType == -1)
        {
            if (identity.x() == 0)
                return true;
        }
        else if (patchType == -2)
        {
            if (identity.y() == 0)
                return true;
        }
        else if (patchType == -3)
        {
            if (identity.z() == 0)
                return true;
        }
        // Check maximum boundaries
        else if (patchType == 1)
        {
            if (identity.x() == maxX)
                isBoundary = true;
        }
        else if (patchType == 2)
        {
            if (identity.y() == maxY)
                isBoundary = true;
        }
        else if (patchType == 3)
        {
            if (identity.z() == maxZ)
                isBoundary = true;
        }
        // Check if patchType is not in {-1, -2, -3, 0, 1, 2, 3}
        else if (patchType != 0)
        {
            // Since -1, -2, -3, 1, 2, 3 are already handled, this catches all other non-zero values
            isBoundary = true;
        }

        return isBoundary;
    }
}
