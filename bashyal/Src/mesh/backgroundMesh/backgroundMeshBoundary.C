#include "backgroundMesh.H"

namespace Bashyal
{
    bool backgroundMesh::isFaceMeshBoundary(const Foam::Vector<int> &identity, int patchType)
    {
        // Check minimum boundaries
        if (patchType == -1)
        {
            if (identity.x() == 0)
                return true;

            return false;
        }
        else if (patchType == -2)
        {
            if (identity.y() == 0)
                return true;

            return false;
        }
        else if (patchType == -3)
        {
            if (identity.z() == 0)
                return true;

            return false;
        }
        // Check maximum boundaries
        else if (patchType == 1)
        {
            if (identity.x() == dim_[0]-1)
                return true;

            return false;
        }
        else if (patchType == 2)
        {
            if (identity.y() == dim_[1]-1)
                return true;

            return false;
        }
        else if (patchType == 3)
        {
            if (identity.z() == dim_[2]-1)
                return true;

            return false;
        }
        // Check if patchType is not in {-1, -2, -3, 0, 1, 2, 3}
        else if (patchType != 0)
        {
            // Since -1, -2, -3, 1, 2, 3 are already handled, this catches all other non-zero values
            return true;
        }

        return false;
    }
}
