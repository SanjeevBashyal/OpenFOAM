#include "boundary.H"

namespace Bashyal
{
    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    boundary::boundary(const Foam::pointField &points, const Foam::faceList &faces)
        : geomObject(points, faces)
    {
        // Additional boundary-specific initialization can be added here if needed
    }

    boundary::boundary(const Foam::fileName &stlFile)
        : geomObject(stlFile)
    {
        // Additional boundary-specific initialization can be added here if needed
    }

    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    boundary::~boundary()
    {
    }

}
