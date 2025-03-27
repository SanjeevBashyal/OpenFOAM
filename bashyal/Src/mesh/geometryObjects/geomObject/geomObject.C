#include "geomObject.H"

namespace Bashyal
{
    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
    // Empty constructor
    geomObject::geomObject()
        : points_(), // Default-initialized to empty pointField
          faces_()   // Default-initialized to empty faceList
    {
    }

    geomObject::geomObject(const Foam::pointField &points, const Foam::faceList &faces)
        : points_(points),
          faces_(faces)
    {
        // Basic validation
        if (points_.empty() || faces_.empty())
        {
            FatalErrorInFunction
                << "Boundary initialized with empty points or faces"
                << abort(Foam::FatalError);
        }
    }

    geomObject::geomObject(const Foam::fileName &stlFile)
    {
        // Use triSurface to read the STL file
        Foam::triSurface surface(stlFile);

        // Assign points and faces from the surface
        points_ = surface.points();
        faces_ = surface.faces();

        // Basic validation
        if (points_.empty() || faces_.empty())
        {
            FatalErrorInFunction
                << "Boundary initialized with empty points or faces from STL file: " << stlFile
                << abort(Foam::FatalError);
        }
    }

    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    geomObject::~geomObject()
    {
    }

    // * * * * * * * * * * * * * * * * Methods  * * * * * * * * * * * * * * * //

    void geomObject::writeToSTL(const Foam::fileName &stlFile) const
    {
        faceOperations::writeToSTL(points_, faces_, stlFile);
    }

}
