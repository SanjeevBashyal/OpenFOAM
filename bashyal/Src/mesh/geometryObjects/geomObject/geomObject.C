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

    Bashyal::geomObject::geomObject(const Foam::fileName &stlFile)
    {
        // Use triSurface to read the STL file
        Foam::triSurface surface(stlFile);

        // Assign points from the surface
        points_ = surface.points();

        // Get the list of triFaces using the faces() method
        const Foam::List<Foam::labelledTri> &trie = surface.surfFaces();

        // Convert triSurface faces to faceList
        faces_.setSize(trie.size());
        forAll(trie, i)
        {
            const Foam::triFace &tri = trie[i];
            faces_[i] = Foam::face(3);
            faces_[i][0] = tri[0];
            faces_[i][1] = tri[1];
            faces_[i][2] = tri[2];
        }

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

    Foam::boundBox geomObject::createBoundBox()
    {
        if (points_.size() > 0)
        {
            return Foam::boundBox(points_);
        }
        else
        {
            FatalErrorInFunction
                << "Cannot create bounding box from empty points_"
                << abort(Foam::FatalError);
        }
        return Foam::boundBox(); // Return an empty bounding box if points_ is empty
    }

}
