#ifndef geomObject_H
#define geomObject_H

#include "quickInclude.H"
#include "faceOperations.H"

namespace Bashyal
{
    class geomObject
    {
    public:
        // Private data members
        Foam::pointField points_; // List of 3D points defining the boundary vertices
        Foam::faceList faces_;    // List of faces, each referencing indices into points_

    public:
        // Constructors
        geomObject();
        geomObject(const Foam::pointField &points, const Foam::faceList &faces);
        geomObject(const Foam::fileName &stlFile); // New constructor for STL file

        // Destructor
        ~geomObject();

        // Access methods
        const Foam::pointField &points() const { return points_; }
        const Foam::faceList &faces() const { return faces_; }

        void writeToSTL(const Foam::fileName &stlFile) const;
    };
}

#endif
