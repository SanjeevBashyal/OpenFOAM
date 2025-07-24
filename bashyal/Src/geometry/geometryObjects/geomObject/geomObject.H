#ifndef geomObject_H
#define geomObject_H

#include "quickInclude.H"
#include "triSurface.H"
#include "boundBox.H"
#include "debugClass.H"
#include "faceOperations.H"
#include "indexedFaceSet.H"

namespace Bashyal
{
    class geomObject : public debugClass, public Foam::particleModels::indexedFaceSet
    {
    public:
        // Constructors
        geomObject();
        geomObject(const Foam::pointField &points, const Foam::faceList &faces);
        geomObject(const Foam::fileName &stlFile); // New constructor for STL file

        // Destructor
        ~geomObject();

        // Access methods
        // points() and faces() are inherited from indexedFaceSet

        void writeToSTL(const Foam::fileName &stlFile) const;
        Foam::boundBox createBoundBox();
        void writeVtp(const Foam::fileName &filename) const;
    };
}

#endif
