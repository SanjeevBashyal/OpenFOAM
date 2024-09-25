#ifndef backgroundMesh_H
#define backgroundMesh_H

#include "quickInclude.H"
#include "block.H"

namespace Bashyal
{
    class backgroundMesh
    {
    private:
        /* data */
    public:
        double s_;
        Foam::pointField vertices_;
        Foam::polyMesh* meshPtr_;
        Foam::block* blockPtr_;

        Foam::List<Foam::face> faceListPtr_;
        Foam::List<Foam::cell> cellListPtr_;

        Foam::faceListList boundaryFaces_;
        Foam::wordList boundaryPatchNemes_;
        Foam::wordList boundaryPatchTypes_;
        Foam::word defaultBoundaryPatchName_;
        Foam::word defaultBoundaryPatchType_;
        Foam::wordList boundaryPatchPhysicalTypes_;
        
    public:
        explicit backgroundMesh(Foam::Time* runTime, double size);
        Foam::pointField createVertices();
        void createFacesAndCells();
        void createBoundaryFacesAndPatches();
        ~backgroundMesh() = default;
    };

}

#endif
