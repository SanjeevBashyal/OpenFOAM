#ifndef backgroundMesh_H
#define backgroundMesh_H

#include "quickInclude.H"
#include "quickMesh.H"
#include "quickWriter.H"
#include "block.H"

using namespace Foam;

namespace Bashyal
{
    class backgroundMesh
    {
    private:
        /* data */
    public:
        double s_;
        Foam::Time *runTime_;
        
        Foam::pointField vertices_;
        Foam::polyMesh *meshPtr_;
        Foam::block *blockPtr_;

        

        Foam::List<Foam::face> faceListPtr_;
        Foam::List<Foam::cell> cellListPtr_;

        Foam::faceListList boundaryFaces_;
        Foam::wordList boundaryPatchNemes_;
        Foam::wordList boundaryPatchTypes_;
        Foam::word defaultBoundaryPatchName_;
        Foam::word defaultBoundaryPatchType_;
        Foam::wordList boundaryPatchPhysicalTypes_;

    public:
        explicit backgroundMesh(Foam::Time *runTime, double size);
        Foam::pointField createVertices();
        void createFacesAndCells();
        void createBoundaryFacesAndPatches();
        ~backgroundMesh() = default;

        static void writePolyMeshPlain(const std::string &meshDir, const pointField &points, const faceList &faces, const labelList &owners, const labelList &neighbours);
        void writePolyMeshFromOwnerNeighbour(const std::string &meshDir, const pointField &points, const faceList &faces, const labelList &owners, const labelList &neighbours);
    };

}

#endif
