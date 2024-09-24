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
    public:
        explicit backgroundMesh(Foam::Time* runTime, double size);
        Foam::pointField createVertices();
        void createFacesAndCells();
        ~backgroundMesh() = default;
    };

}

#endif
