#ifndef cubeAggregate_H
#define cubeAggregate_H

#include "quickInclude.H"
#include "aggregate.H"

#include "triSurface.H"
#include "triSurfaceMesh.H"

using namespace Foam;

namespace Bashyal
{
    class cubeAggregate : public aggregate
    {
    public:
        float s1_;
        float s2_;
        float s_;

    private:
        /* data */
    public:
        cubeAggregate();
        explicit cubeAggregate(float s1, float s2);
        // cubeAggregate *getPtr() { return this; }

        void createFaces();
        Foam::List<Foam::face> createQuadFaces();
        void createTriSurface();
        Foam::List<Foam::labelledTri> createTriangularFacesFromPoints(const Foam::pointField &points);

        void hit(scalar index);
        List<pointField> facePoints_;

        ~cubeAggregate() = default;
    };

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //