#ifndef cubeAggregate_H
#define cubeAggregate_H

#include "quickInclude.H"
#include "aggregate.H"

#include "triSurface.H"
#include "triSurfaceMesh.H"

namespace Bashyal
{
    class cubeAggregate : public aggregate
    {
    public:
        float s1_;
        float s2_;
        float s_;
        Foam::label identifier_;

    private:
        /* data */
    public:
        cubeAggregate();
        explicit cubeAggregate(float s1, float s2, Foam::label identifier);
        // cubeAggregate *getPtr() { return this; }

        void createFaces();
        Foam::List<Foam::face> createQuadFaces();
        void createTriSurface();
        Foam::List<Foam::labelledTri> createTriangularFacesFromPoints(const Foam::pointField &points);

        void hit(Foam::scalar index);
        Foam::List<Foam::pointField> facePoints_;

        ~cubeAggregate() = default;

    };

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //