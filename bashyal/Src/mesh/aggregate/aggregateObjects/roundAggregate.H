#ifndef roundAggregate_H
#define roundAggregate_H

#include "quickInclude.H"
#include "aggregate.H"

#include "triSurface.H"
#include "triSurfaceMesh.H"

namespace Bashyal
{
    class roundAggregate : public aggregate
    {
    public:
        float radius_;
        int resolution_;    // Number of divisions in the sphere
        Foam::label identifier_;

    public:
        roundAggregate();
        explicit roundAggregate(float radius, int resolution, Foam::label identifier);
        ~roundAggregate() = default;

        void createFaces();                       // Generate quadrilateral faces
        Foam::List<Foam::face> createQuadFaces(); // Generate faces for spherical structure
        void createTriSurface();                 // Create triangulated surface for STL output
        Foam::pointField createSpherePoints(Foam::scalar radius, int resolution); // Generate points for sphere
    };
}

#endif

// ************************************************************************* //
