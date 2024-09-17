#ifndef cubeAggregate_H
#define cubeAggregate_H

#include "point.H"
#include "scalar.H"
#include "pointField.H"
#include "List.H"

#include "Time.H"
#include "Random.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"

namespace Bashyal
{
    class cubeAggregate
    {
        public:
            float s1;
            float s2;
            float s;

            // Foam::scalar size_;    // Size of the bounding box for the aggregate
            // Foam::pointField points_;  // Collection of points to form the aggregate surface
            // Foam::triSurface surface_; // Triangulated surface mesh


        private:
            /* data */
        public:
            explicit cubeAggregate(float s1, float s2);
            Foam::List<Foam::labelledTri> createFacesFromPoints(const Foam::pointField &points);
            ~cubeAggregate()=default;
    };

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //