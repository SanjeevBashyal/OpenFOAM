#ifndef cubeAggregate_H
#define cubeAggregate_H

#include "quickInclude.H"

#include "triSurface.H"
#include "triSurfaceMesh.H"

namespace Bashyal
{
    class cubeAggregate
    {
        public:
            float s1_;
            float s2_;
            float s_;

            // Foam::scalar size_;    // Size of the bounding box for the aggregate
            Foam::point centroid_;
            Foam::pointField points_;  // Collection of points to form the aggregate surface
            Foam::triSurface surface_; // Triangulated surface mesh


        private:
            /* data */
        public:
            explicit cubeAggregate(float s1, float s2);
            Foam::List<Foam::labelledTri> createFacesFromPoints(const Foam::pointField &points);
            void translate(Foam::vector translationVector);
            void rotate(Foam::tensor rotationMatrix);
            ~cubeAggregate()=default;
    };

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //