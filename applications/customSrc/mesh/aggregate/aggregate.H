#ifndef aggregate_H
#define aggregate_H

#include "aggregate.H"
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
    class aggregate
    {
        public:
            float r1;
            float r2;
            int nPoints;

            // Foam::scalar size_;    // Size of the bounding box for the aggregate
            // Foam::pointField points_;  // Collection of points to form the aggregate surface
            // Foam::triSurface surface_; // Triangulated surface mesh


        private:
            /* data */
        public:
            explicit aggregate(float r1, float r2, int nPoints);
            Foam::point sphericalToCartesian(Foam::scalar r, Foam::scalar theta, Foam::scalar phi);
            Foam::scalar distance(const Foam::point &p1, const Foam::point &p2);
            Foam::faceList createFacesFromPoints(const Foam::pointField &points);
            ~aggregate()=default;
    };

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //