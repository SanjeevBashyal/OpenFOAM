#ifndef cubeAggregate_H
#define cubeAggregate_H

#include "quickInclude.H"

#include "triSurface.H"
#include "triSurfaceMesh.H"

using namespace Foam;

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
            Foam::scalar xRotation_;
            Foam::scalar yRotation_;
            Foam::scalar zRotation_;
            Foam::pointField localPoints_;  // Collection of points to form the aggregate surface
            Foam::pointField globalPoints_;  // Collection of points to form the aggregate surface
            Foam::triSurface surface_; // Triangulated surface mesh


        private:
            /* data */
        public:
            explicit cubeAggregate(float s1, float s2);

            void createSurface();
            Foam::List<Foam::labelledTri> createFacesFromPoints(const Foam::pointField &points);
            
            void translate(Foam::vector translationVector);
            pointField translatePoints(pointField points, Foam::vector translationVector);

            void rotate(scalar alpha, scalar beta, scalar gamma);
            pointField rotatePoints(Foam::tensor rotationMatrix);
            tensor rotationMatrixFromAngles();
            void locate();

            ~cubeAggregate()=default;
    };

    tensor operator*(const tensor& A, const tensor& B)
    {
        tensor C = tensor::zero;  // Initialize result tensor to zero

        for (int i = 0; i < 3; ++i)  // Loop over rows
        {
            for (int j = 0; j < 3; ++j)  // Loop over columns
            {
                for (int k = 0; k < 3; ++k)  // Summing over the inner dimension
                {
                    C(i, j) += A(i, k) * B(k, j);
                }
            }
        }
        return C;
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //