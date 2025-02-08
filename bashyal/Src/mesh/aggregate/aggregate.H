#ifndef aggregate_H
#define aggregate_H

#include "quickInclude.H"
#include "baseClass.H"

#include "triSurface.H"

using namespace Foam;

namespace Bashyal
{
    class aggregate : public baseClass
    {
    public:
        float backgroundFinenessIndex_ = 1; // 1 for 0.1, 2 for 0.01

        Foam::point centroid_;
        Foam::scalar xRotation_ = 0;
        Foam::scalar yRotation_ = 0;
        Foam::scalar zRotation_ = 0;

        Foam::pointField localPoints_;  // Collection of points to form the aggregate surface
        Foam::pointField globalPoints_; // Collection of points to form the aggregate surface
        Foam::boundBox boundBox_;
        Foam::boundBox roundedBoundBox_;

        Foam::faceList faces_;
        Foam::triSurface surface_; // Triangulated surface mesh

    private:
        /* data */
    public:
        aggregate();
        ~aggregate() = default;

        void translate(Foam::vector translationVector);
        void position(Foam::vector positionThis);
        pointField translatePoints(pointField points, Foam::vector translationVector);

        void rotate(scalar alpha, scalar beta, scalar gamma);
        pointField rotatePoints(Foam::tensor rotationMatrix);
        tensor rotationMatrixFromAngles();
        void locate();

        point getCentroid() { return centroid_; }
        boundBox getBoundBox();
        point roundPoint(point value);
        point floorPoint(point value);
        point ceilPoint(point value);
        scalar roundToRequiredDecimal(scalar value);
        scalar ceilToRequiredDecimal(scalar value);
        scalar floorToRequiredDecimal(scalar value);

        Foam::point sphericalToCartesian(Foam::scalar r, Foam::scalar theta, Foam::scalar phi);
        Foam::scalar distance(const Foam::point &p1, const Foam::point &p2);
    };

    inline tensor operator*(const tensor &A, const tensor &B)
    {
        tensor C = tensor::zero; // Initialize result tensor to zero

        for (int i = 0; i < 3; ++i) // Loop over rows
        {
            for (int j = 0; j < 3; ++j) // Loop over columns
            {
                for (int k = 0; k < 3; ++k) // Summing over the inner dimension
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