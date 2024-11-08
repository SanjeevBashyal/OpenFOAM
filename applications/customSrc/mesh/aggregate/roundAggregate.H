#ifndef roundAggregate_H
#define roundAggregate_H

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

        float backgroundFinenessIndex_ = 1;     //1 for 0.1, 2 for 0.01

        // Foam::scalar size_;    // Size of the bounding box for the aggregate
        Foam::point centroid_;
        Foam::scalar xRotation_;
        Foam::scalar yRotation_;
        Foam::scalar zRotation_;

        Foam::pointField localPoints_;  // Collection of points to form the aggregate surface
        Foam::pointField globalPoints_; // Collection of points to form the aggregate surface
        Foam::boundBox boundBox_;
        Foam::boundBox roundedBoundBox_;

        Foam::faceList faces_;
        Foam::triSurface surface_; // Triangulated surface mesh

    private:
        /* data */
    public:
        explicit cubeAggregate(float s1, float s2);

        void createSurface();
        Foam::List<Foam::face> createQuadFaces();
        Foam::List<Foam::labelledTri> createTriangularFacesFromPoints(const Foam::pointField &points);

        void translate(Foam::vector translationVector);
        pointField translatePoints(pointField points, Foam::vector translationVector);

        void rotate(scalar alpha, scalar beta, scalar gamma);
        pointField rotatePoints(Foam::tensor rotationMatrix);
        tensor rotationMatrixFromAngles();
        void locate();

        point getCentroid() {return centroid_;}
        boundBox getBoundBox();
        point roundPoint(point value); 
        point floorPoint(point value); 
        point ceilPoint(point value); 
        scalar roundToRequiredDecimal(scalar value);
        scalar ceilToRequiredDecimal(scalar value);
        scalar floorToRequiredDecimal(scalar value);
        
        void hit(scalar index);
        List<pointField> facePoints_;

        ~cubeAggregate() = default;
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