#ifndef aggregate_H
#define aggregate_H

#include "quickInclude.H"
#include "debugClass.H"
#include "backgroundBlock.H"
#include "triSurface.H"

#include <CGAL/Nef_polyhedron_3.h>
// Note: <CGAL/Exact_predicates_exact_constructions_kernel.h> is included via "foamCGALConverter.H"

namespace Bashyal
{
    class backgroundBlock;
    class aggregate : public debugClass
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

        Foam::faceList faces_;
        Foam::triSurface surface_; // Triangulated surface mesh

        // Add Nef_polyhedron_3 member with Epeck kernel
        CGAL::Nef_polyhedron_3<CGAL::Epeck> nef_;

    private:
        /* data */
    public:
        aggregate();
        ~aggregate() = default;

        void translate(Foam::vector translationVector);
        void position(Foam::vector positionThis);
        Foam::pointField translatePoints(Foam::pointField points, Foam::vector translationVector);

        void rotate(Foam::scalar alpha, Foam::scalar beta, Foam::scalar gamma);
        Foam::pointField rotatePoints(Foam::tensor rotationMatrix);
        Foam::tensor rotationMatrixFromAngles();
        void locate();

        Foam::point getCentroid() { return centroid_; }
        Foam::boundBox getBoundBox();
        Foam::point roundPoint(Foam::point value);
        Foam::point floorPoint(Foam::point value);
        Foam::point ceilPoint(Foam::point value);
        Foam::scalar roundToRequiredDecimal(Foam::scalar value);
        Foam::scalar ceilToRequiredDecimal(Foam::scalar value);
        Foam::scalar floorToRequiredDecimal(Foam::scalar value);

        Foam::point sphericalToCartesian(Foam::scalar r, Foam::scalar theta, Foam::scalar phi);
        Foam::scalar distance(const Foam::point &p1, const Foam::point &p2);

        void createTriSurface();

        // Add method declaration
        void generateNefPolyhedron();

        virtual void intersectWithBlock(backgroundBlock& block) const = 0;
    };

    inline Foam::tensor operator*(const Foam::tensor &A, const Foam::tensor &B)
    {
        Foam::tensor C = Foam::tensor::zero; // Initialize result tensor to zero

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
