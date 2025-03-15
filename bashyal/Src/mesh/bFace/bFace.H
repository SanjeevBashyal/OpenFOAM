#ifndef bFace_H
#define bFace_H

#include "face.H"
#include "pointField.H"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Define CGAL types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;

class bFace : public face
{
private:
    Plane_3 plane;

public:
    // Constructor
    bFace(const face& f, const pointField& meshPoints);

    // Method to check if a point lies on the plane
    bool liesOnPlane(const Foam::vector& point) const;
};

#endif
