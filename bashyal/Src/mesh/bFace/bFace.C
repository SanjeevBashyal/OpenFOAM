#include "bFace.H"
#include "FatalError.H"

bFace::bFace(const face& f, const pointField& meshPoints)
    : face(f)
{
    // Compute the normal vector of the face
    const Foam::vector normal = f.normal(meshPoints);
    
    // Get the points of the face
    const Foam::pointField facePoints = f.points(meshPoints);
    
    // Check if the face has at least 3 points
    if (facePoints.size() < 3)
    {
        FatalErrorIn("bFace::bFace(const face&, const pointField&)")
            << "Face has less than 3 points"
            << exit(FatalError);
    }
    
    // Convert the first point to CGAL Point_3
    Point_3 p(facePoints[0].x(), facePoints[0].y(), facePoints[0].z());
    
    // Convert the normal vector to CGAL Vector_3
    Vector_3 n(normal.x(), normal.y(), normal.z());
    
    // Initialize the plane with a point and the normal direction
    plane = Plane_3(p, CGAL::Direction_3(n));
}

bool bFace::liesOnPlane(const Foam::vector& point) const
{
    // Convert the input point to CGAL Point_3
    Point_3 p(point.x(), point.y(), point.z());
    
    // Check if the point lies on the plane
    return plane.has_on(p);
}