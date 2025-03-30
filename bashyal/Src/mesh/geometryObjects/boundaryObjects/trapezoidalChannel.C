#include "trapezoidalChannel.H"

trapezoidalChannel::trapezoidalChannel(
    Foam::scalar length,
    Foam::scalar bottomWidth,
    Foam::scalar topWidth,
    Foam::scalar height,
    Foam::scalar slope
)
    : geomObject()  // Call the empty constructor of geomObject
{
    // Define the points based on trapezoidal geometry
    points_ = Foam::pointField(8);  // Allocate space for 8 points
    points_[0] = Foam::point(0, -topWidth / 2, height);              // Top-left front
    points_[1] = Foam::point(0, -bottomWidth / 2, 0);                // Bottom-left front
    points_[2] = Foam::point(0, bottomWidth / 2, 0);                 // Bottom-right front
    points_[3] = Foam::point(0, topWidth / 2, height);               // Top-right front
    points_[4] = Foam::point(length, -topWidth / 2, height - slope * length);  // Top-left back
    points_[5] = Foam::point(length, -bottomWidth / 2, -slope * length);       // Bottom-left back
    points_[6] = Foam::point(length, bottomWidth / 2, -slope * length);        // Bottom-right back
    points_[7] = Foam::point(length, topWidth / 2, height - slope * length);   // Top-right back

    // Define the faces (6 faces for a trapezoidal prism)
    faces_ = Foam::faceList(6);  // Allocate space for 6 faces
    faces_[0] = Foam::face({0, 1, 2, 3});  // Front face
    faces_[1] = Foam::face({4, 5, 6, 7});  // Back face
    faces_[2] = Foam::face({1, 2, 6, 5});  // Bottom face
    faces_[3] = Foam::face({0, 3, 7, 4});  // Top face
    faces_[4] = Foam::face({0, 1, 5, 4});  // Left side
    faces_[5] = Foam::face({3, 2, 6, 7});  // Right side
}