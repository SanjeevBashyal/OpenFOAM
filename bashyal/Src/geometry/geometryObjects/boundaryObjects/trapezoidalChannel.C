#include "trapezoidalChannel.H"
#include "pointField.H"
#include "faceList.H"
#include "point.H"
#include "boundBox.H"
#include "List.H"
#include "error.H"

namespace Bashyal
{

// --- Constructor ---

trapezoidalChannel::trapezoidalChannel()
    : boundary(),
      length_(1.0),
      bottomWidth_(0.2),
      topWidth_(0.4),
      height_(0.3),
      slope_(0.01)
{
    initializeDefaultPatchTypes();
    develop();
    // generateNefPolyhedron(); // Uncomment if required by base class
}

// --- Private Methods ---

void trapezoidalChannel::initializeDefaultPatchTypes()
{
    patchTypes_.setSize(6);
    patchTypes_[0] = 10; // Face 0: Inlet (Front face)
    patchTypes_[1] = 11; // Face 1: Outlet (Back face)
    patchTypes_[2] = 13; // Face 2: Bottom face (Wall)
    patchTypes_[3] = 15; // Face 3: Top face (SymmetryPlane)
    patchTypes_[4] = 13; // Face 4: Left side face (Wall)
    patchTypes_[5] = 13; // Face 5: Right side face (Wall)
}

void trapezoidalChannel::develop()
{
    // Basic validation
    if (length_ <= 0 || bottomWidth_ <= 0 || topWidth_ <= 0 || height_ <= 0)
    {
        WarningIn("trapezoidalChannel::develop()")
            << "Developing trapezoidal channel geometry with non-positive dimensions." << Foam::endl
            << "  Length: " << length_
            << ", Bottom Width: " << bottomWidth_
            << ", Top Width: " << topWidth_
            << ", Height: " << height_ << Foam::endl;
    }

    // Define the 8 vertices of the trapezoidal prism
    Foam::pointField points(8);
    Foam::scalar zBackBottom = -slope_ * length_;
    Foam::scalar zBackTop = height_ - slope_ * length_;

    // Front face (x=0)
    points[0] = Foam::point(0, -topWidth_ / 2, height_);        // Top-left front
    points[1] = Foam::point(0, -bottomWidth_ / 2, 0);           // Bottom-left front
    points[2] = Foam::point(0, bottomWidth_ / 2, 0);            // Bottom-right front
    points[3] = Foam::point(0, topWidth_ / 2, height_);         // Top-right front

    // Back face (x=length)
    points[4] = Foam::point(length_, -topWidth_ / 2, zBackTop); // Top-left back
    points[5] = Foam::point(length_, -bottomWidth_ / 2, zBackBottom); // Bottom-left back
    points[6] = Foam::point(length_, bottomWidth_ / 2, zBackBottom);  // Bottom-right back
    points[7] = Foam::point(length_, topWidth_ / 2, zBackTop);  // Top-right back

    // Define the 6 faces using vertex indices
    Foam::faceList faces(6);
    faces[0] = Foam::face({0, 1, 2, 3}); // Front face (inlet)
    faces[1] = Foam::face({4, 5, 6, 7}); // Back face (outlet)
    faces[2] = Foam::face({1, 2, 6, 5}); // Bottom face
    faces[3] = Foam::face({0, 3, 7, 4}); // Top face
    faces[4] = Foam::face({0, 1, 5, 4}); // Left side face
    faces[5] = Foam::face({3, 2, 6, 7}); // Right side face

    // Assign to base class
    *static_cast<Foam::particleModels::indexedFaceSet*>(this) = Foam::particleModels::indexedFaceSet(points, faces);
}

// --- Public Member Functions ---

void trapezoidalChannel::initialize(
    Foam::scalar length,
    Foam::scalar bottomWidth,
    Foam::scalar topWidth,
    Foam::scalar height,
    Foam::scalar slope
)
{
    length_ = length;
    bottomWidth_ = bottomWidth;
    topWidth_ = topWidth;
    height_ = height;
    slope_ = slope;
    develop();
}

const Foam::List<int>& trapezoidalChannel::getPatchTypes() const
{
    return patchTypes_;
}

} // End namespace Bashyal