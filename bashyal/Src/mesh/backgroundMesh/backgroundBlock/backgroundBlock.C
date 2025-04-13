// #include "backgroundBlock.H"
#include "backgroundMesh.H"
#include "patchTypes.H"

using namespace Foam;
namespace Bashyal
{
    backgroundBlock::backgroundBlock(backgroundMesh *ref, const Vector<int> identity, const boundBox &bounds, label blockID)
        : ref_(ref), identity_(identity), bounds_(bounds), ncells_(1), blockID_(blockID)
    {
        this->generateCubeGeometry();
    }

    void backgroundBlock::generateCubeGeometry()
    {
        const point &minPt = bounds_.min(); // Minimum bound point
        const point &maxPt = bounds_.max(); // Maximum bound point

        // Define the 8 vertices of the cube (corner points)
        points_.setSize(8);
        points_[0] = minPt;                                  // (xmin, ymin, zmin)
        points_[1] = point(minPt.x(), maxPt.y(), minPt.z()); // (xmin, ymax, zmin)
        points_[2] = point(maxPt.x(), maxPt.y(), minPt.z()); // (xmax, ymax, zmin)
        points_[3] = point(maxPt.x(), minPt.y(), minPt.z()); // (xmax, ymin, zmin)
        points_[4] = point(minPt.x(), minPt.y(), maxPt.z()); // (xmin, ymin, zmax)
        points_[5] = point(minPt.x(), maxPt.y(), maxPt.z()); // (xmin, ymax, zmax)
        points_[6] = maxPt;                                  // (xmax, ymax, zmax)
        points_[7] = point(maxPt.x(), minPt.y(), maxPt.z()); // (xmax, ymin, zmax)

        // Define the 6 faces of the cube (each face with 4 vertices in counterclockwise order)
        faces_.setSize(6);
        patches_.setSize(6);

        // Bottom face (zmin): counterclockwise when viewed from below
        faces_[0] = face({0, 1, 2, 3});
        patches_[0] = patchType::XY_Zmin;

        // Top face (zmax): counterclockwise when viewed from above
        faces_[1] = face({4, 7, 6, 5});
        patches_[1] = patchType::XY_Zmax;

        // Left face (xmin): counterclockwise when viewed from the left
        faces_[2] = face({0, 4, 5, 1});
        patches_[2] = patchType::YZ_Xmin;

        // Right face (xmax): counterclockwise when viewed from the right
        faces_[3] = face({3, 2, 6, 7});
        patches_[3] = patchType::YZ_Xmax;

        // Front face (ymin): counterclockwise when viewed from the front
        faces_[4] = face({0, 3, 7, 4});
        patches_[4] = patchType::XZ_Ymin;

        // Back face (ymax): counterclockwise when viewed from the back
        faces_[5] = face({1, 5, 6, 2});
        patches_[5] = patchType::XZ_Ymax;

        owners_ = labelList{0, 0, 0, 0, 0, 0};           // All faces owned by a single cell
        neighbours_ = labelList{-1, -1, -1, -1, -1, -1}; // No neighbor cells (external boundary)
        nboundaries_ = 6;

        generateNefPolyhedron(); // Generate Nef polyhedron for the cube
    }

    void backgroundBlock::reset()
    {
        this->generateCubeGeometry();
        this->edited_ = false;
        this->dead_ = false;
    }

    bool backgroundBlock::contains(const point &pt) const
    {
        const point &minPt = bounds_.min();
        const point &maxPt = bounds_.max();

        return (pt.x() >= minPt.x() && pt.x() <= maxPt.x()) &&
               (pt.y() >= minPt.y() && pt.y() <= maxPt.y()) &&
               (pt.z() >= minPt.z() && pt.z() <= maxPt.z());
    }

    void backgroundBlock::triangulateFaces(
        const pointField &points, // Input points (not modified)
        const faceList &faces,    // Input faces to triangulate
        pointField &outPoints,    // Output points (same as input)
        faceList &outFaces        // Output triangulated faces
    )
    {
        // Assign input points to output (no modification occurs)
        outPoints = points;

        // Step 1: Calculate total number of triangles
        label totalTriangles = 0;
        forAll(faces, faceI)
        {
            const face &f = faces[faceI];
            if (f.size() >= 3) // Only process faces with 3 or more vertices
            {
                totalTriangles += f.size() - 2; // n-2 triangles per n-sided polygon
            }
        }

        // Step 2: Preallocate outFaces with the total number of triangles
        outFaces.setSize(totalTriangles);

        // Step 3: Triangulate each face and store in outFaces
        label triIndex = 0;
        forAll(faces, faceI)
        {
            const face &f = faces[faceI];
            if (f.size() < 3)
                continue; // Skip invalid faces
            // Fan triangulation: connect vertex 0 to vertices i and i+1
            for (label i = 1; i < f.size() - 1; ++i)
            {
                face tri(3);       // Create a triangle (face with 3 vertices)
                tri[0] = f[0];     // First vertex of the fan
                tri[1] = f[i];     // Current vertex
                tri[2] = f[i + 1]; // Next vertex
                outFaces[triIndex++] = tri;
            }
        }
    }

    // void backgroundBlock::write(const std::string &meshDir) const
    // {
    //     this->ref_->writePolyMeshPlain(meshDir, points_, faces_, owners_, neighbours_, );
    // }

    backgroundBlock::~backgroundBlock()
    {
    }
}
