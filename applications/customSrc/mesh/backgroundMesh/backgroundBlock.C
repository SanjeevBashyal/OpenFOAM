#include "backgroundBlock.H"

namespace Bashyal
{
    backgroundBlock::backgroundBlock(backgroundMesh* ref, const label identity, const boundBox &bounds)
        : ref_(ref), identity_(identity), bounds_(bounds)
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
        points_[1] = point(maxPt.x(), minPt.y(), minPt.z()); // (xmax, ymin, zmin)
        points_[2] = point(maxPt.x(), maxPt.y(), minPt.z()); // (xmax, ymax, zmin)
        points_[3] = point(minPt.x(), maxPt.y(), minPt.z()); // (xmin, ymax, zmin)
        points_[4] = point(minPt.x(), minPt.y(), maxPt.z()); // (xmin, ymin, zmax)
        points_[5] = point(maxPt.x(), minPt.y(), maxPt.z()); // (xmax, ymin, zmax)
        points_[6] = maxPt;                                  // (xmax, ymax, zmax)
        points_[7] = point(minPt.x(), maxPt.y(), maxPt.z()); // (xmin, ymax, zmax)

        // Define the 6 faces of the cube (each face with 4 vertices in counterclockwise order)
        faces_.setSize(6);

        // Bottom face (zmin): counterclockwise when viewed from below
        faces_[0] = face({0, 1, 2, 3}); 

        // Top face (zmax): counterclockwise when viewed from above
        faces_[1] = face({4, 7, 6, 5}); 

        // Left face (xmin): counterclockwise when viewed from the left
        faces_[2] = face({0, 4, 5, 1}); 

        // Right face (xmax): counterclockwise when viewed from the right
        faces_[3] = face({2, 6, 7, 3}); 

        // Front face (ymin): counterclockwise when viewed from the front
        faces_[4] = face({0, 3, 7, 4}); 

        // Back face (ymax): counterclockwise when viewed from the back
        faces_[5] = face({1, 5, 6, 2}); 


        owners_ = labelList{0, 0, 0, 0, 0, 0};           // All faces owned by a single cell
        neighbours_ = labelList{-1, -1, -1, -1, -1, -1}; // No neighbor cells (external boundary)
    }

    bool backgroundBlock::contains(const point &pt) const
    {
        const point &minPt = bounds_.min();
        const point &maxPt = bounds_.max();

        return (pt.x() >= minPt.x() && pt.x() <= maxPt.x()) &&
               (pt.y() >= minPt.y() && pt.y() <= maxPt.y()) &&
               (pt.z() >= minPt.z() && pt.z() <= maxPt.z());
    }

    void backgroundBlock::write(const std::string &meshDir) const
    {
        this->ref_->writePolyMeshPlain(meshDir, points_, faces_, owners_, neighbours_);
    }

    backgroundBlock::~backgroundBlock()
    {
    }
}
