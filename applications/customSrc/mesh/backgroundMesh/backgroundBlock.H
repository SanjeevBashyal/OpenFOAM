#ifndef backgroundBlock_H
#define backgroundBlock_H

#include "quickMesh.H"

using namespace Foam;
namespace Bashyal
{
    class backgroundMesh;

    class backgroundBlock
    {
    public:
        backgroundMesh *ref_;  // Reference to another backgroundBlock (if needed)
        Vector<int> identity_; // Identifier for the block
        boundBox bounds_;      // Bounding box for the block

        label ncells_;
        label blockID_;

        bool dead_ = false;
        bool edited_ = false;

        pointField points_;    // Stores the vertices of the cube
        faceList faces_;       // Stores the six faces of the cube
        labelList owners_;     // Stores owner cells
        labelList neighbours_; // Stores neighbor cells

    private:
        /* data */
    public:
        backgroundBlock(backgroundMesh *ref, const Vector<int> identity, const boundBox &bounds, label blockID);
        void generateCubeGeometry();
        bool contains(const point &pt) const;
        void write(const std::string &meshDir) const;

        label getBlockID() { return blockID_; }
        pointField getPoints() { return points_; }
        faceList getFaces() { return faces_; }
        labelList getOwners() { return owners_; }
        labelList getNeighbours() { return neighbours_; }

        ~backgroundBlock();
    };
}

#endif
