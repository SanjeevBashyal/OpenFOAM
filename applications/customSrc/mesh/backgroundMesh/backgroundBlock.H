#ifndef backgroundBlock_H
#define backgroundBlock_H

#include "quickMesh.H"
#include "pointFaceHit.H"
#include "cubeAggregate.H"

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

        void intersectCube(cubeAggregate &cubeAgg);
        void intersectBlockFace(const face &cubeFace, const pointField &cubePoints);
        void intersectOthersFace(const face &cubeFace, const pointField &cubePoints, face &outputFace, pointField &outputPoints, point targetPoint);
        void findIntersections(const face &face1, const face &face2, const UList<point> &points1, const UList<point> &points2, List<point> &pts, HashTable<pointFaceHit, point> &hitMap);
        void shootRays(const face &srcFace, const UList<point> &srcPoints, const face &tgtFace, const UList<point> &tgtPoints, List<point> &pts, HashTable<pointFaceHit, point> &hitMap);
        void shootEdges(const face &srcFace, const UList<point> &srcPoints, const face &tgtFace, const UList<point> &tgtPoints, List<point> &pts, HashTable<pointFaceHit, point> &hitMap);
        void generateFace(List<point> pts, HashTable<pointFaceHit, point> &hitMap, face &outputFace, pointField &outputPoints, point &targetPoint);

        bool isPointInsideSurface(const point &p);
        int countVerticesInsideBlock(const face &face, const pointField &points);
        int countVerticesInsideSurface(const faceList &faces, const pointField &points);

        ~backgroundBlock();
    };
}

#endif
