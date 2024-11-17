#ifndef backgroundBlock_H
#define backgroundBlock_H

#include "quickMesh.H"
#include "pointFaceHit.H"
#include "baseClass.H"

using namespace Foam;
namespace Bashyal
{
    class backgroundMesh;

    class backgroundBlock : public baseClass
    {
    public:
        backgroundMesh *ref_;  // Reference to backgroundMesh
        Vector<int> identity_; // Identifier for the block
        boundBox bounds_;      // Bounding box for the block

        label ncells_;
        label blockID_;

        bool dead_ = false;
        bool edited_ = false;

        pointField points_; // Stores the vertices of the cube
        faceList faces_;    // Stores the six faces of the cube
        wordList patches_;
        wordList stringPtrs_;
        labelList aggregateFaceMap_;

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
        wordList getPatches() { return patches_; }
        wordList getStringPtrs() { return stringPtrs_; }

        labelList getOwners() { return owners_; }
        labelList getNeighbours() { return neighbours_; }

        void intersectClosedSurface(const faceList &faces, const pointField &points, point insidePoint, word stringPtr);
        void generateBlockFace(List<point> &pts, HashTable<pointFaceHit, point> &hitMap, int faceIndex, const faceList &faces, const Foam::pointField &points, face &outputFace, pointField &outputPoints);
        void intersectSurfaceFace(const face &faceI, const int faceIndex, const pointField &points, point targetPoint, HashTable<pointFaceHit, point> &hitMap, List<List<point>> &blockFaceHitPoints, List<point> &cutFaceHitPointsI);
        void findIntersections(const face &face1, const int faceIndex, const face &face2, const UList<point> &points1, const UList<point> &points2, List<point> &pts, HashTable<pointFaceHit, point> &hitMap, List<List<point>> &blockFaceHitPoints);
        void shootSurfaceRays(const face &srcFace, const int faceIndex, const UList<point> &srcPoints, const face &tgtFace, const UList<point> &tgtPoints, List<point> &pts, HashTable<pointFaceHit, point> &hitMap, List<List<point>> &blockFaceHitPoints);
        void shootBlockRays(const face &srcFace, const UList<point> &srcPoints, const face &tgtFace, const int faceIndex, const UList<point> &tgtPoints, List<point> &pts, HashTable<pointFaceHit, point> &hitMap, List<List<point>> &blockFaceHitPoints);
        void generateCutFace(List<point> pts, HashTable<pointFaceHit, point> &hitMap, int faceIndex, face &outputFace, pointField &outputPoints, point &targetPoint);

        bool isPointInsideSurface(const point &p, const faceList &faces, const pointField &points);
        int countVerticesInsideBlockBounds(const face &face1, const pointField &points);
        int countVerticesInsideSurface(const pointField &pts, const faceList &faces, const pointField &points);

        void addPoints(const pointField &points, Foam::HashTable<label, point> &pointMap, pointField &mergedPoints);
        void addFace(const pointField &points, const face &faceI, Foam::HashTable<label, point> &pointMap, faceList &mergedFaces, pointField &mergedPoints);

        bool pointOnList(const Foam::List<Foam::point> &pointList, const Foam::point &checkPoint, Foam::point &outputPoint, const double tolerance);
        bool foundOnList(const Foam::List<Foam::point> &pointList, const Foam::point &checkPoint, const double tolerance);
        bool arePointsSame(const Foam::point &point1, const Foam::point &point2, const double tolerance);

        ~backgroundBlock();
    };
}

#endif
