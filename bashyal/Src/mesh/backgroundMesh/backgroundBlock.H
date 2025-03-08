#ifndef backgroundBlock_H
#define backgroundBlock_H

#include "debugClass.H"
#include "quickMesh.H"
#include "pointFaceHit.H"
#include "baseClass.H"

namespace Bashyal
{
    class backgroundMesh;
    class backgroundBlock : public debugClass
    {
    public:
        backgroundMesh *ref_;        // Reference to backgroundMesh
        Foam::Vector<int> identity_; // Identifier for the block
        Foam::boundBox bounds_;      // Bounding box for the block

        Foam::label ncells_;
        Foam::label blockID_;

        bool dead_ = false;
        bool edited_ = false;

        Foam::pointField points_; // Stores the vertices of the cube
        Foam::faceList faces_;    // Stores the six faces of the cube
        Foam::wordList patches_;
        Foam::wordList stringPtrs_;
        Foam::labelList aggregateFaceMap_;

        Foam::labelList owners_;     // Stores owner cells
        Foam::labelList neighbours_; // Stores neighbor cells

        Foam::faceList triFaces_; // New attribute to store triangulated faces

    private:
        /* data */
    public:
        backgroundBlock(backgroundMesh *ref, const Foam::Vector<int> identity, const Foam::boundBox &bounds, Foam::label blockID);
        void generateCubeGeometry();
        void reset();
        bool contains(const Foam::point &pt) const;
        void write(const std::string &meshDir) const;

        Foam::label getBlockID() { return blockID_; }
        Foam::pointField getPoints() { return points_; }
        Foam::faceList getFaces() { return faces_; }
        Foam::wordList getPatches() { return patches_; }
        Foam::wordList getStringPtrs() { return stringPtrs_; }

        Foam::labelList getOwners() { return owners_; }
        Foam::labelList getNeighbours() { return neighbours_; }

        void intersectClosedSurface(const Foam::faceList &faces, const Foam::pointField &points, Foam::point insidePoint, Foam::word stringPtr);
        void generateBlockFace(Foam::List<Foam::point> &pts, Foam::HashTable<pointFaceHit, Foam::point> &hitMap, int faceIndex, const Foam::faceList &faces, const Foam::pointField &points, Foam::face &outputFace, Foam::pointField &outputPoints);
        void intersectSurfaceFace(const Foam::face &faceI, const int faceIndex, const Foam::pointField &points, Foam::point targetPoint, Foam::HashTable<pointFaceHit, Foam::point> &hitMap, Foam::List<Foam::List<Foam::point>> &blockFaceHitPoints, Foam::List<Foam::point> &cutFaceHitPointsI);
        void findIntersections(const Foam::face &face1, const int faceIndex, const Foam::face &face2, const Foam::UList<Foam::point> &points1, const Foam::UList<Foam::point> &points2, Foam::List<Foam::point> &pts, Foam::HashTable<pointFaceHit, Foam::point> &hitMap, Foam::List<Foam::List<Foam::point>> &blockFaceHitPoints);
        void shootSurfaceRays(const Foam::face &srcFace, const int faceIndex, const Foam::UList<Foam::point> &srcPoints, const Foam::face &tgtFace, const Foam::UList<Foam::point> &tgtPoints, Foam::List<Foam::point> &pts, Foam::HashTable<pointFaceHit, Foam::point> &hitMap, Foam::List<Foam::List<Foam::point>> &blockFaceHitPoints);
        void shootBlockRays(const Foam::face &srcFace, const Foam::UList<Foam::point> &srcPoints, const Foam::face &tgtFace, const int faceIndex, const Foam::UList<Foam::point> &tgtPoints, Foam::List<Foam::point> &pts, Foam::HashTable<pointFaceHit, Foam::point> &hitMap, Foam::List<Foam::List<Foam::point>> &blockFaceHitPoints);
        void generateCutFace(Foam::List<Foam::point> pts, Foam::HashTable<pointFaceHit, Foam::point> &hitMap, int faceIndex, Foam::face &outputFace, Foam::pointField &outputPoints, Foam::point &targetPoint);

        bool isPointInsideSurface(const Foam::point &p, const Foam::faceList &faces, const Foam::pointField &points);
        int countVerticesInsideBlockBounds(const Foam::face &face1, const Foam::pointField &points);
        int countVerticesInsideSurface(const Foam::pointField &pts, const Foam::faceList &faces, const Foam::pointField &points);

        void addPoints(const Foam::pointField &points, Foam::HashTable<Foam::label, Foam::point> &pointMap, Foam::pointField &mergedPoints);
        void addFace(const Foam::pointField &points, const Foam::face &faceI, Foam::HashTable<Foam::label, Foam::point> &pointMap, Foam::faceList &mergedFaces, Foam::pointField &mergedPoints);

        bool pointOnList(const Foam::List<Foam::point> &pointList, const Foam::point &checkPoint, Foam::point &outputPoint, const double tolerance);
        bool foundOnList(const Foam::List<Foam::point> &pointList, const Foam::point &checkPoint, const double tolerance);
        bool arePointsSame(const Foam::point &point1, const Foam::point &point2, const double tolerance);

        void discretizeFaceToQuads(const Foam::face &inputFace, const Foam::pointField &inputPoints, Foam::pointField &outputPoints, Foam::faceList &outputFaces);
        void decomposeToConvex();
        void triangulateFaces(); // New method declaration to triangulate faces

        ~backgroundBlock();
    };
}

#endif
