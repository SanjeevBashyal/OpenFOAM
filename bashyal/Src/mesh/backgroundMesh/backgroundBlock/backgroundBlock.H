#ifndef backgroundBlock_H
#define backgroundBlock_H

#include "debugClass.H"
#include "quickMesh.H"
#include "boundary.H"
#include "aggregate.H"
#include <CGAL/Nef_polyhedron_3.h>

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
        Foam::label nboundaries_;
        Foam::label globalNCells_ = 0;
        Foam::label blockID_;

        bool dead_ = false;
        bool edited_ = false;
        bool multiple_ = false;

        Foam::pointField points_; // Stores the vertices of the cube
        Foam::faceList faces_;    // Stores the six faces of the cube
        Foam::List<int> patches_;

        Foam::labelList owners_;     // Stores owner cells
        Foam::labelList neighbours_; // Stores neighbor cells

        CGAL::Nef_polyhedron_3<CGAL::Epeck> nef_;
        CGAL::Nef_polyhedron_3<CGAL::Epeck> nefBase_;
        // Foam::faceList triFaces_; // New attribute to store triangulated faces
        Foam::List<const boundary*> intersectedBoundaries_;

    private:
        /* data */
    public:
        // backgroundBlock.C
        backgroundBlock(backgroundMesh *ref, const Foam::Vector<int> identity, const Foam::boundBox &bounds, Foam::label blockID);
        void generateCubeGeometry();
        void reset();
        bool contains(const Foam::point &pt) const;
        void triangulateFaces(const Foam::pointField &points, const Foam::faceList &faces, Foam::pointField &outPoints, Foam::faceList &outFaces);
        ~backgroundBlock();

        Foam::label getGlobalNcells() { return globalNCells_; }
        Foam::label getBlockID() { return blockID_; }
        Foam::pointField getPoints() { return points_; }
        Foam::faceList getFaces() { return faces_; }
        Foam::List<int> getPatches() { return patches_; }
        Foam::labelList getOwners() { return owners_; }
        Foam::labelList getNeighbours() { return neighbours_; }

        void develop();
        void intersectBoundary(const boundary &domainBoundary, bool keepInside);
        void subtractAggregate(const aggregate &agg);

        // backgroundBlockCGALIntersection.C
        void intersectClosedSurfaceCGAL(const Foam::faceList &faces, const Foam::pointField &points, unsigned int identifier);
        void mapNewFacesToBoundaries(const Foam::pointField &newPoints, const Foam::faceList &newFaces, const Foam::labelList &newNeighbours, const Foam::pointField &pointsA, const Foam::faceList &facesA, const Foam::List<int> &patchesA, const Foam::pointField &pointsB, const Foam::faceList &facesB, const Foam::List<int> &patchesB, Foam::List<int> &newPatches);
        void reorderToUpperTriangularNeighboursOnly(const Foam::faceList &faces, const Foam::labelList &owners, const Foam::labelList &neighbours, const Foam::List<int> &patches, Foam::faceList &outFaces, Foam::labelList &outOwners, Foam::labelList &outNeighbours, Foam::List<int> &outPatches);

        // backgroundBlockSearch.C
        bool isPointInsideSurface(const Foam::point &p, const Foam::faceList &faces, const Foam::pointField &points);
        bool pointOnList(const Foam::List<Foam::point> &pointList, const Foam::point &checkPoint, Foam::point &outputPoint, const double tolerance);
        bool arePointsSame(const Foam::point &point1, const Foam::point &point2, const double tolerance);

        void generateNefPolyhedron();

        void write() const; // To add later in backgroundBlockWriter
    };
}

#endif
