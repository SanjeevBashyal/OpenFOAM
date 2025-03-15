#ifndef backgroundMesh_H
#define backgroundMesh_H

#include "debugClass.H"
#include "quickInclude.H"
#include "quickMesh.H"
#include "quickWriter.H"
#include "block.H"
#include "backgroundBlock.H"
#include "cubeAggregate.H"
#include "roundAggregate.H"
#include "cubeAggregates.H"

namespace Bashyal
{
    class backgroundMesh : public debugClass
    {
    private:
        /* data */
    public:
        Foam::point meshMin_;
        Foam::point meshMax_;
        Foam::scalar resolution_;

        Foam::Vector<int> dim_;
        Foam::List<Foam::List<Foam::List<Foam::autoPtr<backgroundBlock>>>> backgroundBlocks_;

        Foam::label cellCount_ = 0;
        Foam::pointField globalPoints_;
        Foam::faceList globalFaces_; // all internal faces
        Foam::labelList globalOwners_;
        Foam::labelList globalNeighbours_;
        // Foam::List<int> globalPatches_;

        Foam::List<Foam::face> boundaryFaces_;    // Stores all boundary faces
        Foam::labelList boundaryOwners_;    // Stores all boundary faces
        Foam::List<int> boundaryPatches_;
        Foam::List<bool> boolBoundaryFaces_; // Stores boundaryFace flag
        
        Foam::dictionary boundaryDict_;

        Foam::HashTable<Foam::label, Foam::point> pointMap_;
        // Foam::HashTable<Foam::label, Foam::face> faceMap_;
        // Foam::HashTable<Foam::label, Foam::face> faceOwnerMap_;    // For tracking boundary faces
        // Foam::HashTable<Foam::label, Foam::face> facePositionMap_; // For tracking face Indices
        // Foam::HashTable<Foam::faceList, Foam::word> stringPtrMap_; // For tracking boundary faces

        double s_;
        Foam::Time *runTime_;

        Foam::pointField vertices_;
        Foam::polyMesh *meshPtr_;
        Foam::block *blockPtr_;

        // Foam::List<Foam::face> faceListPtr_;
        // Foam::List<Foam::cell> cellListPtr_;

    public:
        explicit backgroundMesh(Foam::Time *runTime, Foam::point meshMin, Foam::point meshMax, float resolution);
        void setBoundaryPatchType(Foam::dictionary &boundaryDict);
        Foam::pointField createVertices();
        Foam::Vector<int> countBlocksPerAxis() const;
        Foam::List<Foam::List<Foam::List<Foam::autoPtr<backgroundBlock>>>> createListPointers();
        Foam::Vector<int> getBlockIndexContainingPoint(Foam::point &pt);
        bool contains(const Foam::point &pt);

        ~backgroundMesh() = default;

        void auditMesh();
        void developMesh();
        void addPoints(const Foam::pointField &blockPoints);
        void addFaces(Foam::label globalNcells, const Foam::Vector<int> &identity, const Foam::pointField &blockPoints, const Foam::faceList &blockFaces, const Foam::labelList &owners, const Foam::labelList &neighbours, const Foam::List<int> &blockPatches);

        void reset();
        void resetBlocks();

        Foam::point roundPoint(Foam::point value);
        Foam::scalar roundToSeven(Foam::scalar value);

        void intersectCube(cubeAggregate &cubeAgg);
        void intersectRound(roundAggregate &roundAgg);
        void intersectCubes(cubeAggregates &cubeAggs);

        void createPolyMesh();
        void writeBackgroundMesh(const std::string &meshDir);
        static void writePolyMeshPlain(const std::string &meshDir, const Foam::pointField &points, const Foam::faceList &faces, const Foam::labelList &owners, const Foam::labelList &neighbours, const Foam::labelList &boundaryFaceSizes, const Foam::List<int> &patchNames, const Foam::wordList &patchTypes);
        void writePolyMeshFromOwnerNeighbour(const std::string &meshDir, const Foam::pointField &points, const Foam::faceList &faces, const Foam::labelList &owners, const Foam::labelList &neighbours);

        bool isFaceMeshBoundary(const Foam::Vector<int> &identity, int patchType);

        void createBlockMesh();
    };
}

#endif
