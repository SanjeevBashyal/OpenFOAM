#ifndef backgroundMesh_H
#define backgroundMesh_H

#include "quickInclude.H"
#include "quickMesh.H"
#include "quickWriter.H"
#include "block.H"
#include "backgroundBlock.H"
#include "cubeAggregate.H"
#include "roundAggregate.H"
#include "cubeAggregates.H"

using namespace Foam;

namespace Bashyal
{
    class backgroundMesh
    {
    private:
        /* data */
    public:
        Foam::point meshMin_;
        Foam::point meshMax_;
        Foam::scalar resolution_;

        Vector<int> dim_;
        List<List<List<autoPtr<backgroundBlock>>>> backgroundBlocks_;

        label cellCount_ = 0;
        Foam::pointField globalPoints_;
        faceList globalFaces_;      //all internal faces
        labelList globalOwners_;
        labelList globalNeighbours_;
        
        List<Foam::face> allFaces_; // Stores all faces
        List<bool> boolBoundaryFaces_; // Stores all faces
        List<Foam::word> facePatches_;
        Foam::dictionary boundaryDict_;

        Foam::HashTable<label, point> pointMap_;
        // Foam::HashTable<label, face> faceMap_;
        Foam::HashTable<label, face> faceOwnerMap_; // For tracking boundary faces
        Foam::HashTable<label, face> facePositionMap_; // For tracking face Indices
        Foam::HashTable<faceList, word> stringPtrMap_; // For tracking boundary faces

        double s_;
        Foam::Time *runTime_;

        Foam::pointField vertices_;
        Foam::polyMesh *meshPtr_;
        Foam::block *blockPtr_;

        // Foam::List<Foam::face> faceListPtr_;
        // Foam::List<Foam::cell> cellListPtr_;


    public:
        explicit backgroundMesh(Foam::Time *runTime, point meshMin, point meshMax, float resolution);
        void setBoundaryPatchType(Foam::dictionary &boundaryDict);
        Foam::pointField createVertices();
        Vector<int> countBlocksPerAxis() const;
        List<List<List<autoPtr<backgroundBlock>>>> createListPointers();
        Vector<int> getBlockIndexContainingPoint(point& pt);
        bool contains(const point &pt);
        
        ~backgroundMesh() = default;

        void developMesh();
        void addPoints(const pointField &blockPoints);
        void addFaces(const pointField &blockPoints, const faceList &blockFaces, const wordList &blockPatches, const wordList &stringPtrs);

        void reset();
        point roundPoint(point value);
        scalar roundToSeven(scalar value);

        void intersectCube(cubeAggregate &cubeAgg);
        void intersectRound(roundAggregate &roundAgg);
        void intersectCubes(cubeAggregates &cubeAggs);

        void createPolyMesh();
        void writeBackgroundMesh(const std::string &meshDir);
        static void writePolyMeshPlain(const std::string &meshDir, const pointField &points, const faceList &faces, const labelList &owners, const labelList &neighbours, const Foam::labelList &boundaryFaceSizes, const Foam::wordList &patchNames, const Foam::wordList &patchTypes);
        void writePolyMeshFromOwnerNeighbour(const std::string &meshDir, const pointField &points, const faceList &faces, const labelList &owners, const labelList &neighbours);

        void createBlockMesh();
    };

}

#endif
