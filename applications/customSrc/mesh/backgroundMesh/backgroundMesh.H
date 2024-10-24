#ifndef backgroundMesh_H
#define backgroundMesh_H

#include "quickInclude.H"
#include "quickMesh.H"
#include "quickWriter.H"
#include "block.H"
#include "backgroundBlock.H"

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
        faceList globalFaces_;
        labelList globalOwners_;
        labelList globalNeighbours_;
        DynamicList<Foam::face> boundaryFaces_; // Stores boundary faces temporarily

        Foam::HashTable<label, point> pointMap_;
        Foam::HashTable<label, face> faceMap_;
        Foam::HashTable<label, face> boundaryFaceMap_; // For tracking boundary faces

        double s_;
        Foam::Time *runTime_;

        Foam::pointField vertices_;
        Foam::polyMesh *meshPtr_;
        Foam::block *blockPtr_;

        // Foam::List<Foam::face> faceListPtr_;
        // Foam::List<Foam::cell> cellListPtr_;


    public:
        explicit backgroundMesh(Foam::Time *runTime, point meshMin, point meshMax, float resolution);
        Foam::pointField createVertices();
        Vector<int> countBlocksPerAxis() const;
        List<List<List<autoPtr<backgroundBlock>>>> createListPointers();
        
        ~backgroundMesh() = default;

        void developMesh();
        void addPoints(const pointField &blockPoints);
        void addFaces(const pointField &blockPoints, const faceList &blockFaces, const labelList &owners, const labelList &neighbours);

        void createPolyMesh();
        void writeBackgroundMesh(const std::string &meshDir);
        static void writePolyMeshPlain(const std::string &meshDir, const pointField &points, const faceList &faces, const labelList &owners, const labelList &neighbours);
        void writePolyMeshFromOwnerNeighbour(const std::string &meshDir, const pointField &points, const faceList &faces, const labelList &owners, const labelList &neighbours);

        void createBlockMesh();
    };

}

#endif
