#include "backgroundMesh.H"
#include "quickMesh.H"

using namespace Foam;

namespace Bashyal
{
    void backgroundMesh::createBlockMesh()
    {
        Foam::cellShape bshape;
        Foam::pointField vertices;
        Foam::blockEdgeList edges;
        Foam::blockFaceList faces;
        Foam::labelVector density(2, 1, 1);
        Foam::word model("hex");

        Foam::label verticesIndex[8] = {0, 1, 2, 3, 4, 5, 6, 7};

        Foam::UList<Foam::label> verticesLabel(verticesIndex, 8);

        bshape = Foam::cellShape(
            model,
            verticesLabel);

        Foam::block bgBlock(bshape, vertices_, edges, faces, density);
        this->blockPtr_ = &bgBlock;

        Foam::word regionName("region0");
        Foam::word constantLocation = "constant";

        Foam::IOobject *io = new Foam::IOobject(
            regionName,
            constantLocation,
            *(this->runTime_),
            Foam::IOobject::NO_READ,
            Foam::IOobject::AUTO_WRITE);

        Foam::faceListList boundaryFaces(6);
        Foam::wordList boundaryPatchNemes({"YZ-Xmin", "YZ-Xmax", "XZ-Ymin", "XZ-Ymax", "XY-Zmin", "XY-Zmax"});
        Foam::wordList boundaryPatchTypes({"patch", "patch", "patch", "patch", "patch", "patch"});
        Foam::word defaultBoundaryPatchName("boundaryPlane");
        Foam::word defaultBoundaryPatchType("patch");
        Foam::wordList boundaryPatchPhysicalTypes({"wall", "wall", "wall", "wall", "wall", "wall"});

        for (int i = 0; i < 6; i++)
        {
            Foam::faceList singleBoundaryFace(this->blockPtr_->blockPatches_[i].size());
            for (Foam::label j = 0; j < this->blockPtr_->blockPatches_[i].size(); j++)
            {
                singleBoundaryFace[j] = Foam::face(this->blockPtr_->blockPatches_[i][j]);
            }
            boundaryFaces[i] = singleBoundaryFace;
        }

        Foam::polyMesh mesh(
            *io,
            Foam::pointField(bgBlock.points_), // Use a copy of the block points
            bgBlock.shapes(),
            boundaryFaces,
            boundaryPatchNemes,
            boundaryPatchTypes,
            defaultBoundaryPatchName,
            defaultBoundaryPatchType,
            boundaryPatchPhysicalTypes,
            false);

        // Foam::polyMesh mesh(
        //     *io,
        //     Foam::pointField(backgroundBlock.points_), // Use a copy of the block points
        //     Foam::List<Foam::face>(this->faceListPtr_),
        //     Foam::List<Foam::cell>(this->cellListPtr_),
        //     false);

        this->meshPtr_ = &mesh;

        mesh.removeFiles();
        if (!mesh.write())
        {
            FatalErrorInFunction
                << "Failed writing polyMesh."
                << exit(Foam::FatalError);
        }
    }
}
