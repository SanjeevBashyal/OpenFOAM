#include "backgroundMesh.H"
#include "quickMesh.H"

using namespace Foam;

namespace Bashyal
{
    backgroundMesh::backgroundMesh(Foam::Time *runTime, double size)
        : s_(size),
          runTime_(runTime),
          vertices_(createVertices())
    {
        Foam::cellShape bshape;
        Foam::pointField vertices;
        Foam::blockEdgeList edges;
        Foam::blockFaceList faces;
        Foam::labelVector density(10, 10, 10);
        Foam::word model("hex");

        Foam::label verticesIndex[8] = {0, 1, 2, 3, 4, 5, 6, 7};

        Foam::UList<Foam::label> verticesLabel(verticesIndex, 8);

        bshape = Foam::cellShape(
            model,
            verticesLabel);

        Foam::block backgroundBlock(bshape, vertices_, edges, faces, density);
        this->blockPtr_ = &backgroundBlock;

        Foam::word regionName("region0");
        Foam::word constantLocation = "constant";

        // this->createFacesAndCells();
        this->createBoundaryFacesAndPatches();

        Foam::IOobject *io = new Foam::IOobject(
            regionName,
            constantLocation,
            *runTime,
            Foam::IOobject::NO_READ,
            Foam::IOobject::AUTO_WRITE);

        Foam::polyMesh mesh(
            *io,
            Foam::pointField(backgroundBlock.points_), // Use a copy of the block points
            backgroundBlock.shapes(),
            this->boundaryFaces_,
            this->boundaryPatchNemes_,
            this->boundaryPatchTypes_,
            this->defaultBoundaryPatchName_,
            this->defaultBoundaryPatchType_,
            this->boundaryPatchPhysicalTypes_,
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

        Foam::Info << "Here" << Foam::endl;
    }

    // autoPtr<polyMesh> meshPtr =
    //     blocks.mesh(
    //         IOobject(regionName, meshInstance, runTime));

    // // Set the precision of the points data to 10
    // IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    // Info << nl << "Writing polyMesh with "
    //         << mesh.cellZones().size() << " cellZones";

    inline Foam::pointField backgroundMesh::createVertices()
    {
        Foam::pointField Nodes(8);
        Nodes[0] = Foam::point(0, 0, 0);                      // Point 0
        Nodes[1] = Foam::point(this->s_, 0, 0);               // Point 1
        Nodes[2] = Foam::point(this->s_, this->s_, 0);        // Point 2
        Nodes[3] = Foam::point(0, this->s_, 0);               // Point 3
        Nodes[4] = Foam::point(0, 0, this->s_);               // Point 4
        Nodes[5] = Foam::point(this->s_, 0, this->s_);        // Point 5
        Nodes[6] = Foam::point(this->s_, this->s_, this->s_); // Point 6
        Nodes[7] = Foam::point(0, this->s_, this->s_);        // Point 7
        return Nodes;
    }

    void backgroundMesh::createFacesAndCells()
    {
        // Create a HashSet to store unique faces temporarily
        Foam::label count = 0;

        Foam::HashTable<Foam::label, Foam::face> uniqueFaces;
        Foam::List<Foam::cell> cells(this->blockPtr_->blockCells_.size());
        Foam::List<Foam::face> faceList;

        // Loop over each hexCell in the list
        for (Foam::label i = 0; i < this->blockPtr_->blockCells_.size(); ++i)
        {
            Foam::hexCell &hex = this->blockPtr_->blockCells_[i]; // Get the hex cell

            // Get the list of faces from the current hex cell (returns Foam::List<Foam::face>)
            Foam::List<Foam::face> faces = hex.faces();

            Foam::labelList faceLabels(faces.size());

            // Add faces to the HashSet (only unique faces will be kept)

            for (Foam::label j = 0; j < faces.size(); ++j)
            {
                Foam::face &currentFace = faces[j];
                Foam::face reversedFace = currentFace.reverseFace();

                // Check if face already exists in faceMap
                if (uniqueFaces.found(currentFace))
                {
                    faceLabels[j] = uniqueFaces[currentFace];
                    continue;
                }
                else if (uniqueFaces.found(reversedFace))
                {
                    faceLabels[j] = uniqueFaces[reversedFace];
                    continue;
                }
                else
                {
                    uniqueFaces.insert(faces[j], count);
                    faceLabels[j] = count;
                    faceList.append(faces[j]);
                    count = count + 1;
                }
            }
            cells[i] = Foam::cell(faceLabels);
        }

        // Convert the HashSet to a List

        this->faceListPtr_ = faceList;
        this->cellListPtr_ = cells;
    }

    void backgroundMesh::createBoundaryFacesAndPatches()
    {
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

        this->boundaryFaces_ = boundaryFaces;
        this->boundaryPatchNemes_ = boundaryPatchNemes;
        this->boundaryPatchTypes_ = boundaryPatchTypes;
        this->defaultBoundaryPatchName_ = defaultBoundaryPatchName;
        this->defaultBoundaryPatchType_ = defaultBoundaryPatchType;
        this->boundaryPatchPhysicalTypes_ = boundaryPatchPhysicalTypes;
    }
}
