#include "backgroundMesh.H"
#include "patchTypes.H"
#include <string>
using namespace Foam;
namespace Bashyal
{
    void backgroundMesh::createPolyMesh()
    {
        // Ensure faces and owners have the same size
        if (globalFaces_.size() != globalOwners_.size())
        {
            FatalErrorInFunction
                << "Mismatch between number of faces and owners!" << exit(FatalError);
        }

        // Construct pointField (vertices)
        pointField meshPoints = globalPoints_; // Already contains all points

        // Construct faceList (faces)
        faceList meshFaces = globalFaces_; // Already contains all faces

        // Construct owner list
        labelList meshOwners = globalOwners_;

        // Construct neighbour list (some faces may not have neighbours)
        labelList meshNeighbours = globalNeighbours_;

        // Handle boundary faces separately: add them to the meshFaces but with no neighbour
        Foam::label boundaryFaceCount = 0;
        for (const auto &faceI : boundaryFaces_)
        {
            meshFaces.append(faceI);                               // Add the boundary face
            meshOwners.append(boundaryOwners_[boundaryFaceCount]); // Boundary faces have no owner
            boundaryFaceCount++;
            // meshNeighbours.append(-1);      // Boundary faces have no neighbour
        }

        // Construct the mesh data (note: this does not write anything to files)
        // polyMesh constructor accepts points, faces, owners, neighbours, and boundary data
        Foam::word regionName("region0");
        Foam::word constantLocation = "constant";

        Foam::IOobject *io = new Foam::IOobject(
            regionName,
            constantLocation,
            *(this->runTime_),
            Foam::IOobject::NO_READ,
            Foam::IOobject::AUTO_WRITE);

        this->meshPtr_ = new Foam::polyMesh(
            *io,
            Foam::pointField(meshPoints), // Use a copy of the block points
            faceList(meshFaces),
            labelList(meshOwners),
            labelList(meshNeighbours),
            false);
    }

    void backgroundMesh::writeBackgroundMesh(const std::string &meshDir)
    {
        // Step 1: Construct points, faces, owners, neighbours lists
        pointField meshPoints = globalPoints_;
        faceList meshFaces = globalFaces_;
        labelList meshOwners = globalOwners_;
        labelList meshNeighbours = globalNeighbours_;
        List<int> meshPatches = boundaryPatches_;

        // Step 2: Initialize containers for sorted boundary faces and patch information
        List<int> patchNames;
        HashTable<Foam::List<std::pair<Foam::face, int>>, int> boundaryFacesMap;

        // Step 3: Sort boundary faces by patch name in globalPatches_
        int i = 0;
        for (const auto &boundaryFace : boundaryFaces_)
        {
            int facePatch = meshPatches[i];

            if (boundaryFacesMap.found(facePatch))
            {
                boundaryFacesMap[facePatch].append(std::make_pair(boundaryFace, boundaryOwners_[i]));
            }
            else
            {
                Foam::List<std::pair<Foam::face, int>> newPatchFaces;
                newPatchFaces.append(std::make_pair(boundaryFace, boundaryOwners_[i]));
                boundaryFacesMap.insert(facePatch, newPatchFaces);
                patchNames.append(facePatch);
            }
            i++;
        }

        // Step 2: Add boundary faces to the meshFaces and mark them with no neighbour
        i = 0;

        labelList patchBoundarySizes;
        patchBoundarySizes.setSize(patchNames.size());

        wordList patchTypes;
        patchTypes.setSize(patchNames.size());
        for (const int patchName : patchNames)
        {
            Foam::List<std::pair<Foam::face, int>> &boundaryFacesI = boundaryFacesMap[patchName];
            patchBoundarySizes[i] = boundaryFacesI.size();
            patchTypes[i] = boundaryDict_.getOrDefault<Foam::word>(std::to_string(patchName), "patch");
            i++;

            for (const auto &p : boundaryFacesI)
            {
                meshFaces.append(p.first);
                meshOwners.append(p.second); // Boundary face, no owner
                // meshNeighbours.append(-1); // Boundary face, no neighbour
            }
        }

        // Step 3: Call the writePolyMeshPlain function to write the mesh data to files
        writePolyMeshPlain(meshDir, meshPoints, meshFaces, meshOwners, meshNeighbours, patchBoundarySizes, patchNames, patchTypes);
    }

    // Static method to write polyMesh data
    void backgroundMesh::writePolyMeshPlain(
        const std::string &meshDir,
        const pointField &points,
        const faceList &faces,
        const labelList &owners,
        const labelList &neighbours,
        const Foam::labelList &boundaryFaceSizes,
        const Foam::List<int> &patchInts,
        const Foam::wordList &patchTypes)
    {
        mkDir(meshDir.c_str(), 0777); // Create the mesh directory if it doesn't exist

        // Write points file
        std::string outputDir = "constant/polyMesh";
        std::string pointsFileName = meshDir + "/points";
        OFstream pointsFile(pointsFileName);
        pointsFile << "FoamFile\n"
                   << "{\n"
                   << "    version     2.0;\n"
                   << "    format      ascii;\n"
                   << "    class       vectorField;\n"
                   << "    location    " << outputDir << ";\n"
                   << "    object      points;\n"
                   << "}\n\n";
        pointsFile << points.size() << "\n(\n";
        for (const auto &pt : points)
        {
            pointsFile << pt << "\n";
        }
        pointsFile << ")\n";

        // Write faces file
        std::string facesFileName = meshDir + "/faces";
        OFstream facesFile(facesFileName);
        facesFile << "FoamFile\n"
                  << "{\n"
                  << "    version     2.0;\n"
                  << "    format      ascii;\n"
                  << "    class       faceList;\n"
                  << "    location    " << outputDir << ";\n"
                  << "    object      faces;\n"
                  << "}\n\n";
        facesFile << faces.size() << "\n(\n";
        for (const auto &face : faces)
        {
            facesFile << face.size() << "(";
            for (const auto &vert : face)
            {
                facesFile << vert << " ";
            }
            facesFile << ")\n";
        }
        facesFile << ")\n";

        // Write owner file
        std::string ownerFileName = meshDir + "/owner";
        OFstream ownerFile(ownerFileName);
        ownerFile << "FoamFile\n"
                  << "{\n"
                  << "    version     2.0;\n"
                  << "    format      ascii;\n"
                  << "    class       labelList;\n"
                  << "    location    " << outputDir << ";\n"
                  << "    object      owner;\n"
                  << "}\n\n";
        ownerFile << owners.size() << "\n(\n";
        for (const auto &owner : owners)
        {
            ownerFile << owner << "\n";
        }
        ownerFile << ")\n";

        // Write neighbour file
        std::string neighbourFileName = meshDir + "/neighbour";
        OFstream neighbourFile(neighbourFileName);
        neighbourFile << "FoamFile\n"
                      << "{\n"
                      << "    version     2.0;\n"
                      << "    format      ascii;\n"
                      << "    class       labelList;\n"
                      << "    location    " << outputDir << ";\n"
                      << "    object      neighbour;\n"
                      << "}\n\n";
        neighbourFile << neighbours.size() << "\n(\n";
        for (const auto &neighbour : neighbours)
        {
            neighbourFile << neighbour << "\n";
        }
        neighbourFile << ")\n";

        // Write boundary file
        OFstream boundaryFile(meshDir + "/boundary");
        boundaryFile << "FoamFile\n"
                     << "{\n"
                     << "    version     2.0;\n"
                     << "    format      ascii;\n"
                     << "    class       polyBoundaryMesh;\n"
                     << "    location    " << outputDir << ";\n"
                     << "    object      boundary;\n"
                     << "}\n\n";
        boundaryFile << patchInts.size() << "\n(\n";

        // Start face index for each patch
        label startFace = neighbours.size();

        std::map<int, Foam::word> patchNamesMap = {
            {patchType::YZ_Xmin, "Xmin"},
            {patchType::YZ_Xmax, "Xmax"},
            {patchType::XZ_Ymin, "Ymin"},
            {patchType::XZ_Ymax, "Ymax"},
            {patchType::XY_Zmin, "Zmin"},
            {patchType::XY_Zmax, "Zmax"},
            {1000, "aggregate"}};

        for (Foam::label i = 0; i < patchInts.size(); ++i)
        {
            boundaryFile << patchNamesMap[patchInts[i]] << "\n{\n";
            boundaryFile << "    type " << patchTypes[i] << ";\n";
            boundaryFile << "    nFaces " << boundaryFaceSizes[i] << ";\n";
            boundaryFile << "    startFace " << startFace << ";\n";
            boundaryFile << "}\n";
            startFace += boundaryFaceSizes[i];
        }
        boundaryFile << ")\n";

        Info << "Mesh written to " << meshDir << nl;
    }

    void backgroundMesh::writePolyMeshFromOwnerNeighbour(
        const std::string &meshDir,
        const pointField &points,
        const faceList &faces,
        const labelList &owners,
        const labelList &neighbours)
    {
        Foam::word regionName("region0");
        Foam::word constantLocation = "constant";

        Foam::IOobject *io = new Foam::IOobject(
            regionName,
            constantLocation,
            *(this->runTime_),
            Foam::IOobject::NO_READ,
            Foam::IOobject::AUTO_WRITE);

        Foam::polyMesh mesh(
            *io,
            Foam::pointField(points), // Use a copy of the block points
            faceList(faces),
            labelList(owners),
            labelList(neighbours),
            false);

        mesh.write();
    }

}
