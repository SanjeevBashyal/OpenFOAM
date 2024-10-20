#include "backgroundMesh.H"

namespace Bashyal
{
    // Static method to write polyMesh data
    void backgroundMesh::writePolyMeshPlain(
        const std::string &meshDir,
        const pointField &points,
        const faceList &faces,
        const labelList &owners,
        const labelList &neighbours)
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
        boundaryFile << "1\n(\n";
        boundaryFile << "cube\n{\n";
        boundaryFile << "    type patch;\n";
        boundaryFile << "    nFaces " << faces.size() << ";\n";
        boundaryFile << "    startFace 0;\n";
        boundaryFile << "}\n)\n";

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
