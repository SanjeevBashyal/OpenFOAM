#include "debugClass.H"
#include <map>
#include <sstream>
#include <iomanip>

namespace Bashyal
{
    void debugClass::writePolyhedronAsPolyMesh(Polyhedron& polyhedron, const std::string& folderName)
    {
        // Define the base path
        const std::string basePath = "/usr/lib/openfoam/openfoam2312/run/debug/";

        // Construct the full path for the polyMesh directory
        std::string meshDir = basePath + folderName;

        // Define the path for the .foam file (e.g., baseName_0.foam)
        std::string foamFilePath = meshDir + "/" + folderName + ".foam";

        // Create a blank .foam file in the root folder
        Foam::OFstream foamFile(foamFilePath);
        
        // 1. Extract points from polyhedron vertices
        std::map<const Polyhedron::Vertex*, Foam::label> vertexIndexMap;
        Foam::pointField points(polyhedron.size_of_vertices());
        Foam::label idx = 0;
        for (auto vit = polyhedron.vertices_begin(); vit != polyhedron.vertices_end(); ++vit)
        {
            const auto& point = vit->point();
            points[idx] = Foam::vector(
                CGAL::to_double(point.x()),
                CGAL::to_double(point.y()),
                CGAL::to_double(point.z())
            );
            vertexIndexMap[&(*vit)] = idx;
            ++idx;
        }

        // 2. Extract faces from polyhedron facets
        Foam::faceList faces(polyhedron.size_of_facets());
        Foam::label faceIdx = 0;
        for (auto fit = polyhedron.facets_begin(); fit != polyhedron.facets_end(); ++fit)
        {
            Foam::labelList faceVertices;
            auto circ = fit->facet_begin();
            do
            {
                faceVertices.append(vertexIndexMap[&(*circ->vertex())]);
                ++circ;
            } while (circ != fit->facet_begin());
            faces[faceIdx] = Foam::face(faceVertices);
            ++faceIdx;
        }

        // 3. Set owners: all faces owned by cell 0 (single cell)
        Foam::labelList owners(faces.size(), 0);

        // 4. Set neighbours: empty since all faces are boundary faces
        Foam::labelList neighbours(0);

        // 5. Define boundary: one patch containing all faces
        Foam::word patchName = "aggregate";
        Foam::word patchType = "patch";
        Foam::label nFaces = faces.size();
        Foam::label startFace = 0;

        // 6. Write polyMesh files directly in the specified directory
        std::string outputDir = meshDir + "/constant/polyMesh";
        Foam::mkDir(outputDir.c_str(), 0777); // Create directory if it doesn't exist

        // Write points file
        {
            Foam::OFstream outFile(outputDir + "/points");
            outFile << "FoamFile\n"
                    << "{\n"
                    << "    version     2.0;\n"
                    << "    format      ascii;\n"
                    << "    class       vectorField;\n"
                    << "    location    \"constant/polyMesh\";\n"
                    << "    object      points;\n"
                    << "}\n\n";
            outFile << points.size() << "\n(\n";
            for (const auto& pt : points)
            {
                outFile << "(" << pt.x() << " " << pt.y() << " " << pt.z() << ")\n";
            }
            outFile << ")\n";
        }

        // Write faces file
        {
            Foam::OFstream outFile(outputDir + "/faces");
            outFile << "FoamFile\n"
                    << "{\n"
                    << "    version     2.0;\n"
                    << "    format      ascii;\n"
                    << "    class       faceList;\n"
                    << "    location    \"constant/polyMesh\";\n"
                    << "    object      faces;\n"
                    << "}\n\n";
            outFile << faces.size() << "\n(\n";
            for (const auto& face : faces)
            {
                outFile << face.size() << "(";
                for (Foam::label vert : face)
                {
                    outFile << vert << " ";
                }
                outFile << ")\n";
            }
            outFile << ")\n";
        }

        // Write owner file
        {
            Foam::OFstream outFile(outputDir + "/owner");
            outFile << "FoamFile\n"
                    << "{\n"
                    << "    version     2.0;\n"
                    << "    format      ascii;\n"
                    << "    class       labelList;\n"
                    << "    location    \"constant/polyMesh\";\n"
                    << "    object      owner;\n"
                    << "}\n\n";
            outFile << owners.size() << "\n(\n";
            for (const auto& label : owners)
            {
                outFile << label << "\n";
            }
            outFile << ")\n";
        }

        // Write neighbour file
        {
            Foam::OFstream outFile(outputDir + "/neighbour");
            outFile << "FoamFile\n"
                    << "{\n"
                    << "    version     2.0;\n"
                    << "    format      ascii;\n"
                    << "    class       labelList;\n"
                    << "    location    \"constant/polyMesh\";\n"
                    << "    object      neighbour;\n"
                    << "}\n\n";
            outFile << neighbours.size() << "\n(\n";
            for (const auto& label : neighbours)
            {
                outFile << label << "\n";
            }
            outFile << ")\n";
        }

        // Write boundary file
        {
            Foam::OFstream outFile(outputDir + "/boundary");
            outFile << "FoamFile\n"
                    << "{\n"
                    << "    version     2.0;\n"
                    << "    format      ascii;\n"
                    << "    class       polyBoundaryMesh;\n"
                    << "    location    \"constant/polyMesh\";\n"
                    << "    object      boundary;\n"
                    << "}\n\n";
            outFile << "1\n(\n"; // One boundary patch
            outFile << patchName << "\n{\n";
            outFile << "    type " << patchType << ";\n";
            outFile << "    nFaces " << nFaces << ";\n";
            outFile << "    startFace " << startFace << ";\n";
            outFile << "}\n";
            outFile << ")\n";
        }

        Foam::Info << "PolyMesh written for polyhedron at " << outputDir << Foam::endl;
    }

    void debugClass::writePoly(std::vector<Polyhedron>& polyhedrons)
    {
        std::string baseName = "polyhedron";
        // Iterate over each polyhedron in the vector
        for (size_t i = 0; i < polyhedrons.size(); ++i)
        {
            // Generate a unique directory name using baseName and index
            std::string folderName = baseName + "_" + std::to_string(i);
            
            // Call the existing method to write polyMesh files for this polyhedron
            writePolyhedronAsPolyMesh(polyhedrons[i], folderName);
        }
    }
}
