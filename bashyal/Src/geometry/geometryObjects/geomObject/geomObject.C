#include "geomObject.H"

namespace Bashyal
{
    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
    // Empty constructor
    geomObject::geomObject()
        : Foam::particleModels::indexedFaceSet() {}

    geomObject::geomObject(const Foam::pointField &points, const Foam::faceList &faces)
        : Foam::particleModels::indexedFaceSet(points, faces)
    {
        if (vertices().empty() || faces().empty())
        {
            FatalErrorInFunction
                << "Boundary initialized with empty points or faces"
                << abort(Foam::FatalError);
        }
    }

    Bashyal::geomObject::geomObject(const Foam::fileName &stlFile)
    {
        // Use triSurface to read the STL file
        Foam::triSurface surface(stlFile);

        // Assign points from the surface
        pointField pts = surface.points();

        // Get the list of triFaces using the faces() method
        const Foam::List<Foam::labelledTri> &trie = surface.surfFaces();

        // Convert triSurface faces to faceList
        faceList fcs(trie.size());
        forAll(trie, i)
        {
            const Foam::triFace &tri = trie[i];
            fcs[i] = Foam::face(3);
            fcs[i][0] = tri[0];
            fcs[i][1] = tri[1];
            fcs[i][2] = tri[2];
        }
        // Assign to base class
        *static_cast<Foam::particleModels::indexedFaceSet*>(this) = Foam::particleModels::indexedFaceSet(pts, fcs);

        // Basic validation
        if (vertices().empty() || faces().empty())
        {
            FatalErrorInFunction
                << "Boundary initialized with empty points or faces from STL file: " << stlFile
                << abort(Foam::FatalError);
        }
    }

    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    geomObject::~geomObject()
    {
    }

    // * * * * * * * * * * * * * * * * Methods  * * * * * * * * * * * * * * * //

    void geomObject::writeToSTL(const Foam::fileName &stlFile) const
    {
        faceOperations::writeToSTL(vertices(), faces(), stlFile);
    }

    void Bashyal::geomObject::writeVtp(const Foam::fileName &filename) const
    {
        Foam::OFstream vtpFile(filename);
        if (!vtpFile.good())
        {
            FatalErrorIn("writeVtp") << "Cannot open file " << filename << exit(FatalError);
        }

        const Foam::pointField& vertices = this->vertices();
        const Foam::faceList& faces = this->faces();

        // VTK XML Header
        vtpFile << "<?xml version=\"1.0\"?>" << endl;
        vtpFile << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
        vtpFile << "  <PolyData>" << endl;

        // Piece defines the geometry. We have one object, so one piece.
        vtpFile << "    <Piece NumberOfPoints=\"" << vertices.size()
                << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\""
                << faces.size() << "\">" << endl;

        // 1. Write Points (Vertices)
        vtpFile << "      <Points>" << endl;
        vtpFile << "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
        for (const Foam::point& pt : vertices)
        {
            vtpFile << "          " << pt.x() << " " << pt.y() << " " << pt.z() << endl;
        }
        vtpFile << "        </DataArray>" << endl;
        vtpFile << "      </Points>" << endl;

        // 2. Write Polygons (Faces)
        vtpFile << "      <Polys>" << endl;
        // a) Connectivity: a flat list of all vertex indices for all faces
        vtpFile << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;
        for (const Foam::face& f : faces)
        {
            vtpFile << "          ";
            for (const Foam::label vIdx : f)
            {
                vtpFile << vIdx << " ";
            }
            vtpFile << endl;
        }
        vtpFile << "        </DataArray>" << endl;

        // b) Offsets: the cumulative count of vertices per face
        vtpFile << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
        vtpFile << "          ";
        Foam::label offset = 0;
        for (const Foam::face& f : faces)
        {
            offset += f.size();
            vtpFile << offset << " ";
        }
        vtpFile << endl;
        vtpFile << "        </DataArray>" << endl;
        vtpFile << "      </Polys>" << endl;

        // VTK XML Footer
        vtpFile << "    </Piece>" << endl;
        vtpFile << "  </PolyData>" << endl;
        vtpFile << "</VTKFile>" << endl;
    }

    Foam::boundBox geomObject::createBoundBox()
    {
        if (vertices().size() > 0)
        {
            return Foam::boundBox(vertices());
        }
        else
        {
            FatalErrorInFunction
                << "Cannot create bounding box from empty vertices()"
                << abort(Foam::FatalError);
        }
        return Foam::boundBox(); // Return an empty bounding box if vertices() is empty
    }

}
