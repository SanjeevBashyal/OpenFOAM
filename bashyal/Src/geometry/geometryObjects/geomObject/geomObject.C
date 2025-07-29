#include "geomObject.H"
#include "pointField.H"
#include "faceList.H"
#include "error.H"
#include "Ostream.H"

namespace Bashyal
{
    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
    // Empty constructor
    geomObject::geomObject()
        : Foam::particleModels::indexedFaceSet() {}

    geomObject::geomObject(const Foam::pointField &points, const Foam::faceList &faces)
        : Foam::particleModels::indexedFaceSet(points, faces)
    {
        if (vertices().empty() || this->faces().empty())
        {
            FatalErrorInFunction
                << "Boundary initialized with empty points or faces"
                << Foam::abort(Foam::FatalError);
        }
    }

    Bashyal::geomObject::geomObject(const Foam::fileName &stlFile)
    {
        // Use triSurface to read the STL file
        Foam::triSurface surface(stlFile);

        // Assign points from the surface
        Foam::pointField pts = surface.points();

        // Get the list of triFaces using the faces() method
        const Foam::List<Foam::labelledTri> &trie = surface.surfFaces();

        // Convert triSurface faces to faceList
        Foam::faceList fcs(trie.size());
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
        if (vertices().empty() || this->faces().empty())
        {
            FatalErrorInFunction
                << "Boundary initialized with empty points or faces from STL file: " << stlFile
                << Foam::abort(Foam::FatalError);
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
            FatalErrorIn("writeVtp") << "Cannot open file " << filename << Foam::exit(Foam::FatalError);
        }

        const Foam::pointField& vertices = this->vertices();
        const Foam::faceList& faces = this->faces();

        // VTK XML Header
        vtpFile << "<?xml version=\"1.0\"?>" << Foam::endl;
        vtpFile << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">" << Foam::endl;
        vtpFile << "  <PolyData>" << Foam::endl;

        // Piece defines the geometry. We have one object, so one piece.
        vtpFile << "    <Piece NumberOfPoints=\"" << vertices.size()
                << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\""
                << faces.size() << "\">" << Foam::endl;

        // 1. Write Points (Vertices)
        vtpFile << "      <Points>" << Foam::endl;
        vtpFile << "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << Foam::endl;
        for (const Foam::point& pt : vertices)
        {
            vtpFile << "          " << pt.x() << " " << pt.y() << " " << pt.z() << Foam::endl;
        }
        vtpFile << "        </DataArray>" << Foam::endl;
        vtpFile << "      </Points>" << Foam::endl;

        // 2. Write Polygons (Faces)
        vtpFile << "      <Polys>" << Foam::endl;
        // a) Connectivity: a flat list of all vertex indices for all faces
        vtpFile << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << Foam::endl;
        for (const Foam::face& f : faces)
        {
            vtpFile << "          ";
            for (const Foam::label vIdx : f)
            {
                vtpFile << vIdx << " ";
            }
            vtpFile << Foam::endl;
        }
        vtpFile << "        </DataArray>" << Foam::endl;

        // b) Offsets: the cumulative count of vertices per face
        vtpFile << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << Foam::endl;
        vtpFile << "          ";
        Foam::label offset = 0;
        for (const Foam::face& f : faces)
        {
            offset += f.size();
            vtpFile << offset << " ";
        }
        vtpFile << Foam::endl;
        vtpFile << "        </DataArray>" << Foam::endl;
        vtpFile << "      </Polys>" << Foam::endl;

        // VTK XML Footer
        vtpFile << "    </Piece>" << Foam::endl;
        vtpFile << "  </PolyData>" << Foam::endl;
        vtpFile << "</VTKFile>" << Foam::endl;
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
                << Foam::abort(Foam::FatalError);
        }
        return Foam::boundBox(); // Return an empty bounding box if vertices() is empty
    }

}
