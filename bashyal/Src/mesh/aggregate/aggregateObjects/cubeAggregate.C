#include "cubeAggregate.H"
#include "point.H"
#include "scalar.H"
#include "pointField.H"
#include "faceList.H"
#include "List.H"
#include "labelledTri.H" // Includes Foam::labelledTri

#include "Time.H"
#include "IOobject.H"
#include "OFstream.H"
#include "fileName.H"
#include "Random.H"
#include "transform.H"

#include "triSurface.H"
#include "triSurfaceMesh.H"

using namespace Foam;

namespace Bashyal
{
    cubeAggregate::cubeAggregate()
    {
    }

    cubeAggregate::cubeAggregate(float s1, float s2, Foam::label identifier) // to produce size between s1 and s2
        : identifier_(identifier)
    {

        // Create an instance of the Random class with automatic seeding
        Foam::Random randomGen(12);

        Foam::pointField Nodes(8);

        this->s1_ = s1;
        this->s2_ = s2;

        Foam::scalar r = randomGen.sample01<Foam::scalar>();
        Foam::scalar s = s1 + (s2 - s1) * r;

        // Define the 8 vertices of the cube
        Nodes[0] = Foam::point(-s / 2, -s / 2, -s / 2); // Point 0
        Nodes[1] = Foam::point(-s / 2, s / 2, -s / 2);  // Point 1
        Nodes[2] = Foam::point(s / 2, s / 2, -s / 2);   // Point 2
        Nodes[3] = Foam::point(s / 2, -s / 2, -s / 2);  // Point 3
        Nodes[4] = Foam::point(-s / 2, -s / 2, s / 2);  // Point 4
        Nodes[5] = Foam::point(-s / 2, s / 2, s / 2);   // Point 5
        Nodes[6] = Foam::point(s / 2, s / 2, s / 2);    // Point 6
        Nodes[7] = Foam::point(s / 2, -s / 2, s / 2);   // Point 7
        this->localPoints_ = Nodes;

        this->s_ = s;
        this->centroid_ = Foam::point(0, 0, 0);
    }

    void cubeAggregate::createFaces()
    {
        this->createQuadFaces();
    }

    void cubeAggregate::createTriSurface()
    {
        Foam::List<Foam::labelledTri> triangles = this->createTriangularFacesFromPoints(this->points_);
        Foam::triSurface surface(triangles, this->points_);

        Foam::fileName outputFile("/usr/lib/openfoam/openfoam2312/run/debug/cube.stl"); // Second argument is the output file name

        surface.write(outputFile);
        this->surface_ = surface;
    }

    Foam::List<Foam::face> cubeAggregate::createQuadFaces()
    {
        Foam::List<Foam::face> faces(6);
        faces[0] = Foam::face({0, 1, 2, 3}); // Bottom Face
        faces[1] = Foam::face({4, 7, 6, 5}); // Top Face
        faces[2] = Foam::face({0, 4, 5, 1}); // Left Face
        faces[3] = Foam::face({3, 2, 6, 7}); // Right Face
        faces[4] = Foam::face({0, 3, 7, 4}); // Front Face
        faces[5] = Foam::face({1, 5, 6, 2}); // Back Face

        this->faces_ = faces;
        return faces;
    }

    // Function to find the nearest points and form triangles
    Foam::List<Foam::labelledTri> cubeAggregate::createTriangularFacesFromPoints(const Foam::pointField &points)
    {
        // Define the faces of the cube using triangular faces
        Foam::List<Foam::labelledTri> triangles(12);

        triangles[0] = Foam::labelledTri(0, 1, 2); // Triangle 1 of the bottom face
        triangles[1] = Foam::labelledTri(0, 2, 3); // Triangle 2 of the bottom face

        // Top face (4, 5, 6, 7) split into two triangles
        triangles[2] = Foam::labelledTri(4, 7, 6); // Triangle 1 of the top face
        triangles[3] = Foam::labelledTri(4, 6, 5); // Triangle 2 of the top face

        // Front face (0, 1, 5, 4) split into two triangles
        triangles[4] = Foam::labelledTri(0, 4, 5); // Triangle 1 of the front face
        triangles[5] = Foam::labelledTri(0, 5, 1); // Triangle 2 of the front face

        // Back face (2, 3, 7, 6) split into two triangles
        triangles[6] = Foam::labelledTri(3, 2, 6); // Triangle 1 of the back face
        triangles[7] = Foam::labelledTri(3, 6, 7); // Triangle 2 of the back face

        // Left face (0, 3, 7, 4) split into two triangles
        triangles[8] = Foam::labelledTri(0, 3, 7); // Triangle 1 of the left face
        triangles[9] = Foam::labelledTri(0, 7, 4); // Triangle 2 of the left face

        // Right face (1, 2, 6, 5) split into two triangles
        triangles[10] = Foam::labelledTri(1, 5, 6); // Triangle 1 of the right face
        triangles[11] = Foam::labelledTri(1, 6, 2); // Triangle 2 of the right face

        return triangles;
    }

}
