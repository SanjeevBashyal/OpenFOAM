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
    cubeAggregate::cubeAggregate(float s1, float s2) // to produce size between s1 and s2
    :
    xRotation_(0),
    yRotation_(0),
    zRotation_(0)
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
        Nodes[1] = Foam::point(s / 2, -s / 2, -s / 2);  // Point 1
        Nodes[2] = Foam::point(s / 2, s / 2, -s / 2);   // Point 2
        Nodes[3] = Foam::point(-s / 2, s / 2, -s / 2);  // Point 3
        Nodes[4] = Foam::point(-s / 2, -s / 2, s / 2);  // Point 4
        Nodes[5] = Foam::point(s / 2, -s / 2, s / 2);   // Point 5
        Nodes[6] = Foam::point(s / 2, s / 2, s / 2);    // Point 6
        Nodes[7] = Foam::point(-s / 2, s / 2, s / 2);   // Point 7
        this->localPoints_ = Nodes;

        this->s_ = s;
        this->centroid_ = Foam::point(0, 0, 0);
    }

    void cubeAggregate::createSurface()
    {
        Foam::List<Foam::labelledTri> triangles = this->createFacesFromPoints(this->globalPoints_);
        Foam::triSurface surface(triangles, this->globalPoints_);

        Foam::fileName outputFile("/usr/lib/openfoam/openfoam2312/run/aggregate/cube.stl"); // Second argument is the output file name

        surface.write(outputFile);
        this->surface_ = surface;
    }

    // Function to find the nearest points and form triangles
    Foam::List<Foam::labelledTri> cubeAggregate::createFacesFromPoints(const Foam::pointField &points)
    {
        // Define the faces of the cube using triangular faces
        Foam::List<Foam::labelledTri> triangles(12);

        triangles[0] = Foam::labelledTri(0, 1, 2); // Triangle 1 of the bottom face
        triangles[1] = Foam::labelledTri(0, 2, 3); // Triangle 2 of the bottom face

        // Top face (4, 5, 6, 7) split into two triangles
        triangles[2] = Foam::labelledTri(4, 5, 6); // Triangle 1 of the top face
        triangles[3] = Foam::labelledTri(4, 6, 7); // Triangle 2 of the top face

        // Front face (0, 1, 5, 4) split into two triangles
        triangles[4] = Foam::labelledTri(0, 1, 5); // Triangle 1 of the front face
        triangles[5] = Foam::labelledTri(0, 5, 4); // Triangle 2 of the front face

        // Back face (2, 3, 7, 6) split into two triangles
        triangles[6] = Foam::labelledTri(2, 3, 7); // Triangle 1 of the back face
        triangles[7] = Foam::labelledTri(2, 7, 6); // Triangle 2 of the back face

        // Left face (0, 3, 7, 4) split into two triangles
        triangles[8] = Foam::labelledTri(0, 3, 7); // Triangle 1 of the left face
        triangles[9] = Foam::labelledTri(0, 7, 4); // Triangle 2 of the left face

        // Right face (1, 2, 6, 5) split into two triangles
        triangles[10] = Foam::labelledTri(1, 2, 6); // Triangle 1 of the right face
        triangles[11] = Foam::labelledTri(1, 6, 5); // Triangle 2 of the right face

        return triangles;
    }

    void cubeAggregate::translate(Foam::vector translationVector)
    {
        this->centroid_ += translationVector;
    }

    pointField cubeAggregate::translatePoints(pointField points, Foam::vector translationVector)
    {
        for (Foam::label i = 0; i < points.size(); i++)
        {
            points[i] += translationVector;
        }
        return points;
    }

    void cubeAggregate::rotate(scalar alpha, scalar beta, scalar gamma)
    {
        this->xRotation_ += alpha;
        this->yRotation_ += beta;
        this->zRotation_ += gamma;
    }

    pointField cubeAggregate::rotatePoints(Foam::tensor rotationMatrix)
    {
        return Foam::transform(rotationMatrix, this->localPoints_);
    }

    tensor cubeAggregate::rotationMatrixFromAngles()
    {
        // Convert angles from degrees to radians (if necessary)
        scalar xRad = this->xRotation_ * M_PI / 180.0;
        scalar yRad = this->yRotation_ * M_PI / 180.0;
        scalar zRad = this->zRotation_ * M_PI / 180.0;

        // Rotation matrix for rotation around X-axis
        tensor Rx(
            1, 0, 0,
            0, std::cos(xRad), -std::sin(xRad),
            0, std::sin(xRad), std::cos(xRad));

        // Rotation matrix for rotation around Y-axis
        tensor Ry(
            std::cos(yRad), 0, std::sin(yRad),
            0, 1, 0,
            -std::sin(yRad), 0, std::cos(yRad));

        // Rotation matrix for rotation around Z-axis
        tensor Rz(
            std::cos(zRad), -std::sin(zRad), 0,
            std::sin(zRad), std::cos(zRad), 0,
            0, 0, 1);

        // Combined rotation matrix: R = Rz * Ry * Rx
        return Rz * Ry * Rx;
    }

    void cubeAggregate::locate()
    {
        pointField rotatedPoints = this->rotatePoints(this->rotationMatrixFromAngles());
        vector translationVector(this->centroid_[0], this->centroid_[1], this->centroid_[2]);
        this->globalPoints_ = this->translatePoints(rotatedPoints, translationVector);
    }

}
