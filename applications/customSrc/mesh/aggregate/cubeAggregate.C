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
        : xRotation_(0),
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

    void cubeAggregate::createSurface()
    {
        Foam::List<Foam::labelledTri> triangles = this->createTriangularFacesFromPoints(this->globalPoints_);
        Foam::triSurface surface(triangles, this->globalPoints_);

        Foam::fileName outputFile("/usr/lib/openfoam/openfoam2312/run/test/cube.stl"); // Second argument is the output file name

        surface.write(outputFile);
        this->surface_ = surface;
        this->createQuadFaces();
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

    boundBox cubeAggregate::getBoundBox()
    {
        this->boundBox_ = boundBox(this->globalPoints_);
        point min = this->floorPoint(this->boundBox_.min());
        point max = this->ceilPoint(this->boundBox_.max());
        // this->roundedBoundBox_ = boundBox(min, max);
        return boundBox(min, max);
    }

    scalar cubeAggregate::roundToRequiredDecimal(scalar value)
    {
        return std::round(value * std::pow(10.0, this->backgroundFinenessIndex_)) / std::pow(10.0, this->backgroundFinenessIndex_);
    }

    scalar cubeAggregate::floorToRequiredDecimal(scalar value)
    {
        return std::floor(value * std::pow(10.0, this->backgroundFinenessIndex_)) / std::pow(10.0, this->backgroundFinenessIndex_);
    }

    scalar cubeAggregate::ceilToRequiredDecimal(scalar value)
    {
        return std::ceil(value * std::pow(10.0, this->backgroundFinenessIndex_)) / std::pow(10.0, this->backgroundFinenessIndex_);
    }

    point cubeAggregate::roundPoint(point value)
    {
        return point(this->roundToRequiredDecimal(value[0]), this->roundToRequiredDecimal(value[1]), this->roundToRequiredDecimal(value[2]));
    }

    point cubeAggregate::floorPoint(point value)
    {
        return point(this->floorToRequiredDecimal(value[0]), this->floorToRequiredDecimal(value[1]), this->floorToRequiredDecimal(value[2]));
    }

    point cubeAggregate::ceilPoint(point value)
    {
        return point(this->ceilToRequiredDecimal(value[0]), this->ceilToRequiredDecimal(value[1]), this->ceilToRequiredDecimal(value[2]));
    }

    void cubeAggregate::hit(scalar index)
    {
        point min = this->roundedBoundBox_.min();
        point max = this->roundedBoundBox_.max();
        int count = 0;

        for (face face_0 : this->faces_)
        {
            pointField facePoints = face_0.points(this->globalPoints_);
            for (scalar i = min[0]; i < max[0]; i = i + 0.1)
            {
                for (scalar j = min[1]; j < max[1]; j = j + 0.1)
                {
                    point p(i, j, min[2]);
                    vector v(0, 0, 1);
                    pointHit hit_0 = face_0.ray(p, v, this->globalPoints_);
                    if (hit_0.hit())
                    {
                        facePoints.append(hit_0.point());
                    }
                }

                for (scalar k = min[2]; k < max[2]; k = k + 0.1)
                {
                    point p(i, min[1], k);
                    vector v(0, 1, 0);
                    pointHit hit_0 = face_0.ray(p, v, this->globalPoints_);
                    if (hit_0.hit())
                    {
                        facePoints.append(hit_0.point());
                    }
                }
            }
            for (scalar j = min[1]; j < max[1]; j = j + 0.1)
            {
                for (scalar k = min[2]; k < max[2]; k = k + 0.1)
                {
                    point p(min[0], j, k);
                    vector v(1, 0, 0);
                    pointHit hit_0 = face_0.ray(p, v, this->globalPoints_);
                    if (hit_0.hit())
                    {
                        facePoints.append(hit_0.point());
                    }
                }
            }
            this->facePoints_.append(facePoints);
            count++;
        }
    }

}
