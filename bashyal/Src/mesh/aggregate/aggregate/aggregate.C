#include "aggregate.H"
#include "point.H"
#include "scalar.H"
#include "pointField.H"
#include "faceList.H"
#include "List.H"

#include "Time.H"
#include "fileName.H"
#include "Random.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"

using namespace Foam;
namespace Bashyal
{
    aggregate::aggregate()
    {
    }

    Foam::point aggregate::sphericalToCartesian(Foam::scalar r, Foam::scalar theta, Foam::scalar phi)
    {
        // Using OpenFOAM's point structure for Cartesian coordinates
        Foam::scalar x = r * std::sin(theta) * std::cos(phi);
        Foam::scalar y = r * std::sin(theta) * std::sin(phi);
        Foam::scalar z = r * std::cos(theta);

        // Return the result as a point
        return Foam::point(x, y, z);
    }

    Foam::scalar aggregate::distance(const Foam::point &p1, const Foam::point &p2)
    {
        return mag(p2 - p1);
    }

    void aggregate::translate(Foam::vector translationVector)
    {
        this->centroid_ += translationVector;
    }

    void aggregate::position(Foam::vector positionThis)
    {
        this->centroid_ = positionThis;
    }

    pointField aggregate::translatePoints(pointField points, Foam::vector translationVector)
    {
        for (Foam::label i = 0; i < points.size(); i++)
        {
            points[i] += translationVector;
        }
        return points;
    }

    void aggregate::rotate(scalar alpha, scalar beta, scalar gamma)
    {
        this->xRotation_ += alpha;
        this->yRotation_ += beta;
        this->zRotation_ += gamma;
    }

    pointField aggregate::rotatePoints(Foam::tensor rotationMatrix)
    {
        return Foam::transform(rotationMatrix, this->localPoints_);
    }

    tensor aggregate::rotationMatrixFromAngles()
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

    void aggregate::locate()
    {
        pointField rotatedPoints = this->rotatePoints(this->rotationMatrixFromAngles());
        vector translationVector(this->centroid_[0], this->centroid_[1], this->centroid_[2]);
        this->points_ = this->translatePoints(rotatedPoints, translationVector);
    }

    boundBox aggregate::getBoundBox()
    {
        this->boundBox_ = boundBox(this->points_);
        point min = this->floorPoint(this->boundBox_.min());
        point max = this->ceilPoint(this->boundBox_.max());
        // this->roundedBoundBox_ = boundBox(min, max);
        return boundBox(min, max);
    }

}
