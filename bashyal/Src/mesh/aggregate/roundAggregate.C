#include "roundAggregate.H"
#include "point.H"
#include "scalar.H"
#include "pointField.H"
#include "faceList.H"
#include "labelledTri.H"

#include "Random.H"
#include "transform.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"

using namespace Foam;

namespace Bashyal
{
    roundAggregate::roundAggregate()
    {
    }

    roundAggregate::roundAggregate(float radius, int resolution, Foam::label identifier)
        : radius_(radius),
          resolution_(resolution), identifier_(identifier)
    {
        this->localPoints_ = createSpherePoints(radius, resolution); // Generate sphere points with default resolution
        this->centroid_ = Foam::point(0, 0, 0);                      // Initialize centroid
    }

    void roundAggregate::createFaces()
    {
        this->createQuadFaces();
    }

    Foam::List<Foam::face> roundAggregate::createQuadFaces()
    {
        Foam::List<Foam::face> faces;
        int pointsPerRing = resolution_;

        // Create triangular faces for the top (first ring)
        for (int j = 0; j < resolution_; ++j)
        {
            int p1 = 0;                         // Top point (north pole)
            int p2 = j + 1;                     // Current point on the first ring
            int p3 = (j + 1) % resolution_ + 1; // Next point on the first ring (wrap around)
            faces.append(Foam::face({p1, p2, p3}));
        }

        // Create quadrilateral faces for intermediate rings
        int offset = 1;
        for (int i = 1; i < resolution_ - 1; ++i)
        {
            for (int j = 0; j < resolution_; ++j)
            {
                // Compute indices for four vertices of the quadrilateral
                int p1 = offset + j;
                int p2 = offset + (j + 1) % resolution_;
                int p3 = (offset + pointsPerRing) + (j + 1) % resolution_;
                int p4 = (offset + pointsPerRing) + j;

                // Create the face
                faces.append(Foam::face({p1, p2, p3, p4}));
            }
            offset = offset + pointsPerRing;
        }

        // Create triangular faces for the bottom (last ring)
        int bottomOffset = (resolution_ - 2) * pointsPerRing + 1;
        for (int j = 0; j < resolution_; ++j)
        {
            int p1 = bottomOffset + resolution_;           // Bottom point (south pole)
            int p2 = bottomOffset + j;                     // Current point on the last ring
            int p3 = bottomOffset + (j + 1) % resolution_; // Next point on the last ring (wrap around)
            faces.append(Foam::face({p1, p2, p3}));
        }

        this->faces_ = faces;
        return faces;
    }

    Foam::pointField roundAggregate::createSpherePoints(Foam::scalar radius, int resolution)
    {
        // Total number of points
        int pointsPerRing = resolution;                         // Points in each ring (including wrap-around)
        int totalPoints = (resolution - 1) * pointsPerRing + 2; // Intermediate + top/bottom poles

        Foam::pointField points(totalPoints);

        // Add the north pole
        points[0] = Foam::point(0, 0, radius);

        // Generate points for intermediate rings
        for (int i = 1; i < resolution; ++i) // Loop over latitude divisions (excluding poles)
        {
            // Polar angle (theta) from 0 to π
            Foam::scalar theta = Foam::constant::mathematical::pi * i / resolution;

            Foam::scalar sinTheta = Foam::sin(theta);
            Foam::scalar cosTheta = Foam::cos(theta);

            for (int j = 0; j < resolution; ++j) // Loop over longitude divisions
            {
                // Azimuthal angle (phi) from 0 to 2π
                Foam::scalar phi = 2 * Foam::constant::mathematical::pi * j / resolution;

                Foam::scalar x = radius * sinTheta * Foam::cos(phi);
                Foam::scalar y = radius * sinTheta * Foam::sin(phi);
                Foam::scalar z = radius * cosTheta;

                int index = (i - 1) * pointsPerRing + j + 1;
                points[index] = Foam::point(x, y, z);
            }
        }

        // Add the south pole
        points[totalPoints - 1] = Foam::point(0, 0, -radius);

        return points;
    }

    void roundAggregate::createTriSurface()
    {
        Foam::List<Foam::labelledTri> triangles;
        for (const auto &face : this->faces_)
        {
            if (face.size() == 4)
            {
                // Split each quadrilateral face into two triangles
                triangles.append(Foam::labelledTri(face[0], face[1], face[2]));
                triangles.append(Foam::labelledTri(face[0], face[2], face[3]));
            }
            else
            {
                triangles.append(Foam::labelledTri(face[0], face[1], face[2]));
            }
        }

        Foam::triSurface surface(triangles, this->localPoints_);

        Foam::fileName outputFile("/usr/lib/openfoam/openfoam2312/run/test/roundAggregate.stl"); // Output file
        surface.write(outputFile);

        this->surface_ = surface;
    }
}
