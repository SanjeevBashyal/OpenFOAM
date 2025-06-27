#include "wellRoundedAggregate.H"
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
    // Constructor
    wellRoundedAggregate::wellRoundedAggregate(float s1, float s2, float roundnessFactor, int resolution, Foam::label identifier)
        : s1_(s1), s2_(s2), roundnessFactor_(roundnessFactor), resolution_(resolution), identifier_(identifier)
    {
        // Validate inputs
        if (s1 <= 0 || s2 <= 0 || s1 > s2 || roundnessFactor < 0 || roundnessFactor > 1 || resolution < 3)
        {
            FatalErrorInFunction << "Invalid parameters: s1=" << s1 << ", s2=" << s2 
                                 << ", roundnessFactor=" << roundnessFactor << ", resolution=" << resolution 
                                 << nl << exit(FatalError);
        }

        // Random number generator
        Foam::Random randomGen(12); // Fixed seed for reproducibility

        // Random size between s1 and s2
        Foam::scalar r = randomGen.sample01<Foam::scalar>();
        Foam::scalar s = s1 + (s2 - s1) * r;

        // Scaling factors based on roundness factor
        Foam::scalar a = 1.0 - roundnessFactor; // 0 when rf=1 (sphere), 1 when rf=0 (max variation)
        Foam::scalar sx = 1.0 + a * (2.0 * randomGen.sample01<Foam::scalar>() - 1.0); // Range [1-a, 1+a]
        Foam::scalar sy = 1.0 + a * (2.0 * randomGen.sample01<Foam::scalar>() - 1.0);
        Foam::scalar sz = 1.0 + a * (2.0 * randomGen.sample01<Foam::scalar>() - 1.0);

        // Generate and store points
        this->localPoints_ = generatePoints(s, sx, sy, sz, resolution);

        // Generate and store faces
        this->createFaces();

        // Initialize centroid
        this->centroid_ = Foam::point(0, 0, 0);
    }

    // Generate points for the well-rounded aggregate
    Foam::pointField wellRoundedAggregate::generatePoints(Foam::scalar s, Foam::scalar sx, Foam::scalar sy, Foam::scalar sz, int resolution)
    {
        // Generate points on a unit sphere
        Foam::pointField points = createSpherePoints(1.0, resolution);

        // Apply scaling along each axis
        for (auto& p : points)
        {
            p.x() *= sx;
            p.y() *= sy;
            p.z() *= sz;
        }

        // Compute bounding box and scale to size s
        Foam::boundBox bb(points);
        Foam::scalar maxExtent = max(max(bb.max().x() - bb.min().x(), bb.max().y() - bb.min().y()), bb.max().z() - bb.min().z());
        Foam::scalar scaleFactor = s / maxExtent;

        for (auto& p : points)
        {
            p *= scaleFactor;
        }

        return points;
    }

    // Create points on a sphere (adapted from roundAggregate)
    Foam::pointField wellRoundedAggregate::createSpherePoints(Foam::scalar radius, int resolution)
    {
        int pointsPerRing = resolution;
        int totalPoints = (resolution - 1) * pointsPerRing + 2; // Poles + intermediate rings
        Foam::pointField points(totalPoints);

        // North pole
        points[0] = Foam::point(0, 0, radius);

        // Intermediate rings
        for (int i = 1; i < resolution; ++i)
        {
            Foam::scalar theta = Foam::constant::mathematical::pi * i / resolution;
            Foam::scalar sinTheta = Foam::sin(theta);
            Foam::scalar cosTheta = Foam::cos(theta);

            for (int j = 0; j < resolution; ++j)
            {
                Foam::scalar phi = 2 * Foam::constant::mathematical::pi * j / resolution;
                Foam::scalar x = radius * sinTheta * Foam::cos(phi);
                Foam::scalar y = radius * sinTheta * Foam::sin(phi);
                Foam::scalar z = radius * cosTheta;
                int index = (i - 1) * pointsPerRing + j + 1;
                points[index] = Foam::point(x, y, z);
            }
        }

        // South pole
        points[totalPoints - 1] = Foam::point(0, 0, -radius);

        return points;
    }

    // Create faces
    void wellRoundedAggregate::createFaces()
    {
        this->faces_ = createQuadFaces();
    }

    // Generate quadrilateral faces
    Foam::List<Foam::face> wellRoundedAggregate::createQuadFaces()
    {
        Foam::List<Foam::face> faces;
        int pointsPerRing = resolution_;

        // Top triangular faces (north pole)
        for (int j = 0; j < resolution_; ++j)
        {
            int p1 = 0;                         // North pole
            int p2 = j + 1;                     // Current point on first ring
            int p3 = (j + 1) % resolution_ + 1; // Next point on first ring
            faces.append(Foam::face({p1, p2, p3}));
        }

        // Intermediate quadrilateral faces
        int offset = 1;
        for (int i = 1; i < resolution_ - 1; ++i)
        {
            for (int j = 0; j < resolution_; ++j)
            {
                int p1 = offset + j;
                int p2 = offset + (j + 1) % resolution_;
                int p3 = offset + pointsPerRing + (j + 1) % resolution_;
                int p4 = offset + pointsPerRing + j;
                faces.append(Foam::face({p1, p2, p3, p4}));
            }
            offset += pointsPerRing;
        }

        // Bottom triangular faces (south pole)
        int bottomOffset = (resolution_ - 2) * pointsPerRing + 1;
        int southPole = bottomOffset + resolution_;
        for (int j = 0; j < resolution_; ++j)
        {
            int p1 = southPole;                    // South pole
            int p2 = bottomOffset + j;             // Current point on last ring
            int p3 = bottomOffset + (j + 1) % resolution_; // Next point
            faces.append(Foam::face({p1, p2, p3}));
        }

        return faces;
    }

}
