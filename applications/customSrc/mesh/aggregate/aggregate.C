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
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
typedef Delaunay::Point Point;

Bashyal::aggregate::aggregate(float r1, float r2, int n)
{

    // Create an instance of the Random class with automatic seeding
    Foam::Random randomGen(12);

    Foam::pointField Nodes(n);
    this->r1 = r1;
    Foam::scalar r = 0;
    Foam::scalar phi = 0;
    Foam::scalar theta = 0;
    float addTheta = (2 * 3.141592) / n;
    float addPhi = (3.141592) / n;

    std::vector<Point> points;

    for (int i = 0; i < n; i++)
    {
        Foam::scalar randomPhi = randomGen.sample01<Foam::scalar>();
        Foam::scalar randomTheta = randomGen.sample01<Foam::scalar>();
        Foam::scalar randomR = randomGen.sample01<Foam::scalar>();

        r = r1 + (r2 - r1) * randomR;

        theta = theta + randomTheta * addTheta;
        phi = phi + randomPhi * addPhi;

        // Generate a random scalar between 0 and 1 using the scalar01() method

        Foam::point xyz = this->sphericalToCartesian(r, theta, phi);

        points.push_back(Point(xyz[0], xyz[1], xyz[2]));

        Nodes.append(xyz);
        Foam::faceList faces = this->createFacesFromPoints(Nodes);
        // Foam::triSurface surface(faces, Nodes);
        // Foam::fileName outputFile = argList::globalArgs()[1];  // Second argument is the output file name
        // surface.write(outputFile);
    }

    Delaunay dt;
    dt.insert(points.begin(), points.end());
    int c = 0;
}

Foam::point Bashyal::aggregate::sphericalToCartesian(Foam::scalar r, Foam::scalar theta, Foam::scalar phi)
{
    // Using OpenFOAM's point structure for Cartesian coordinates
    Foam::scalar x = r * sin(theta) * cos(phi);
    Foam::scalar y = r * sin(theta) * sin(phi);
    Foam::scalar z = r * cos(theta);

    // Return the result as a point
    return Foam::point(x, y, z);
}

Foam::scalar Bashyal::aggregate::distance(const Foam::point &p1, const Foam::point &p2)
{
    return mag(p2 - p1);
}

// Function to find the nearest points and form triangles
Foam::faceList Bashyal::aggregate::createFacesFromPoints(const Foam::pointField &points)
{
    Foam::faceList faces;

    // We will create a very simple triangulation by connecting each point with its two closest neighbors
    for (Foam::label i = 0; i < points.size(); i++)
    {
        std::vector<std::pair<Foam::scalar, Foam::label>> distances;

        // Calculate the distance of point 'i' from every other point
        for (Foam::label j = 0; j < points.size(); j++)
        {
            if (i != j)
            {
                distances.push_back(std::make_pair(distance(points[i], points[j]), j));
            }
        }

        // Sort the distances to find the nearest two neighbors
        std::sort(distances.begin(), distances.end());

        // Form a face from the point and its two nearest neighbors
        Foam::face newFace(3); // Triangular face
        newFace[0] = i;
        newFace[1] = distances[0].second; // First nearest neighbor
        newFace[2] = distances[1].second; // Second nearest neighbor

        // Add the face to the list
        faces.append(newFace);
    }

    return faces;
}

namespace Bashyal
{
    int hello(int a)
    {
        return a + 1;
    }
}
