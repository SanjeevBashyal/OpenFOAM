#include "particle.H"
#include "timeRegistry.H"
#include "OFstream.H"
#include "quaternion.H" 
#include "point.H"
#include "vector.H"
#include "faceList.H"
#include "pointField.H"
#include <cmath>

using namespace Bashyal;
using namespace Foam;

int main(int argc, char *argv[])
{
    // Simulation parameters
    scalar endTime = 1.0;
    scalar dt = 0.1;
    vector linearVelocity(0.5, 0.2, 0.0); // m/s
    vector angularVelocity(0, 0, 1.57);   // rad/s (90 deg/s around Z)
    point initialPosition(0, 0, 0);
    quaternion initialOrientation = quaternion::I;

    // Create a time registry
    timeRegistry registry(dt, endTime, "VTK_particle");

    // Create a particle (assuming a cube for geometry)
    pointField vertices(8);
    vertices[0] = point(-0.5, -0.5, -0.5);
    vertices[1] = point(0.5, -0.5, -0.5);
    vertices[2] = point(0.5, 0.5, -0.5);
    vertices[3] = point(-0.5, 0.5, -0.5);
    vertices[4] = point(-0.5, -0.5, 0.5);
    vertices[5] = point(0.5, -0.5, 0.5);
    vertices[6] = point(0.5, 0.5, 0.5);
    vertices[7] = point(-0.5, 0.5, 0.5);

    faceList faces(6);
    faces[0] = face({0, 3, 2, 1});
    faces[1] = face({4, 5, 6, 7});
    faces[2] = face({0, 1, 5, 4});
    faces[3] = face({2, 3, 7, 6});
    faces[4] = face({1, 2, 6, 5});
    faces[5] = face({3, 0, 4, 7});

    // Create a boundary object first, then the particle
    boundary cubeBoundary(vertices, faces);
    particle p(cubeBoundary, 1.0, tensor::I, initialPosition);

    // Simulation loop
    point currentPosition = initialPosition;
    quaternion currentOrientation = initialOrientation;
    for (scalar t = 0; t <= endTime + SMALL; t += dt)
    {
        // Transform particle geometry
        pointField worldVertices(vertices.size());
        forAll(vertices, i)
        {
            worldVertices[i] = currentPosition + currentOrientation.transform(vertices[i]);
        }
        // Create a geomObject for the current state
        geomObject g(worldVertices, faces);
        // Add to registry
        registry.addObject(g, "particle");
        // Advance time
        currentPosition += linearVelocity * dt;
        vector dTheta = angularVelocity * dt;
        quaternion deltaRotation(dTheta);
        currentOrientation = deltaRotation * currentOrientation;
        currentOrientation.normalise();
        registry.advanceTime();
    }
    // Write all time series files
    registry.writeTimeSeries();
    Info << "Simulation finished. Output written to VTK_particle/timeSeries.pvd" << endl;
    return 0;
} 