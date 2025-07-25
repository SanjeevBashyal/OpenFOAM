#include "particle.H"
#include "timeRegistry.H"
#include <OpenFOAM/OFstream.H>
#include <OpenFOAM/quaternion.H>
#include <OpenFOAM/point.H>
#include <OpenFOAM/vector.H>
#include <OpenFOAM/faceList.H>
#include <OpenFOAM/pointField.H>
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
    quaternion initialOrientation = quaternion::one;

    // Create a time registry
    timeRegistry registry(dt, endTime, "VTK_particle");

    // Create a particle (assuming a cube for geometry)
    pointField vertices = {
        {-0.5, -0.5, -0.5}, {0.5, -0.5, -0.5}, {0.5, 0.5, -0.5}, {-0.5, 0.5, -0.5},
        {-0.5, -0.5, 0.5}, {0.5, -0.5, 0.5}, {0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}
    };
    faceList faces = {
        {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
        {2, 3, 7, 6}, {1, 2, 6, 5}, {3, 0, 4, 7}
    };
    particle p(initialPosition, initialOrientation, vertices, faces);

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
        currentOrientation.normalize();
        registry.advanceTime();
    }
    // Write all time series files
    registry.writeTimeSeries();
    Info << "Simulation finished. Output written to VTK_particle/timeSeries.pvd" << endl;
    return 0;
} 