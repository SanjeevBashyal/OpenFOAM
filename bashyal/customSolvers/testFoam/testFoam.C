#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

#include "quickInclude.H"
#include "backgroundMesh.H"
#include "backgroundBlock.C"
#include "cubeAggregate.H"
#include "roundAggregate.H"
#include "PSD.H"
#include "cubeAggregates.H"

using namespace Bashyal;

int main(int argc, char *argv[])
{
    argList::addNote("Coupled solver for sediment mix turbulent flows.");

#include "addRegionOption.H"
#include "setRootCase.H" // Defines argList args
#include "createTime.H"
#include "createBackgroundMesh.H"
#include "createAggregates.H"

    // Time step list
    scalar timeStep = 0.0501;  // seconds
    scalar totalTime = 0.9520; // total simulation time in seconds
    scalar currentTime = 0.0;

    // Predefined positions for translation (time step positions)
    List<vector> positions{
        vector(0.0045, 0.3905, 0.0045),
        vector(0.0045, 0.3835, 0.0045),
        vector(0.0045, 0.3731, 0.0045),
        vector(0.0045, 0.3632, 0.0045),
        vector(0.0045, 0.3533, 0.0045),
        vector(0.0045, 0.3433, 0.0045),
        vector(0.0045, 0.3334, 0.0045),
        vector(0.0045, 0.3235, 0.0045),
        vector(0.0045, 0.3136, 0.0045),
        vector(0.0045, 0.3036, 0.0045),
        vector(0.0045, 0.2937, 0.0045),
        vector(0.0045, 0.2838, 0.0045),
        vector(0.0045, 0.2739, 0.0045),
        vector(0.0045, 0.2639, 0.0045),
        vector(0.0045, 0.2543, 0.0045),
        vector(0.0045, 0.2441, 0.0045),
        vector(0.0045, 0.2342, 0.0045),
        vector(0.0045, 0.2242, 0.0045),
        vector(0.0045, 0.2143, 0.0045),
        vector(0.0045, 0.2044, 0.0045)};

    // Ensure the size of positions matches the number of time steps
    if (positions.size() != scalar(totalTime / timeStep + 1))
    {
        Foam::Info << "Size Mismatch" << Foam::endl;
    }

    // Initialize rotation angle
    scalar rotationAngle = 0.0; // Degrees

    while (currentTime < totalTime)
    {
        Info << "Current Time = " << currentTime << nl << endl;

        // Find the corresponding position for the current time step
        label timeIndex = round(currentTime / timeStep);
        vector translationVector = positions[timeIndex];

        // Update cube properties
        a.position(translationVector); // Update position
        a.rotate(rotationAngle, 0, 0);  // Apply rotation
        a.locate();
        a.createFaces();

        // Update and develop mesh
        bMesh.resetBlocks();
        bMesh.intersectCube(a);
        bMesh.developMesh();

        // Write mesh to new folder for the current time step
        word timeFolderName = Foam::name(currentTime);                // Convert time to folder name
        fileName currentMeshDir = runDir / timeFolderName / polyMeshThis; // Mesh folder for this timestep
        Foam::mkDir(currentMeshDir);                                  // Create the directory if it doesn't exist
        bMesh.writeBackgroundMesh(currentMeshDir);                    // Write the mesh

        Foam::Info << "Mesh written to: " << currentMeshDir << Foam::endl;

        // Advance time
        currentTime += timeStep;
    }

    Foam::Info << "Simulation completed successfully!" << Foam::endl;
    return 0;
}
