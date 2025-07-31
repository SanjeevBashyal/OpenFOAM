/*---------------------------------------------------------------------------*\
| =========                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
|  \\    /   O peration     | Version:  2312
|   \\  /    A nd           | Website:  https://openfoam.org
|    \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "vectorField.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading particle data...\n" << endl;

    // Read particle positions
    vectorField positions;
    {
        IOobject io
        (
            "positions",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        positions = vectorField(io);
    }

    // Read particle velocities
    vectorField velocities;
    {
        IOobject io
        (
            "velocities",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        velocities = vectorField(io);
    }

    // Read particle diameters
    scalarField diameters;
    {
        IOobject io
        (
            "diameters",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        diameters = scalarField(io);
    }

    Info<< "Number of particles: " << positions.size() << endl;
    Info<< "Particle positions: " << positions << endl;
    Info<< "Particle velocities: " << velocities << endl;
    Info<< "Particle diameters: " << diameters << endl;

    // Create a cloud for visualization
    fileName cloudDir = runTime.timePath()/"lagrangian"/"particleCloud";
    mkDir(cloudDir);

    // Write particle positions in cloud format
    {
        OFstream os(cloudDir/"positions");
        os << positions.size() << nl;
        forAll(positions, i)
        {
            os << positions[i].x() << ' ' 
               << positions[i].y() << ' ' 
               << positions[i].z() << nl;
        }
    }

    // Write particle properties
    {
        OFstream os(cloudDir/"properties");
        os << "positions" << nl;
        os << "velocities" << nl;
        os << "diameters" << nl;
    }

    Info<< "Particle data written to: " << cloudDir << endl;
    Info<< "You can now visualize this in ParaView by:" << endl;
    Info<< "1. Open ParaView" << endl;
    Info<< "2. Open the case file: " << runTime.caseName() << endl;
    Info<< "3. Enable lagrangian clouds in the pipeline" << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* // 