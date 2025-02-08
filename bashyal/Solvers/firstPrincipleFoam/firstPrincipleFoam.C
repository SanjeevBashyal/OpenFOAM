#include "introHeader.H"

// OpenFOAM Default include
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

// Custom include
#include "quickInclude.H"
#include "backgroundMesh.H"
#include "backgroundBlock.C"
#include "cubeAggregate.H"
#include "roundAggregate.H"
#include "PSD.H"
#include "cubeAggregates.H"

using namespace Bashyal;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote(
        "Coupled solver for sediment mix turbulent flows.");

    // #include "postProcess.H"
    // #include "addCheckCaseOptions.H"
    // #include "setRootCaseLists.H"
    // #include "createTime.H"
    // #include "createMesh.H"
    // #include "createControl.H"
    // #include "createFields.H"
    // #include "initContinuityErrs.H"

#include "addRegionOption.H"
#include "setRootCase.H" //defines argList args
#include "createTime.H"
#include "createBackgroundMesh.H"
#include "createAggregates.H"

#include "createDebug.H"

    // bMesh.intersectCubes(aggregates);

    bMesh.intersectCube(a);
    // bMesh.intersectRound(a);

    bMesh.developMesh();
    bMesh.writeBackgroundMesh(runDir / meshDir0);
    Foam::Info << "Run Successfully" << Foam::endl;

#include "createInitMesh.H"
#include "createControl.H"
#include "createFields.H"
#include "initContinuityErrs.H"

    Info << "\nStarting time loop\n"
         << endl;

    volVectorField constantVectorField(
        IOobject(
            "myVectorField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector("constantVectorField", dimVelocity, vector(1.0, 0.0, 0.0)));

    const dimensionSet dimOne(0, 2, 0, 0, 0, 0, 0);
    const dimensionSet dimPressure2(0, 2, -2, 0, 0, 0, 0);

    volScalarField constantScalarField(
        IOobject(
            "constantScalarField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("constantScalarField", dimPressure2, scalar(-5.0)));

    volScalarField one(
        IOobject(
            "one",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("one", dimOne, scalar(1.0)));

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        fvScalarMatrix pEqn(
            fvm::laplacian(one, p));

        // fvVectorMatrix pEqn
        // (
        //     fvm::ddt(p)
        // );

        // volVectorField gradP = fvc::grad(p);
        // Solve the equation
        solve(pEqn == constantScalarField);

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info << "End\n"
         << endl;

    return 0;

    //     turbulence->validate();

    //     // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //     Info << "\nStarting time loop\n"
    //          << endl;

    //     while (simple.loop())
    //     {
    //         Info << "Time = " << runTime.timeName() << nl << endl;

    //         // --- Pressure-velocity SIMPLE corrector
    //         {
    // // #include "UEqn.H"
    // // #include "pEqn.H"
    //         }

    //         laminarTransport.correct();
    //         turbulence->correct();

    //         runTime.write();

    //         runTime.printExecutionTime(Info);
    //     }

    //     Info << "End\n"
    //          << endl;

    // return 0;
}

// ************************************************************************* //
