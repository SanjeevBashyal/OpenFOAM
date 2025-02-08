/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    simpleFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state solver for incompressible, turbulent flows.

    \heading Solver details
    The solver uses the SIMPLE algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

\*---------------------------------------------------------------------------*/

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
