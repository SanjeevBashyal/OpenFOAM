#include "introHeader.H"

// OpenFOAM Default include
#include "fvCFD_rest.H"
// bashyal include
#include "sediment_rest.H"

using namespace Bashyal;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote(
        "Coupled solver for sediment mix turbulent flows.");

#include "sedimentObjectsCreate.H"
#include "openfoamObjectsCreate.H"

    Info << "Starting time loop" << endl;

    turbulence->validate();

    if (!LTS)
    {
#include "CourantNo.H"
#include "setInitialDeltaT.H"
    }

    Info << "\nStarting time loop\n"
         << endl;

    while (runTime.run())
    {
#include "readDyMControls.H"

        if (LTS)
        {
#include "setRDeltaT.H"
        }
        else
        {
#include "CourantNo.H"
#include "setDeltaT.H"
        }

        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

#include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
#include "meshCourantNo.H"
                    }
                }
            }

#include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
#include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
