#include "quickInclude.H"
#include "backgroundMesh.H"

int main(int argc, char *argv[])
{
    Foam::argList::addNote(
        "backgroundMesh Generator to generate foundation mesh where fluids and sediments will traverse.\n"
        "Work is ongoing.....");

#include "setRootCase.H"
#include "createTime.H"

    Bashyal::backgroundMesh a(&runTime, 1.0);

    Foam::Info << "Run Successfully" << Foam::endl;

    return 0;
}