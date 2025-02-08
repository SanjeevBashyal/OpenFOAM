#include "quickInclude.H"
#include "backgroundMesh.H"
#include "backgroundBlock.C"

using namespace Foam;
using namespace Bashyal;

int main(int argc, char *argv[])
{
    Foam::argList::addNote(
        "backgroundMesh Generator to generate foundation mesh where fluids and sediments will traverse.\n"
        "Work is ongoing.....");

#include "setRootCase.H" //defines argList args
#include "createTime.H"

    Bashyal::backgroundMesh a(&runTime, point(0, 0, 0), point(0.2, 0.1, 0.1), 0.1);
    // a.createBlockMesh();

    a.developMesh();
    Foam::fileName runDir = args.path();
    fileName meshDir0 = "constant/polyMesh";
    a.writeBackgroundMesh(runDir / meshDir0);
    Foam::Info << "Run Successfully" << Foam::endl;

    return 0;
}
