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

#include "setRootCase.H"    //defines argList args
#include "createTime.H"

    Bashyal::backgroundMesh a(&runTime, 1.0);
    boundBox boundBox0 = boundBox(point(0.1, 0.1, 0.1), point(0.2, 0.2, 0.2));
    backgroundBlock b = backgroundBlock(&a, 5, boundBox0);

    Foam::fileName runDir = args.path();
    fileName meshDir0="constant/polyMesh";
    b.write(runDir/meshDir0);
    Foam::Info << "Run Successfully" << Foam::endl;

    return 0;
}