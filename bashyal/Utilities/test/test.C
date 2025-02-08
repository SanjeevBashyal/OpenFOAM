#include "quickInclude.H"
#include "roundAggregate.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    Info << "Here";
    Bashyal::roundAggregate a(0.1,10);
    a.translate(vector(0.2, 0.2, 0.2));
    a.rotate(30, 0, 0);
    a.locate();
    a.createFaces();
    a.createTriSurface();
    a.getBoundBox();
    Info << "Here";
}

