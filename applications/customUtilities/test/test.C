#include "quickInclude.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    Vector<int> v(0, 0, 1);
    Vector<int> dim(4, 2, 3);
    Foam::Vector<int> c = (v + 1) * dim;
    Info << "Here";
}