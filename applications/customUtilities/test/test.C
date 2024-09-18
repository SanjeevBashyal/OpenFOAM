#include "List.H"
#include "OFstream.H"

int main(int argc, char *argv[])
{
    Foam::List<int> a={1,2,3};
    Foam::Info<< a[1]<<Foam::endl;
}