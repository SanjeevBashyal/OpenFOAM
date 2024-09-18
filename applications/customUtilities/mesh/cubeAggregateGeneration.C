#include "Time.H"
#include "IOdictionary.H"
#include "IOPtrList.H"

#include "argList.H"
#include "OSspecific.H"
#include "OFstream.H"

#include "wordPair.H"
#include "slidingInterface.H"
#include "scalar.H"

#include "cubeAggregate.H"


int main(int argc, char *argv[])
{
    Foam::argList::addNote
    (
        "Sediment Generator in given Particle Size Distribution.\n"
        "Work is ongoing....."
    );

    Foam::argList::noParallel();
    Foam::argList::noFunctionObjects();

    Foam::argList::addBoolOption
    (
        "show",
        "Display PSD"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    
    
    #include "getRegionOption.H"
    
        if (!Foam::polyMesh::regionName(regionName).empty())
        {
            Foam::Info<< Foam::nl << "Generating mesh for region " << regionName << Foam::nl;
        }

        // Instance for resulting mesh
        bool useTime = false;
        Foam::word meshInstance(runTime.constant());

        if
        (
            args.readIfPresent("time", meshInstance)
        && runTime.constant() != meshInstance
        )
        {
            // Verify that the value is actually good
            Foam::scalar timeValue;

            useTime = Foam::readScalar(meshInstance, timeValue);
            if (!useTime)
            {
                FatalErrorInFunction
                    << "Bad input value: " << meshInstance
                    << "Should be a scalar or 'constant'"
                    << Foam::nl << Foam::endl
                    << exit(Foam::FatalError);
            }
        }
    


    const Foam::word dictName("aggregateDict");

    Foam::autoPtr<Foam::IOdictionary> aggregateDictPtr;

    {
        Foam::fileName dictPath;
        const Foam::word& regionDir = Foam::polyMesh::regionName(regionName);

        if (args.readIfPresent("dict", dictPath))
        {
            // Dictionary specified on the command-line ...

            if (isDir(dictPath))
            {
                dictPath /= dictName;
            }
        }
        else
        {
            // Assume dictionary is to be found in the system directory

            dictPath = runTime.system()/regionDir/dictName;
        }

        Foam::IOobject aggregateDictIO
        (
            dictPath,
            runTime,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE,
            Foam::IOobject::NO_REGISTER
        );

        if (!aggregateDictIO.typeHeaderOk<Foam::IOdictionary>(true))
        {
            FatalErrorInFunction
                << aggregateDictIO.objectPath() << Foam::nl
                << exit(Foam::FatalError);
        }

        Foam::Info<< "Creating aggregates from "
            << aggregateDictIO.objectRelPath() << Foam::endl;

        aggregateDictPtr = Foam::autoPtr<Foam::IOdictionary>::New(aggregateDictIO);
    }

    const Foam::IOdictionary& aggregateDict = *aggregateDictPtr;
    


    Bashyal::cubeAggregate a(5.0,10.0);
    
}