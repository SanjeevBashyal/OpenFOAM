

const Foam::word dictName("aggregateDict");

Foam::autoPtr<Foam::IOdictionary> aggregateDictPtr;

Foam::word regionName
(
    args.getOrDefault<word>("region", Foam::polyMesh::defaultRegion)
);

    {
        Foam::fileName dictPath;
        const Foam::word &regionDir = Foam::polyMesh::regionName(regionName);

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

            dictPath = runTime.system() / regionDir / dictName;
        }

        Foam::IOobject aggregateDictIO(
            dictPath,
            runTime,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE,
            Foam::IOobject::NO_REGISTER);

        if (!aggregateDictIO.typeHeaderOk<Foam::IOdictionary>(true))
        {
            FatalErrorInFunction
                << aggregateDictIO.objectPath() << Foam::nl
                << exit(Foam::FatalError);
        }

        Foam::Info << "Creating aggregates from "
                   << aggregateDictIO.objectRelPath() << Foam::endl;

        aggregateDictPtr = Foam::autoPtr<Foam::IOdictionary>::New(aggregateDictIO);
    }

    const Foam::IOdictionary &aggregateDict = *aggregateDictPtr;

    if (aggregateDict.found("PSD"))
    {
        Foam::dictionary psdDict = aggregateDict.subDict("PSD");
        if (psdDict.found("Clay(<0.002)"))
        {

            // Foam::dictionary clayDict = psdDict.subDict("Clay(<0.002)");
            Foam::scalar clayFraction = psdDict.getOrDefault("Clay(<0.002)", 2);
            Foam::Info << "Here" << Foam::endl;
        }

        // Print the keys and values in the "PSD" sub-dictionary
        Foam::Info << "Contents of sub-dictionary 'PSD':" << Foam::endl;

        forAllConstIter(Foam::dictionary, psdDict, iter)
        {
            // Print each key and value
            // Foam::Info << "Key: " << iter.key()
            //            << ", Value: " << iter() << Foam::endl;
        }
    }
    else
    {
        Foam::Warning << "Sub-dictionary 'PSD' not found!" << Foam::endl;
    }

    Foam::List<float> PSD(10);
    PSD[0] = 5.0;
    PSD[1] = 10.0;
    PSD[2] = 10.0;
    PSD[3] = 15.0;
    PSD[4] = 15.0;
    PSD[5] = 10.0;
    PSD[6] = 10.0;
    PSD[7] = 10.0;
    PSD[8] = 10.0;
    PSD[9] = 5.0;

    int nParticles = 100;

    // Bashyal::cubeAggregates sediments(PSD, nParticles);

    Bashyal::cubeAggregate a(0.1, 0.1);
    a.translate(vector(0.1, 0.1, 0.1));
    a.rotate(30, 0, 0);
    a.locate();
    a.createSurface();
    a.getBoundBox();

