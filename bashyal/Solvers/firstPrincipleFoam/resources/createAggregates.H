

const Foam::word dictName("aggregateDict");

Foam::autoPtr<Foam::IOdictionary> aggregateDictPtr;

Foam::word regionDir = "";
Foam::fileName dictPath = runTime.system() / regionDir / dictName;

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

const Foam::IOdictionary &aggregateDict = *aggregateDictPtr;

// PSD psd(aggregateDict);

// Foam::List<float> percentagePSD(10);
// percentagePSD[0] = 0;
// percentagePSD[1] = 0;
// percentagePSD[2] = 0;
// percentagePSD[3] = 0;
// percentagePSD[4] = 0;
// percentagePSD[5] = 0;
// percentagePSD[6] = 0;
// percentagePSD[7] = 0;
// percentagePSD[8] = 0;
// percentagePSD[9] = 100;

// PSD psd(percentagePSD);
// int nParticles = 5;
// cubeAggregates aggregates(psd, nParticles);
// aggregates.generateAggregates();

Bashyal::cubeAggregate a(0.16, 0.16, 1000); // size_1, size_2, identifier
a.translate(vector(0.15, 0.15, 0.15));
a.rotate(30, 30, 30);
a.locate();
a.createFaces();
a.getBoundBox();
a.createTriSurface();
