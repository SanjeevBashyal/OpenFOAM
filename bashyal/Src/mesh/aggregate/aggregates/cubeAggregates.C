#include "cubeAggregates.H"

using namespace Foam;
namespace Bashyal
{

    cubeAggregates::cubeAggregates(PSD psd, int nParticles)
        : PSD_(psd),
          nParticles_(nParticles)
    {
    }

    void cubeAggregates::generateAggregates()
    {
        sediments_.setSize(nParticles_);
        labelList numberCSD = PSD_.calculateNumbers(nParticles_);

        Foam::Random randomGen(12);

        // Assuming PSD_.size() gives the number of size classes
        int sedimentCount = 0;
        for (Foam::label sizeClass = 0; sizeClass < PSD_.sizes_.size(); ++sizeClass)
        {
            // Define minimum and maximum sizes based on the PSD class
            // float minSize = PSD_.sizes_[sizeClass][0]; // Assuming PSD_ holds minimum sizes
            // float maxSize = PSD_.sizes_[sizeClass][1];

            for (int j = 0; j < numberCSD[sizeClass]; ++j)
            {
                Foam::scalar r = randomGen.sample01<Foam::scalar>();
                Foam::scalar r2 = randomGen.sample01<Foam::scalar>();

                float minSize = 0.1 + r * 0.15;
                float maxSize = 0.1 + 2 * r * 0.15;
                float orientation = 2 * 3.141592 * r2;
                // Create a random cubeAggregate within the size range
                Bashyal::cubeAggregate *newAggregate = new Bashyal::cubeAggregate(minSize, maxSize, sedimentCount);

                newAggregate->translate(vector(0.2 * (sedimentCount + 1), 0.3, 0.3));
                newAggregate->rotate(r2 / 3.141592 * 180 * sedimentCount, 0, 0);

                // Optionally, adjust location, orientation, etc., if needed
                newAggregate->locate();

                newAggregate->createFaces();

                // Add the created aggregate to the sediments_ list
                sediments_.set(sedimentCount, newAggregate);
                sedimentCount++;
            }
        }
    }

}
