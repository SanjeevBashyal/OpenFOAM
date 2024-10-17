#include "cubeAggregates.H"

namespace Bashyal
{

    cubeAggregates::cubeAggregates(Foam::List<float> PSD, int nParticles)
    {
        this->PSD_ = PSD;
        this->nParticles_ = nParticles;
    }

    // Foam::List<Bashyal::cubeAggregate> cubeAggregates::generateAggregates()
    // {
    //     for (int i = 0; i < this->nParticles_; i++)
    //     {
    //     }
    // }

    
}