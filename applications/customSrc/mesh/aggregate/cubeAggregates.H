#ifndef cubeAggregates_H
#define cubeAggregates_H

#include "List.H"
#include "Random.H"
#include "cubeAggregate.H"

namespace Bashyal
{
    class cubeAggregates
    {
    private:
        /* data */
    public:
        Foam::List<float> PSD_;
        int nParticles_;
        Foam::List<Bashyal::cubeAggregate> sediments_;

        cubeAggregates(Foam::List<float> PSD, int nParticles);
        Foam::List<Bashyal::cubeAggregate> generateAggregates();
        ~cubeAggregates() = default;
    };

    cubeAggregates::cubeAggregates(Foam::List<float> PSD, int nParticles)
    {
        this->PSD_ = PSD;
        this->nParticles_ = nParticles;
    }

    Foam::List<Bashyal::cubeAggregate> cubeAggregates::generateAggregates()
    {
        for (int i = 0; i < this->nParticles_; i++)
        {
        }
    }
}

#endif

// ************************************************************************* //