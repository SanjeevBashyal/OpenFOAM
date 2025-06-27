#ifndef cubeAggregates_H
#define cubeAggregates_H

#include "List.H"
#include "Random.H"
#include "cubeAggregate.H"
#include "PSD.H"

namespace Bashyal
{
    class cubeAggregates
    {
    private:
        /* data */
    public:
        PSD PSD_;
        int nParticles_;
        Foam::PtrList<Bashyal::cubeAggregate> sediments_;

        cubeAggregates(PSD psd, int nParticles);
        void generateAggregates();
        ~cubeAggregates() = default;
    };


}

#endif

// ************************************************************************* //