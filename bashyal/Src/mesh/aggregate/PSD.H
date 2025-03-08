#ifndef PSD_H
#define PSD_H

#include "quickInclude.H"


namespace Bashyal
{
    class PSD
    {
    private:
        /* data */
    public:
        int nCategories_ = 10;
        enum particleClass_
        {
            Clay,
            Silt,
            Fine_Sand,
            Medium_Sand,
            Coarse_Sand,
            Fine_Gravel,
            Medium_Gravel,
            Coarse_Gravel,
            Large_Gravel,
            Cobbels,
        };

        Foam::List<Foam::List<float>> sizes_;
        Foam::List<float> defaultFraction_;
        Foam::List<float> percentagePSD_;
        Foam::List<float> percentageCSD_;

        PSD(/* args */);
        PSD(Foam::IOdictionary aggregateDict);
        PSD(Foam::List<float> percentages);

        void initialize();
        void initializeSizes();              // Function to initialize particle sizes
        void initializeDefaultProportions(); // Function to initialize particle sizes
        void calculateCSD();
        Foam::List<Foam::label> calculateNumbers(int nParticles);

        ~PSD();
    };

    // inline PSD::PSD(Foam::List<float> fraction, int nParticles = 100)
    // {
    //     this->percentagePSD_ = fraction;
    //     this->nParticle_ = nParticles;
    // }

}

#endif

// ************************************************************************* //