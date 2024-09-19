#include "List.H"
#include "scalar.H"
#include "Random.H"
#include <cmath>

namespace Bashyal
{
    class PSD
    {
    private:
        /* data */
    public:
        enum particleClass
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
        int nParticles_;
        Foam::List<Foam::List<float>> sizes_;
        Foam::List<float> defaultFraction_;
        Foam::List<float> fraction_;
        Foam::List<float> cumulativeFraction_;
        Foam::List<float> cumulativeParticleNumber_;

        PSD(/* args */);
        PSD(Foam::List<float> fraction, int nParticles);
        void initializeSizes();              // Function to initialize particle sizes
        void initializeDefaultProportions(); // Function to initialize particle sizes
        void calculateCSD();

        ~PSD() = default;
    };

    PSD::PSD(/* args */)
    {
    }

    inline PSD::PSD(Foam::List<float> fraction, int nParticles = 100)
    {
        this->fraction_ = fraction;
        this->nParticle_ = nParticles;
    }

    void PSD::initializeSizes()
    {
        // Initialize particle size ranges (in millimeters)

        sizes_[particleClass::Clay] = {0.001, 0.002};
        sizes_[particleClass::Silt] = {0.002, 0.075};
        sizes_[particleClass::Fine_Sand] = {0.075, 0.425};
        sizes_[particleClass::Medium_Sand] = {0.425, 2};
        sizes_[particleClass::Coarse_Sand] = {2 - 4.75};
        sizes_[particleClass::Fine_Gravel] = {4.75 - 9.5};
        sizes_[particleClass::Medium_Gravel] = {9.5, 19};
        sizes_[particleClass::Coarse_Gravel] = {19, 37.5};
        sizes_[particleClass::Large_Gravel] = {37.5, 75};
        sizes_[particleClass::Cobbels] = {75, 150};
    }

    void PSD::initializeDefaultProportions()
    {
        defaultFraction_[particleClass::Clay] = 5;
        defaultFraction_[particleClass::Silt] = 10;
        defaultFraction_[particleClass::Fine_Sand] = 10;
        defaultFraction_[particleClass::Medium_Sand] = 10;
        defaultFraction_[particleClass::Coarse_Sand] = 15;
        defaultFraction_[particleClass::Fine_Gravel] = 15;
        defaultFraction_[particleClass::Medium_Gravel] = 10;
        defaultFraction_[particleClass::Coarse_Gravel] = 10;
        defaultFraction_[particleClass::Large_Gravel] = 10;
        defaultFraction_[particleClass::Cobbels] = 5;
    }

    void PSD::calculateCSD()
    {
        float sum = 0;
        for (int i = 0; i < 10; i++)
        {
            cumulativeFraction_[i] = sum;
            sum = sum + fraction_[i];
            cumulativeParticleNumber_[i]=std::round(cumulativeFraction_[i]*nParticles_)
        }
    }

    PSD::~PSD()
    {
    }

}