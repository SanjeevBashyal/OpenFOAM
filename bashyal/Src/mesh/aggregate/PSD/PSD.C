#include "PSD.H"

using namespace Foam;
namespace Bashyal
{
    PSD::PSD(/* args */)
    {
    }

    PSD::PSD(IOdictionary aggregateDict)
    {
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
    }

    PSD::PSD(Foam::List<float> percentages)
    {
        this->initialize();
        this->percentagePSD_ = percentages;
    }

    void PSD::initializeSizes()
    {
        // Initialize particle size ranges (in millimeters)
        sizes_.setSize(nCategories_);
        sizes_[particleClass_::Clay] = {0.001, 0.002};
        sizes_[particleClass_::Silt] = {0.002, 0.075};
        sizes_[particleClass_::Fine_Sand] = {0.075, 0.425};
        sizes_[particleClass_::Medium_Sand] = {0.425, 2};
        sizes_[particleClass_::Coarse_Sand] = {2, 4.75};
        sizes_[particleClass_::Fine_Gravel] = {4.75, 9.5};
        sizes_[particleClass_::Medium_Gravel] = {9.5, 19};
        sizes_[particleClass_::Coarse_Gravel] = {19, 37.5};
        sizes_[particleClass_::Large_Gravel] = {37.5, 75};
        sizes_[particleClass_::Cobbels] = {75, 150};
    }

    void PSD::initializeDefaultProportions()
    {
        defaultFraction_.setSize(nCategories_);
        defaultFraction_[particleClass_::Clay] = 5;
        defaultFraction_[particleClass_::Silt] = 10;
        defaultFraction_[particleClass_::Fine_Sand] = 10;
        defaultFraction_[particleClass_::Medium_Sand] = 10;
        defaultFraction_[particleClass_::Coarse_Sand] = 15;
        defaultFraction_[particleClass_::Fine_Gravel] = 15;
        defaultFraction_[particleClass_::Medium_Gravel] = 10;
        defaultFraction_[particleClass_::Coarse_Gravel] = 10;
        defaultFraction_[particleClass_::Large_Gravel] = 10;
        defaultFraction_[particleClass_::Cobbels] = 5;
    }

    List<label> PSD::calculateNumbers(int nParticles)
    {
        this->calculateCSD();
        labelList numberPSD;
        numberPSD.setSize(percentageCSD_.size());
        for (int i = 0; i < 10; i++)
        {
            numberPSD[i] = std::round(percentageCSD_[i]/100 * nParticles);
        }
        return numberPSD;
    }

    void PSD::calculateCSD()
    {
        float sum = 0;
        percentageCSD_.setSize(nCategories_);
        for (int i = 0; i < nCategories_; i++)
        {
            sum = sum + percentagePSD_[i];
            percentageCSD_[i] = sum;
        }
    }

    void PSD::initialize()
    {
        this->initializeSizes();
        percentagePSD_.setSize(nCategories_);
        percentageCSD_.setSize(nCategories_);
    }

    PSD::~PSD()
    {
    }
}