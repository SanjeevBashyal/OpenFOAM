#ifndef TIME_REGISTRY_H
#define TIME_REGISTRY_H

#include "geomObject.H"
#include <vector>
#include <string>

namespace Bashyal {

class timeRegistry {
    std::vector<std::vector<geomObject>> objectTimeStack_; // [object][time]
    std::vector<std::string> objectNames_;
    Foam::scalar currentTime_;
    Foam::scalar timeStep_;
    Foam::scalar endTime_;
    std::string outputDir_;
    std::string pvdFileName_;
public:
    timeRegistry(Foam::scalar timeStep, Foam::scalar endTime, const std::string& outputDir = "VTK_output");

    void addObject(const geomObject& obj, const std::string& name);
    void advanceTime();
    void setCurrentTime(Foam::scalar t);
    Foam::scalar currentTime() const;
    Foam::scalar timeStep() const;
    Foam::scalar endTime() const;
    void writeTimeSeries() const;
};

}

#endif 
