#ifndef TIME_REGISTRY_H
#define TIME_REGISTRY_H

#include "geomObject.H"
#include <vector>
#include <string>

namespace Bashyal {

class timeRegistry {
    std::vector<std::vector<geomObject>> objectTimeStack_; // [object][time]
    std::vector<std::string> objectNames_;
    scalar currentTime_;
    scalar timeStep_;
    scalar endTime_;
    std::string outputDir_;
    std::string pvdFileName_;
public:
    timeRegistry(scalar timeStep, scalar endTime, const std::string& outputDir = "VTK_output");

    void addObject(const geomObject& obj, const std::string& name);
    void advanceTime();
    void setCurrentTime(scalar t);
    scalar currentTime() const;
    scalar timeStep() const;
    scalar endTime() const;
    void writeTimeSeries() const;
};

}

#endif 