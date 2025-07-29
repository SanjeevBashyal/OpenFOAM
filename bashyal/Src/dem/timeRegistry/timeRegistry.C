#include "timeRegistry.H"
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace Bashyal {

timeRegistry::timeRegistry(Foam::scalar timeStep, Foam::scalar endTime, const std::string& outputDir)
    : currentTime_(0.0), timeStep_(timeStep), endTime_(endTime), outputDir_(outputDir), pvdFileName_(outputDir+"/timeSeries.pvd")
{
    mkdir(outputDir_.c_str(), 0777);
}

void timeRegistry::addObject(const geomObject& obj, const std::string& name) {
    if (objectTimeStack_.size() < objectNames_.size() + 1) {
        objectTimeStack_.emplace_back();
    }
    objectTimeStack_.back().push_back(obj);
    objectNames_.push_back(name);
}

void timeRegistry::advanceTime() {
    currentTime_ += timeStep_;
}

void timeRegistry::setCurrentTime(Foam::scalar t) {
    currentTime_ = t;
}

Foam::scalar timeRegistry::currentTime() const { return currentTime_; }
Foam::scalar timeRegistry::timeStep() const { return timeStep_; }
Foam::scalar timeRegistry::endTime() const { return endTime_; }

void timeRegistry::writeTimeSeries() const {
    // Write a .pvd file referencing all .vtp files for each object
    std::ofstream pvd(pvdFileName_);
    pvd << "<?xml version=\"1.0\"?>\n";
    pvd << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    pvd << "  <Collection>\n";
    for (size_t objIdx = 0; objIdx < objectTimeStack_.size(); ++objIdx) {
        const auto& timeSeries = objectTimeStack_[objIdx];
        const std::string& name = objectNames_[objIdx];
        for (size_t tIdx = 0; tIdx < timeSeries.size(); ++tIdx) {
            std::ostringstream vtpName;
            vtpName << outputDir_ << "/" << name << "_" << std::setw(5) << std::setfill('0') << tIdx << ".vtp";
            timeSeries[tIdx].writeVtp(vtpName.str());
            pvd << "    <DataSet timestep=\"" << tIdx * timeStep_ << "\" group=\"" << name << "\" part=\"0\" file=\"" << name << "_" << std::setw(5) << std::setfill('0') << tIdx << ".vtp\"/>\n";
        }
    }
    pvd << "  </Collection>\n";
    pvd << "</VTKFile>\n";
    pvd.close();
}

} 
