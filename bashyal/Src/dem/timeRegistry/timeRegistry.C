#include "timeRegistry.H"
#include "constants.H"
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace Bashyal {

timeRegistry::timeRegistry(Foam::scalar timeStep, Foam::scalar endTime, const std::string& outputDir)
    : currentTime_(0.0), timeStep_(timeStep), endTime_(endTime), outputDir_(outputDir), pvdFileName_(outputDir+"/timeSeries.pvd"),
      writeInterval_(1.0), lastWriteTime_(-1.0), writeAtStart_(true), writeAtEnd_(true)
{
    mkdir(outputDir_.c_str(), 0777);
}

void timeRegistry::addParticle(particle& p, const std::string& name) {
    particleObjects_[name] = &p;
    vtpFiles_[name] = {};
}

void timeRegistry::advanceTime() {
    currentTime_ += timeStep_;
    // For each particle, update its state
    for (auto& pair : particleObjects_) {
        const std::string& name = pair.first;
        particle* p = pair.second;
        p->update(timeStep_); // Advance the particle's state
        if (shouldWrite()) {
            // Write VTP file for this particle
            std::ostringstream vtpName;
            size_t step = vtpFiles_[name].size();
            vtpName << outputDir_ << "/" << name << "_" << std::setw(5) << std::setfill('0') << step << ".vtp";
            p->writeVtp(vtpName.str());
            vtpFiles_[name].push_back({vtpName.str(), currentTime_});
            lastWriteTime_ = currentTime_;
        }
    }
}

void timeRegistry::setCurrentTime(Foam::scalar t) {
    currentTime_ = t;
}

Foam::scalar timeRegistry::currentTime() const { return currentTime_; }
Foam::scalar timeRegistry::timeStep() const { return timeStep_; }
Foam::scalar timeRegistry::endTime() const { return endTime_; }

void timeRegistry::setWriteInterval(Foam::scalar interval) {
    writeInterval_ = interval;
}

void timeRegistry::setWriteAtStart(bool write) {
    writeAtStart_ = write;
}

void timeRegistry::setWriteAtEnd(bool write) {
    writeAtEnd_ = write;
}

Foam::scalar timeRegistry::writeInterval() const {
    return writeInterval_;
}

bool timeRegistry::shouldWrite() const {
    if (writeAtStart_ && currentTime_ <= timeStep_ + Foam::SMALL) {
        return true;
    }
    if (writeInterval_ > 0 && currentTime_ - lastWriteTime_ >= writeInterval_ - Foam::SMALL) {
        return true;
    }
    return false;
}

void timeRegistry::writePvdFile() const {
    std::ofstream pvd(pvdFileName_);
    pvd << "<?xml version=\"1.0\"?>\n";
    pvd << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    pvd << "  <Collection>\n";
    for (const auto& pair : vtpFiles_) {
        const std::string& name = pair.first;
        const auto& entries = pair.second;
        for (const auto& entry : entries) {
            std::string relFile = entry.filename;
            // Remove outputDir_ prefix if present
            if (relFile.find(outputDir_ + "/") == 0) {
                relFile = relFile.substr(outputDir_.size() + 1);
            }
            pvd << "    <DataSet timestep=\"" << entry.time << "\" group=\"" << name << "\" part=\"0\" file=\"" << relFile << "\"/>\n";
        }
    }
    pvd << "  </Collection>\n";
    pvd << "</VTKFile>\n";
    pvd.close();
}

void timeRegistry::writeTimeSeries() const {
    writePvdFile();
}

} 
