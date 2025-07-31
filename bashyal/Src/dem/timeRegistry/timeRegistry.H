#ifndef TIME_REGISTRY_H
#define TIME_REGISTRY_H

#include "particle.H"
#include <vector>
#include <string>
#include <map>

namespace Bashyal {

class timeRegistry {
    struct VtpEntry {
        std::string filename;
        Foam::scalar time;
    };
    std::map<std::string, particle*> particleObjects_; // Map object names to particle pointers
    std::map<std::string, std::vector<VtpEntry>> vtpFiles_; // Map object names to their written VTP files and times
    Foam::scalar currentTime_;
    Foam::scalar timeStep_;
    Foam::scalar endTime_;
    std::string outputDir_;
    std::string pvdFileName_;
    // Write interval attributes
    Foam::scalar writeInterval_;
    Foam::scalar lastWriteTime_;
    bool writeAtStart_;
    bool writeAtEnd_;
public:
    timeRegistry(Foam::scalar timeStep, Foam::scalar endTime, const std::string& outputDir = "VTK_output");
    void addParticle(particle& p, const std::string& name);
    void advanceTime();
    void setCurrentTime(Foam::scalar t);
    Foam::scalar currentTime() const;
    Foam::scalar timeStep() const;
    Foam::scalar endTime() const;
    // Write interval methods
    void setWriteInterval(Foam::scalar interval);
    void setWriteAtStart(bool write);
    void setWriteAtEnd(bool write);
    Foam::scalar writeInterval() const;
    bool shouldWrite() const;
    void writeTimeSeries() const;
private:
    void writePvdFile() const;
};

}

#endif 
