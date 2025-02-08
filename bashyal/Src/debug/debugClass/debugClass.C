#include "debugClass.H"
#include "OFstream.H"
#include "dictionary.H"

namespace Bashyal
{
    // Constructor
    debugClass::debugClass() {}

    // Destructor
    debugClass::~debugClass() {}

    void debugClass::writePoints(Foam::pointField &points)
    {
        const fileName filePath("/usr/lib/openfoam/openfoam2312/run/debug/points.txt");

        // Create an OpenFOAM OFstream object
        Foam::OFstream outFile(filePath);

        // Write each point in the pointField to the file
        for (const Foam::point &pt : points)
        {
            outFile << pt << Foam::nl;
        }

    }

    // Method to write faces to a file
    void debugClass::writeFaces(faceList& faces)
    {
        const fileName filePath("/usr/lib/openfoam/openfoam2312/run/debug/faces.txt");
        // Create an OpenFOAM OFstream object
        Foam::OFstream outFile(filePath);

        for (const auto& face : faces)
        {
            for (label pointId : face)
            {
                outFile << pointId << " ";
            }
            outFile << '\n';
        }
    }
}
