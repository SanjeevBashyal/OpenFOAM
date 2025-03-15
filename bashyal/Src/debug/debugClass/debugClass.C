#include "debugClass.H"
#include "OFstream.H"
#include "dictionary.H"

using namespace Foam;
namespace Bashyal
{
    void debugClass::writePoints(const Foam::pointField &points)
    {
        // Base directory path
        const Foam::fileName basePath("/usr/lib/openfoam/openfoam2312/run/debug/points.txt");

        // Create an OpenFOAM OFstream object
        Foam::OFstream outFile(basePath);

        // Write each point in the pointField to the file in X,Y,Z format
        for (const Foam::point &pt : points)
        {
            outFile << pt.x() << "," << pt.y() << "," << pt.z() << Foam::nl;
        }
        Foam::Info << "Done" << Foam::endl;
    }

    void debugClass::writePoints(const Foam::List<Foam::point> &points)
    {
        // Base directory path
        const Foam::fileName basePath("/usr/lib/openfoam/openfoam2312/run/debug/points.txt");

        // Create an OpenFOAM OFstream object
        Foam::OFstream outFile(basePath);

        // Write each point in the List<point> to the file in X,Y,Z format
        forAll(points, i)
        {
            const Foam::point &pt = points[i];
            outFile << pt.x() << "," << pt.y() << "," << pt.z() << Foam::nl;
        }
        Foam::Info << "Done" << Foam::endl;
    }

    // Method to write faces to a file
    void debugClass::writeFaces(const Foam::faceList &faces)
    {
        // Base directory path
        const Foam::fileName basePath("/usr/lib/openfoam/openfoam2312/run/debug/faces.txt");

        // // Merge basePath with the provided fileName
        // const Foam::fileName filePath = basePath / fileName;

        // Create an OpenFOAM OFstream object
        Foam::OFstream outFile(basePath);

        // Write each face in the faceList to the file
        for (const auto &face : faces)
        {
            for (Foam::label pointId : face)
            {
                outFile << pointId << " ";
            }
            outFile << '\n';
        }
        Foam::Info << "Done" << Foam::endl;
    }

    void debugClass::testDisplay()
    {
        Foam::Info << "Here in the debugClass" << Foam::endl;
    }
}
