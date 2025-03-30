#include "debugClass.H"
#include "OFstream.H"
#include "dictionary.H"

using namespace Foam;
namespace Bashyal
{
    void debugClass::write(Foam::pointField &points)
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

    void debugClass::write(Foam::List<Foam::point> &points)
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
    void debugClass::write(Foam::faceList &faces)
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

    void debugClass::write(Foam::labelList &labels)
    {
        // Base directory path
        const Foam::fileName basePath("/usr/lib/openfoam/openfoam2312/run/debug/List<int>.txt");

        // Create an OpenFOAM OFstream object
        Foam::OFstream outFile(basePath);

        // Write each label in the labelList to the file, one per line
        for (Foam::label label : labels)
        {
            outFile << label << Foam::nl;
        }
        Foam::Info << "Label list written to " << basePath << Foam::endl;
    }

    void debugClass::testDisplay()
    {
        Foam::Info << "Here in the debugClass" << Foam::endl;
    }
}
