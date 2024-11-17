#include "baseClass.H"

namespace Bashyal
{

    // Method to return a pointer to this object
    baseClass *baseClass::ptr()
    {
        return this;
    }

    // Method to return a reference to this object
    baseClass &baseClass::ref()
    {
        return *this;
    }

    // Method to return a string representation of the pointer to this object
    std::string baseClass::stringPtr() const
    {
        std::ostringstream oss;
        oss << static_cast<const void *>(this);
        return oss.str();
    }

    // Method to write pointField data to a file using OpenFOAM's OFstream
    void baseClass::writeDebug(Foam::pointField &points)
    {
        const fileName filePath("/usr/lib/openfoam/openfoam2312/run/debug/text.txt");

        // Create an OpenFOAM OFstream object
        Foam::OFstream outFile(filePath);

        // Write each point in the pointField to the file
        for (const Foam::point &pt : points)
        {
            outFile << pt << Foam::nl;
        }

        // Close the file (optional, handled by OFstream destructor)
        // outFile.close();
    }

} // namespace Bashyal
