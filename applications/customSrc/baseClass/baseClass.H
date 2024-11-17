#ifndef BASECLASS_H
#define BASECLASS_H

#include "pointField.H"
#include "OFstream.H"
#include <string>
#include <sstream>

using namespace Foam;

namespace Bashyal
{

class baseClass
{
public:
    // Constructor
    baseClass() {}

    // Destructor
    virtual ~baseClass() {}

    // Method to return a pointer to this object
    baseClass* ptr();

    // Method to return a reference to this object
    baseClass& ref();

    // Method to return a string representation of the pointer to this object
    std::string stringPtr() const;

    // Method to write pointField data to a file
    virtual void writeDebug(Foam::pointField& points);
};

} // namespace Bashyal

#endif // BASECLASS_H
