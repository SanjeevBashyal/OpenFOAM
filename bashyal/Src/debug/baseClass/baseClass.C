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
    

} // namespace Bashyal
