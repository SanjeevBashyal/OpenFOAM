#ifndef boundary_H
#define boundary_H

#include "quickInclude.H"
#include "geomObject.H"
#include "faceOperations.H"

namespace Bashyal
{
    class boundary : public geomObject
    {
    private:
        // Add boundary-specific private members here if needed

    public:
        // Constructors
        boundary(const Foam::pointField& points, const Foam::faceList& faces);
        boundary(const Foam::fileName& stlFile);

        // Destructor
        ~boundary();
    };
}

#endif
