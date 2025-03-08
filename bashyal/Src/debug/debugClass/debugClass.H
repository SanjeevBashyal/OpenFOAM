#ifndef DEBUGCLASS_H
#define DEBUGCLASS_H

#include "baseClass.H"
#include "quickInclude.H"
#include "quickMesh.H"
#include "OFstream.H"
#include <string>
#include <sstream>

namespace Bashyal
{
    class debugClass : public baseClass
    {
    public:
        // Constructor
        debugClass() {}

        // Destructor
        virtual ~debugClass() {}

        // Methods for writing geometry data
        virtual void writePoints(const Foam::pointField &points);
        virtual void writePoints(const Foam::List<Foam::point>& points);
        virtual void writeFaces(const Foam::faceList &faces);

        void testDisplay();

        // Methods for writing FVM matrix data
    };
}

#endif
