#ifndef DEBUGCLASS_H
#define DEBUGCLASS_H

#include "quickInclude.H"
#include "quickMesh.H"
#include "OFstream.H"
#include <string>
#include <sstream>

using namespace Foam;
namespace Bashyal
{
    class debugClass
    {
    public:
        // Constructor and Destructor
        debugClass();
        ~debugClass();

        // Methods for writing geometry data
        void writePoints(pointField& points);
        void writeFaces(faceList& faces);

        // Methods for writing FVM matrix data
    };
}

#endif
