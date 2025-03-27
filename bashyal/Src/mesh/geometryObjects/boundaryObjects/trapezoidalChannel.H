#ifndef trapezoidalChannel_H
#define trapezoidalChannel_H

#include "geomObject.H"

namespace Bashyal
{
    class trapezoidalChannel : public geomObject
    {
    public:
        // Parameterized constructor
        trapezoidalChannel(
            Foam::scalar length,
            Foam::scalar bottomWidth,
            Foam::scalar topWidth,
            Foam::scalar height,
            Foam::scalar slope
        );
    };
}

#endif