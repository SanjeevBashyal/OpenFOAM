#ifndef aggregate_H
#define aggregate_H

#include "Vector.H"

namespace Bashyal
{
    class aggregate
    {
        public:
            float r1;
            float r2;
            int n;


        private:
            /* data */
        public:
            explicit aggregate(float r1, float r2, int n);
            Foam::Vector<float> getSurfacePoint(float phi, float theta);
            ~aggregate()=default;
    };

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //