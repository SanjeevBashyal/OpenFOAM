#ifndef pointFaceHit_H
#define pointFaceHit_H

#include "quickInclude.H"
#include "quickMesh.H"

using namespace Foam;

namespace Bashyal
{

    class pointFaceHit
    {
    private:
        /* data */

    public:
        point pt_;
        pointHit hit_;
        point previousPoint_;
        point nextPoint_;
        label face1_ = -1;
        label face2_ = -1;

        bool isFace_ = false;
        bool isEdge_ = false;
        bool isInside_ = false;

    public:
        pointFaceHit();
        pointFaceHit(point pt);
        pointFaceHit(pointHit hit);
        pointFaceHit(pointHit hit, point startPoint);

        void setFaceFlag();
        void setEdgeFlag();
        void setInsideFlag();

        ~pointFaceHit();
    };

}

#endif
