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
        point rayVertex1_;
        point rayVertex2_;
        point nextPoint_;

        label blockFace1_ = -1;
        label blockFace2_ = -1;
        label cutFace1_ = -1;
        label cutFace2_ = -1;

        bool isBlockFace_ = false;
        bool isBlockEdge_ = false;
        bool isBlockInside_ = false;

        bool isCutFace_ = false;
        bool isCutEdge_ = false;
        bool isCutInside_ = false;

    public:
        pointFaceHit();
        pointFaceHit(point pt);
        pointFaceHit(pointHit hit);
        pointFaceHit(pointHit hit, point startPoint);

        void setBlockFaceFlag();
        void setBlockEdgeFlag();
        void setBlockInsideFlag();

        void setCutFaceFlag();
        void setCutEdgeFlag();
        void setCutInsideFlag();

        ~pointFaceHit();
    };

}

#endif
