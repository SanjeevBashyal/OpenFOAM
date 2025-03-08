#ifndef pointFaceHit_H
#define pointFaceHit_H

#include "quickInclude.H"
#include "quickMesh.H"

namespace Bashyal
{
    class pointFaceHit
    {
    private:
        /* data */

    public:
        Foam::point pt_;
        Foam::pointHit hit_;
        Foam::point rayVertex1_;
        Foam::point rayVertex2_;
        Foam::point nextPoint_;

        Foam::label blockFace1_ = -1;
        Foam::label blockFace2_ = -1;
        Foam::label cutFace1_ = -1;
        Foam::label cutFace2_ = -1;

        bool isBlockFace_ = false;
        bool isBlockEdge_ = false;
        bool isBlockInside_ = false;

        bool isCutFace_ = false;
        bool isCutEdge_ = false;
        bool isCutInside_ = false;

    public:
        pointFaceHit();
        pointFaceHit(Foam::point pt);
        pointFaceHit(Foam::pointHit hit);
        pointFaceHit(Foam::pointHit hit, Foam::point startPoint);

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
