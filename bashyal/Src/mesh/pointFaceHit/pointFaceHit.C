#include "pointFaceHit.H"

using namespace Foam;
namespace Bashyal
{
    // Default constructor
    pointFaceHit::pointFaceHit()
        : hit_()
    {
    }

    pointFaceHit::pointFaceHit(point pt)
        : pt_(pt)
    {
    }

    pointFaceHit::pointFaceHit(pointHit hit)
        : pt_(hit.point()), hit_(hit)
    {
    }

    pointFaceHit::pointFaceHit(pointHit hit, point startPoint)
        : pt_(hit.point()), hit_(hit), rayVertex1_(startPoint)
    {
    }

    void pointFaceHit::setBlockFaceFlag()
    {
        this->isBlockFace_ = true;
    }

    void pointFaceHit::setBlockEdgeFlag()
    {
        this->isBlockEdge_ = true;
    }

    void pointFaceHit::setBlockInsideFlag()
    {
        this->isBlockInside_ = true;
    }

    void pointFaceHit::setCutFaceFlag()
    {
        this->isCutFace_ = true;
    }

    void pointFaceHit::setCutEdgeFlag()
    {
        this->isCutEdge_ = true;
    }

    void pointFaceHit::setCutInsideFlag()
    {
        this->isCutInside_ = true;
    }

    pointFaceHit::~pointFaceHit()
    {
    }
}
