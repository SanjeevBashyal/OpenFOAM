#include "pointFaceHit.H"

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
        : pt_(hit.point()), hit_(hit), previousPoint_(startPoint)
    {
    }

    void pointFaceHit::setFaceFlag()
    {
        this->isFace_ = true;
    }

    void pointFaceHit::setEdgeFlag()
    {
        this->isEdge_ = true;
    }

    void pointFaceHit::setInsideFlag()
    {
        this->isInside_ = true;
    }

    pointFaceHit::~pointFaceHit()
    {
    }
}
