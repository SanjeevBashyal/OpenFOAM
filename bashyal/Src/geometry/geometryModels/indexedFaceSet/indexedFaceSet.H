#ifndef indexedFaceSet_H
#define indexedFaceSet_H

#include "pointField.H"
#include "faceList.H"
#include "Ostream.H"
#include "implicitPlanes.H"

namespace Foam
{
    namespace particleModels
    {
        class halfEdgeMesh; // Forward declaration
        class implicitPlanes; // Forward declaration

        class indexedFaceSet
        {
            pointField vertices_;
            faceList faces_;

        public:
            indexedFaceSet();
            indexedFaceSet(const pointField &vertices, const faceList &faces);

            const pointField &vertices() const;
            const faceList &faces() const;

            halfEdgeMesh toHalfEdgeMesh() const;
            implicitPlanes toImplicitPlanes() const;
            void print(Ostream &os) const;
        };
    }
}

#endif
