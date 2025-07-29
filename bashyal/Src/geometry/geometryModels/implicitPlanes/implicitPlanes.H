#ifndef ImplicitPlanes_H
#define ImplicitPlanes_H

#include "vector.H"
#include "point.H"
#include "List.H"
#include "Ostream.H"

namespace Foam
{
    namespace particleModels
    {
        class indexedFaceSet; // Forward declaration

        class implicitPlanes
        {
        public:
            struct Plane
            {
                vector normal;
                scalar distance;
            };

        private:
            List<Plane> planes_;
            point centroid_;

        public:
            implicitPlanes();
            implicitPlanes(const List<Plane> &planes, const point &centroid);

            bool isInside(const point &p) const;
            indexedFaceSet toIndexedFaceSet() const;
            void print(Ostream &os) const;
        };
    }
}

#endif
