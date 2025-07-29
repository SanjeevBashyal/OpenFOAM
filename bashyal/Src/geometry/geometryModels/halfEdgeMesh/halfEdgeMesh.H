#ifndef HALFEDGEMESH_H
#define HALFEDGEMESH_H

#include "foamCGALConverter.H"
#include "Ostream.H"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// Add typedefs for Kernel and CGAL types
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using CGALPolyhedron = CGAL::Polyhedron_3<Kernel>;
using CGALPoint = Kernel::Point_3;
using CGALVector = Kernel::Vector_3;
using CGALPlane_3 = Kernel::Plane_3;

namespace Foam
{
    namespace particleModels
    {
        class indexedFaceSet;
        class implicitPlanes;

        class halfEdgeMesh
        {
            CGALPolyhedron polyhedron_;

        public:
            halfEdgeMesh();
            explicit halfEdgeMesh(const CGALPolyhedron &mesh);

            const CGALPolyhedron &cgalMesh() const;
            indexedFaceSet toIndexedFaceSet() const;
            implicitPlanes toImplicitPlanes() const;
            void print(Ostream &os) const;
        };
    }
}

#endif
