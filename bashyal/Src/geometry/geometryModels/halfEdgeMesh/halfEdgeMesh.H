#ifndef HALFEDGEMESH_H
#define HALFEDGEMESH_H

#include "foamCGALConverter.H"
#include "Ostream.H"

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
        explicit halfEdgeMesh(const CgalPolyhedron& mesh);

        const CgalPolyhedron& cgalMesh() const;
        indexedFaceSet toIndexedFaceSet() const;
        implicitPlanes toImplicitPlanes() const;
        void print(Ostream& os) const;
    };
}
}

#endif 