#include "indexedFaceSet.H"
#include "halfEdgeMesh.H"
#include "foamCGALConverter.H"
#include "implicitPlanes.H"

namespace Foam
{
    namespace particleModels
    {
        indexedFaceSet::indexedFaceSet() = default;
        indexedFaceSet::indexedFaceSet(const pointField &vertices, const faceList &faces)
            : vertices_(vertices), faces_(faces) {}

        const pointField &indexedFaceSet::vertices() const { return vertices_; }
        const faceList &indexedFaceSet::faces() const { return faces_; }

        halfEdgeMesh indexedFaceSet::toHalfEdgeMesh() const
        {
            CGALPolyhedron cgalPoly;
            cgalPoly = FoamCGALConverter<CGAL::Exact_predicates_exact_constructions_kernel>::toCGALPolyhedron(vertices_, faces_);
            return halfEdgeMesh(cgalPoly);
        }

        implicitPlanes indexedFaceSet::toImplicitPlanes() const
        {
            return toHalfEdgeMesh().toImplicitPlanes();
        }

        void indexedFaceSet::print(Ostream &os) const
        {
            os << "--- indexedFaceSet ---" << endl
               << "Vertices: " << vertices_.size() << nl
               << "Faces: " << faces_.size() << endl;
        }
    }
}