#include "halfEdgeMesh.H"
#include "indexedFaceSet.H"
#include "implicitPlanes.H"
#include <map>
#include "foamCGALConverter.H"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using CGALPolyhedron = CGAL::Polyhedron_3<Kernel>;
using CGALPoint = Kernel::Point_3;
using CGALVector = Kernel::Vector_3;
using CGALPlane_3 = Kernel::Plane_3;

// Helper: using converter for static functions
using Converter = FoamCGALConverter<Kernel>;

// Helper: check convexity (simple version)
bool isConvex(const CGALPolyhedron& poly)
{
    // For each facet, check all other vertices are on the correct side of the facet plane
    for (auto f = poly.facets_begin(); f != poly.facets_end(); ++f)
    {
        CGALPlane_3 plane = f->plane();
        for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
        {
            if (!(f->halfedge()->vertex()->point() == v->point())) {
                if (plane.oriented_side(v->point()) == CGAL::ON_POSITIVE_SIDE)
                    return false;
            }
        }
    }
    return true;
}

namespace Foam
{
    namespace particleModels
    {
        halfEdgeMesh::halfEdgeMesh() = default;
        halfEdgeMesh::halfEdgeMesh(const CGALPolyhedron &mesh)
            : polyhedron_(mesh)
        {
            if (!polyhedron_.is_valid() || !polyhedron_.is_closed())
            {
                WarningIn("halfEdgeMesh::halfEdgeMesh")
                    << "Constructed with an invalid or non-closed CGAL polyhedron."
                    << endl;
            }
        }

        const CGALPolyhedron &halfEdgeMesh::cgalMesh() const { return polyhedron_; }

        indexedFaceSet halfEdgeMesh::toIndexedFaceSet() const
        {
            pointField vertices;
            faceList faces;
            std::map<CGALPolyhedron::Vertex_const_handle, label> vertexMap;
            vertices.setSize(polyhedron_.size_of_vertices());
            label vIdx = 0;
            for (auto v = polyhedron_.vertices_begin(); v != polyhedron_.vertices_end(); ++v)
            {
                vertices[vIdx] = Converter::toFoamPoint(v->point());
                vertexMap[v] = vIdx++;
            }
            faces.setSize(polyhedron_.size_of_facets());
            label fIdx = 0;
            for (auto f = polyhedron_.facets_begin(); f != polyhedron_.facets_end(); ++f)
            {
                face &currentFace = faces[fIdx++];
                currentFace.setSize(f->size());
                label i = 0;
                auto h_circ = f->facet_begin();
                do
                {
                    currentFace[i++] = vertexMap[h_circ->vertex()];
                } while (++h_circ != f->facet_begin());
            }
            return indexedFaceSet(vertices, faces);
        }

        implicitPlanes halfEdgeMesh::toImplicitPlanes() const
        {
            if (!isConvex(polyhedron_))
            {
                FatalErrorIn("halfEdgeMesh::toImplicitPlanes")
                    << "Cannot convert a non-convex mesh to an implicit plane representation."
                    << exit(FatalError);
            }
            List<implicitPlanes::Plane> planes(polyhedron_.size_of_facets());
            point centroid = point::zero;
            for (auto v = polyhedron_.vertices_begin(); v != polyhedron_.vertices_end(); ++v)
            {
                centroid += Converter::toFoamPoint(v->point());
            }
            centroid /= polyhedron_.size_of_vertices();
            label pIdx = 0;
            for (auto f = polyhedron_.facets_begin(); f != polyhedron_.facets_end(); ++f)
            {
                CGALPlane_3 cgalPlane = f->plane();
                CGALPoint p_on_plane = f->halfedge()->vertex()->point();
                if (cgalPlane.oriented_side(Converter::toCGALPoint(centroid)) == CGAL::ON_POSITIVE_SIDE)
                {
                    cgalPlane = cgalPlane.opposite();
                }
                vector normal = Converter::toFoamVector(cgalPlane.orthogonal_vector());
                normal /= mag(normal);
                scalar distance = normal & Converter::toFoamPoint(p_on_plane);
                scalar centroid_offset = normal & centroid;
                planes[pIdx++] = {normal, distance - centroid_offset};
            }
            return implicitPlanes(planes, centroid);
        }

        void halfEdgeMesh::print(Ostream &os) const
        {
            os << "--- halfEdgeMesh ---" << endl
               << "Valid: " << (polyhedron_.is_valid() ? "Yes" : "No") << nl
               << "Vertices: " << polyhedron_.size_of_vertices() << nl
               << "Halfedges: " << polyhedron_.size_of_halfedges() << nl
               << "Faces: " << polyhedron_.size_of_facets() << endl;
        }
    }
}
