#include "halfEdgeMesh.H"
#include "indexedFaceSet.H"
#include "implicitPlanes.H"
#include <map>

namespace Foam
{
    namespace particleModels
    {
        halfEdgeMesh::halfEdgeMesh() = default;
        halfEdgeMesh::halfEdgeMesh(const CgalPolyhedron &mesh)
            : mesh_(mesh)
        {
            if (!mesh_.is_valid() || !mesh_.is_closed())
            {
                WarningIn("halfEdgeMesh::halfEdgeMesh")
                    << "Constructed with an invalid or non-closed CGAL polyhedron."
                    << endl;
            }
        }

        const CGALPolyhedron &halfEdgeMesh::cgalMesh() const { return mesh_; }

        indexedFaceSet halfEdgeMesh::toIndexedFaceSet() const
        {
            pointField vertices;
            faceList faces;
            std::map<CgalPolyhedron::Vertex_const_handle, label> vertexMap;
            vertices.setSize(mesh_.size_of_vertices());
            label vIdx = 0;
            for (auto v = mesh_.vertices_begin(); v != mesh_.vertices_end(); ++v)
            {
                vertices[vIdx] = toFoamPoint(v->point());
                vertexMap[v] = vIdx++;
            }
            faces.setSize(mesh_.size_of_facets());
            label fIdx = 0;
            for (auto f = mesh_.facets_begin(); f != mesh_.facets_end(); ++f)
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
            if (!mesh_.is_convex())
            {
                FatalErrorIn("halfEdgeMesh::toImplicitPlanes")
                    << "Cannot convert a non-convex mesh to an implicit plane representation."
                    << exit(FatalError);
            }
            List<implicitPlanes::Plane> planes(mesh_.size_of_facets());
            point centroid = point::zero;
            for (auto v = mesh_.vertices_begin(); v != mesh_.vertices_end(); ++v)
            {
                centroid += toFoamPoint(v->point());
            }
            centroid /= mesh_.size_of_vertices();
            label pIdx = 0;
            for (auto f = mesh_.facets_begin(); f != mesh_.facets_end(); ++f)
            {
                CgalPlane_3 cgalPlane = f->plane();
                CgalPoint_3 p_on_plane = f->halfedge()->vertex()->point();
                if (cgalPlane.oriented_side(toCgalPoint(centroid)) == CGAL::ON_POSITIVE_SIDE)
                {
                    cgalPlane = cgalPlane.opposite();
                }
                vector normal = toFoamVector(cgalPlane.orthogonal_vector());
                normal /= mag(normal);
                scalar distance = normal & toFoamPoint(p_on_plane);
                scalar centroid_offset = normal & centroid;
                planes[pIdx++] = {normal, distance - centroid_offset};
            }
            return implicitPlanes(planes, centroid);
        }

        void halfEdgeMesh::print(Ostream &os) const
        {
            os << "--- halfEdgeMesh ---" << endl
               << "Valid: " << (mesh_.is_valid() ? "Yes" : "No") << nl
               << "Vertices: " << mesh_.size_of_vertices() << nl
               << "Halfedges: " << mesh_.size_of_halfedges() << nl
               << "Faces: " << mesh_.size_of_facets() << endl;
        }
    }
}