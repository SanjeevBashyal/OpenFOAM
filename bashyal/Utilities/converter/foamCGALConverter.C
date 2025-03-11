#include "foamCGALConverter.H"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h> // For Epeck kernel
#include <CGAL/Polyhedron_incremental_builder_3.h> // For building Polyhedron_3

// --- Point Conversion ---
template <typename Kernel>
typename FoamCGALConverter<Kernel>::CGALPoint 
FoamCGALConverter<Kernel>::toCGALPoint(const Foam::point& p) {
    return CGALPoint(p.x(), p.y(), p.z());
}

template <typename Kernel>
Foam::point 
FoamCGALConverter<Kernel>::toFoamPoint(const CGALPoint& p) {
    double x = CGAL::to_double(p.x());
    double y = CGAL::to_double(p.y());
    double z = CGAL::to_double(p.z());
    return Foam::point(x, y, z);
}

// --- Vector Conversion ---
template <typename Kernel>
typename FoamCGALConverter<Kernel>::CGALVector 
FoamCGALConverter<Kernel>::toCGALVector(const Foam::vector& v) {
    return CGALVector(v.x(), v.y(), v.z());
}

template <typename Kernel>
Foam::vector 
FoamCGALConverter<Kernel>::toFoamVector(const CGALVector& v) {
    double x = CGAL::to_double(v.x());
    double y = CGAL::to_double(v.y());
    double z = CGAL::to_double(v.z());
    return Foam::vector(x, y, z);
}

// --- Face Conversion ---
template <typename Kernel>
std::vector<typename FoamCGALConverter<Kernel>::CGALPoint> 
FoamCGALConverter<Kernel>::toCGALFace(const Foam::face& face, 
                                      const Foam::pointField& points) {
    std::vector<CGALPoint> cgalFace;
    for (const Foam::label& idx : face) {
        cgalFace.push_back(toCGALPoint(points[idx]));
    }
    return cgalFace;
}

template <typename Kernel>
Foam::face 
FoamCGALConverter<Kernel>::toFoamFace(const std::vector<CGALPoint>& cgalFace, 
                                      Foam::pointField& points) {
    Foam::face foamFace(cgalFace.size());
    for (size_t i = 0; i < cgalFace.size(); ++i) {
        Foam::point pt = toFoamPoint(cgalFace[i]);
        bool found = false;
        for (Foam::label j = 0; j < points.size(); ++j) {
            if (Foam::mag(points[j] - pt) < 1e-10) { // Tolerance for point matching
                foamFace[i] = j;
                found = true;
                break;
            }
        }
        if (!found) {
            foamFace[i] = points.size();
            points.append(pt);
        }
    }
    return foamFace;
}

template <typename Kernel>
std::vector<std::vector<typename FoamCGALConverter<Kernel>::CGALPoint>> 
FoamCGALConverter<Kernel>::toCGALFaces(const Foam::faceList& faces, 
                                       const Foam::pointField& points) {
    std::vector<std::vector<CGALPoint>> cgalFaces;
    for (const Foam::face& face : faces) {
        cgalFaces.push_back(toCGALFace(face, points));
    }
    return cgalFaces;
}

// --- Polyhedron Conversion ---
template <typename Kernel>
typename FoamCGALConverter<Kernel>::CGALPolyhedron 
FoamCGALConverter<Kernel>::toCGALPolyhedron(const Foam::pointField& points, 
                                            const Foam::faceList& faces) {
    CGALPolyhedron polyhedron;
    // Use incremental builder to construct the polyhedron
    CGAL::Polyhedron_incremental_builder_3<typename CGALPolyhedron::HalfedgeDS> builder(polyhedron.hds(), true);
    builder.begin_surface(points.size(), faces.size());
    // Add vertices
    for (const Foam::point& p : points) {
        builder.add_vertex(toCGALPoint(p));
    }
    // Add faces
    for (const Foam::face& face : faces) {
        builder.begin_facet();
        for (Foam::label idx : face) {
            builder.add_vertex_to_facet(idx);
        }
        builder.end_facet();
    }
    builder.end_surface();
    return polyhedron;
}

template <typename Kernel>
void FoamCGALConverter<Kernel>::toFoamPolyhedron(const CGALPolyhedron& polyhedron, 
                                                 Foam::pointField& points, 
                                                 Foam::faceList& faces) {
    points.clear();
    faces.clear();
    // Add points
    for (auto v = polyhedron.vertices_begin(); v != polyhedron.vertices_end(); ++v) {
        points.push_back(toFoamPoint(v->point()));
    }
    // Add faces
    for (auto f = polyhedron.facets_begin(); f != polyhedron.facets_end(); ++f) {
        Foam::face foamFace;
        auto circ = f->facet_begin();
        do {
            // Assign vertex index based on its position in the vertex list
            size_t idx = std::distance(polyhedron.vertices_begin(), circ->vertex());
            foamFace.push_back(idx);
        } while (++circ != f->facet_begin());
        faces.push_back(foamFace);
    }
}

// --- Mesh Conversion ---
template <typename Kernel>
typename FoamCGALConverter<Kernel>::CGALMesh 
FoamCGALConverter<Kernel>::toCGALMesh(const Foam::pointField& points, 
                                      const Foam::faceList& faces) {
    CGALMesh mesh;
    std::vector<typename CGALMesh::Vertex_index> vertex_indices;
    for (const Foam::point& p : points) {
        vertex_indices.push_back(mesh.add_vertex(toCGALPoint(p)));
    }
    for (const Foam::face& face : faces) {
        std::vector<typename CGALMesh::Vertex_index> face_vertices;
        for (const Foam::label& idx : face) {
            face_vertices.push_back(vertex_indices[idx]);
        }
        mesh.add_face(face_vertices);
    }
    return mesh;
}

template <typename Kernel>
void FoamCGALConverter<Kernel>::toFoamMesh(const CGALMesh& mesh, 
                                           Foam::pointField& points, 
                                           Foam::faceList& faces) {
    points.clear();
    faces.clear();
    for (const auto& v : mesh.vertices()) {
        points.push_back(toFoamPoint(mesh.point(v)));
    }
    for (const auto& f : mesh.faces()) {
        Foam::face of_face;
        for (const auto& v : mesh.vertices_around_face(mesh.halfedge(f))) {
            of_face.push_back(v.idx());
        }
        faces.push_back(of_face);
    }
}

// Explicit instantiation for a common kernel
#include <CGAL/Simple_cartesian.h>
template class FoamCGALConverter<CGAL::Simple_cartesian<double>>;

// Explicit instantiation for CGAL::Epeck
template class FoamCGALConverter<CGAL::Epeck>;
