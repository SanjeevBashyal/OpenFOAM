#include "foamCGALConverter.H"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h> // For Epeck kernel

template <typename Kernel>
typename FoamCGALConverter<Kernel>::CGALPoint 
FoamCGALConverter<Kernel>::toCGALPoint(const Foam::point& p) {
    return CGALPoint(p.x(), p.y(), p.z());
}

template <typename Kernel>
Foam::point 
FoamCGALConverter<Kernel>::toFoamPoint(const CGALPoint& p) {
    // Convert exact coordinates to double
    double x = CGAL::to_double(p.x());
    double y = CGAL::to_double(p.y());
    double z = CGAL::to_double(p.z());
    return Foam::point(x, y, z);
}

template <typename Kernel>
typename FoamCGALConverter<Kernel>::CGALVector 
FoamCGALConverter<Kernel>::toCGALVector(const Foam::vector& v) {
    return CGALVector(v.x(), v.y(), v.z());
}

template <typename Kernel>
Foam::vector 
FoamCGALConverter<Kernel>::toFoamVector(const CGALVector& v) {
    // Convert exact coordinates to double
    double x = CGAL::to_double(v.x());
    double y = CGAL::to_double(v.y());
    double z = CGAL::to_double(v.z());
    return Foam::vector(x, y, z);
}

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
std::vector<std::vector<typename FoamCGALConverter<Kernel>::CGALPoint>> 
FoamCGALConverter<Kernel>::toCGALFaces(const Foam::faceList& faces, 
                                       const Foam::pointField& points) {
    std::vector<std::vector<CGALPoint>> cgalFaces;
    for (const Foam::face& face : faces) {
        cgalFaces.push_back(toCGALFace(face, points));
    }
    return cgalFaces;
}

template <typename Kernel>
typename FoamCGALConverter<Kernel>::CGALMesh 
FoamCGALConverter<Kernel>::toCGALMesh(const Foam::pointField& points, 
                                      const Foam::faceList& faces) {
    CGALMesh mesh;
    std::vector<typename CGALMesh::Vertex_index> vertex_indices;
    // Add vertices
    for (const Foam::point& p : points) {
        vertex_indices.push_back(mesh.add_vertex(toCGALPoint(p)));
    }
    // Add faces
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
    // Add points
    for (const auto& v : mesh.vertices()) {
        points.push_back(toFoamPoint(mesh.point(v)));
    }
    // Add faces
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