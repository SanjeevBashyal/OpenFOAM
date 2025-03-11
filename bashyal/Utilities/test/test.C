#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/subdivision_method_3.h>
#include <vector>
#include <iostream>

// Define kernel and polyhedron types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef Kernel::Point_3 Point_3;

// Function to create a simple cube polyhedron
Polyhedron create_cube() {
    Polyhedron cube;
    std::string cube_str = 
        "OFF\n"
        "8 6 0\n"
        "-1 -1 -1\n"  // Vertex 0
        "1 -1 -1\n"   // Vertex 1
        "1 1 -1\n"    // Vertex 2
        "-1 1 -1\n"   // Vertex 3
        "-1 -1 1\n"   // Vertex 4
        "1 -1 1\n"    // Vertex 5
        "1 1 1\n"     // Vertex 6
        "-1 1 1\n"    // Vertex 7
        "4 0 1 2 3\n" // Face 0 (bottom)
        "4 4 5 6 7\n" // Face 1 (top)
        "4 0 1 5 4\n" // Face 2 (front)
        "4 2 3 7 6\n" // Face 3 (back)
        "4 1 2 6 5\n" // Face 4 (right)
        "4 0 3 7 4\n"; // Face 5 (left)

    std::istringstream iss(cube_str);
    iss >> cube;
    if (!cube.is_valid()) {
        std::cerr << "Cube creation failed!" << std::endl;
        exit(1);
    }
    return cube;
}

// Function to refine a polyhedron into hexahedrons
std::vector<Polyhedron> refine_to_hexahedrons(const Polyhedron& input_polyhedron, int refinement_level) {
    // Step 1: Convert the polyhedron to a surface mesh
    Surface_mesh mesh;
    CGAL::copy_face_graph(input_polyhedron, mesh);

    // Step 2: Subdivide the surface mesh using Loop subdivision
    for (int i = 0; i < refinement_level; ++i) {
        CGAL::Subdivision_method_3::Loop_subdivision(mesh);
    }

    // Step 3: Generate hexahedrons by "extruding" each face
    std::vector<Polyhedron> hexahedrons;
    for (auto face : mesh.faces()) {
        // Collect vertices of the current face
        std::vector<Point_3> outer_vertices;
        auto h = mesh.halfedge(face);
        auto h_start = h;
        do {
            outer_vertices.push_back(mesh.point(mesh.target(h)));
            h = mesh.next(h);
        } while (h != h_start);

        // For simplicity, assume quadrilateral faces after subdivision
        if (outer_vertices.size() != 4) {
            std::cerr << "Non-quadrilateral face encountered, skipping..." << std::endl;
            continue;
        }

        // Create inner vertices by offsetting towards the center (simplified extrusion)
        Point_3 centroid(0, 0, 0);
        for (const auto& pt : outer_vertices) {
            centroid = Point_3(
                centroid.x() + pt.x() / 4.0,
                centroid.y() + pt.y() / 4.0,
                centroid.z() + pt.z() / 4.0
            );
        }
        std::vector<Point_3> inner_vertices;
        for (const auto& pt : outer_vertices) {
            inner_vertices.push_back(Point_3(
                pt.x() + 0.5 * (centroid.x() - pt.x()),
                pt.y() + 0.5 * (centroid.y() - pt.y()),
                pt.z() + 0.5 * (centroid.z() - pt.z())
            ));
        }

        // Define the 8 vertices of the hexahedron (outer and inner)
        std::vector<Point_3> hex_vertices = {
            outer_vertices[0], outer_vertices[1], outer_vertices[2], outer_vertices[3], // Bottom face
            inner_vertices[0], inner_vertices[1], inner_vertices[2], inner_vertices[3]   // Top face
        };

        // Create a hexahedron as a polyhedron
        Polyhedron hex;
        std::ostringstream oss;
        oss << "OFF\n"
            << "8 6 0\n";
        for (const auto& v : hex_vertices) {
            oss << v.x() << " " << v.y() << " " << v.z() << "\n";
        }
        oss << "4 0 1 2 3\n"  // Bottom face
            << "4 4 5 6 7\n"  // Top face
            << "4 0 1 5 4\n"  // Front face
            << "4 2 3 7 6\n"  // Back face
            << "4 1 2 6 5\n"  // Right face
            << "4 0 3 7 4\n"; // Left face

        std::istringstream iss(oss.str());
        iss >> hex;
        if (hex.is_valid()) {
            hexahedrons.push_back(hex);
        }
    }

    return hexahedrons;
}

int main() {
    // Step 1: Create a sample polyhedron (cube)
    Polyhedron cube = create_cube();
    std::cout << "Original cube created with " << cube.size_of_facets() << " faces." << std::endl;

    // Step 2: Refine the cube into hexahedrons
    int refinement_level = 1; // Number of subdivision steps
    std::vector<Polyhedron> hexahedrons = refine_to_hexahedrons(cube, refinement_level);

    // Step 3: Output the results
    std::cout << "Number of hexahedrons generated: " << hexahedrons.size() << std::endl;
    for (size_t i = 0; i < hexahedrons.size(); ++i) {
        std::cout << "Hexahedron " << i << " has " << hexahedrons[i].size_of_facets() << " faces." << std::endl;
    }

    return 0;
}