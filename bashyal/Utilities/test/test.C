#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_complex_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <utility>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Mesh_3::Mesh_polyhedron_items<int>> Polyhedron;
typedef CGAL::Polyhedral_complex_mesh_domain_3<Kernel, Polyhedron> Mesh_domain;
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

int main() {
    // Load outer surface (S1)
    Polyhedron S1;
    std::ifstream in1("outer.off");
    if (!in1 || !(in1 >> S1) || !S1.is_closed()) {
        std::cerr << "Error: Unable to load outer surface or it is not closed." << std::endl;
        return 1;
    }
    if (!S1.is_valid()) {
        std::cerr << "Error: Outer surface is not a valid polyhedron." << std::endl;
        return 1;
    }

    // Load inner surface (S2)
    Polyhedron S2;
    std::ifstream in2("inner.off");
    if (!in2 || !(in2 >> S2) || !S2.is_closed()) {
        std::cerr << "Error: Unable to load inner surface or it is not closed." << std::endl;
        return 1;
    }
    if (!S2.is_valid()) {
        std::cerr << "Error: Inner surface is not a valid polyhedron." << std::endl;
        return 1;
    }

    std::cout << "Polyhedra loaded and validated successfully." << std::endl;

    // Define polyhedra and subdomain pairs
    std::vector<Polyhedron> polyhedra = {S1, S2};
    // Corrected subdomain pairs: S1 (outer) {0, 1}, S2 (inner) {1, -1}
    std::vector<std::pair<int, int>> subdomain_pairs = {{0, 1}, {1, -1}};
    // S1: outside = 0, inside = 1
    // S2: outside (between S1 and S2) = 1, inside = -1
    // Result:
    // - Outside S1: 0
    // - Between S1 and S2: 1 (from S1) + 1 (from S2 outside) = 2 (meshed)
    // - Inside S2: 1 (from S1) + (-1) (from S2 inside) = 0 (not meshed)

    // Create mesh domain
    Mesh_domain domain(polyhedra.begin(), polyhedra.end(), subdomain_pairs.begin(), subdomain_pairs.end());

    // Define mesh criteria
    Mesh_criteria criteria(
        CGAL::parameters::facet_angle = 25,
        CGAL::parameters::facet_size = 0.1,
        CGAL::parameters::facet_distance = 0.01,
        CGAL::parameters::cell_radius_edge_ratio = 3,
        CGAL::parameters::cell_size = 0.1
    );

    // Generate the mesh
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

    // Check if tetrahedra exist
    if (c3t3.number_of_cells() == 0) {
        std::cerr << "Error: No tetrahedra generated in the mesh!" << std::endl;
        return 1;
    }
    std::cout << "Number of tetrahedra: " << c3t3.number_of_cells() << std::endl;

    // Output to file
    std::ofstream medit_file("output.mesh");
    if (!medit_file) {
        std::cerr << "Error: Unable to open output file." << std::endl;
        return 1;
    }
    c3t3.output_to_medit(medit_file);
    std::cout << "Mesh generated and saved to 'output.mesh'." << std::endl;

    return 0;
}