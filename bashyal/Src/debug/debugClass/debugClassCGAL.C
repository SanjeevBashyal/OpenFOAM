#include "debugClass.H"

namespace Bashyal
{
    void debugClass::write(Polyhedron &polyhedron)
    {
        const Foam::fileName basePath("/usr/lib/openfoam/openfoam2312/run/debug/CGALpolyhedron.txt");
        Foam::OFstream outFile(basePath);

        std::map<const Polyhedron::Vertex *, int> vertexIndexMap;
        int index = 0;

        outFile << "Vertices: " << polyhedron.size_of_vertices() << "\n";

        // Write vertices
        for (auto vit = polyhedron.vertices_begin(); vit != polyhedron.vertices_end(); ++vit)
        {
            const auto &point = vit->point();
            outFile << CGAL::to_double(point.x()) << " " << CGAL::to_double(point.y()) << " " << CGAL::to_double(point.z()) << "\n";
            vertexIndexMap[&(*vit)] = index++;
        }

        outFile << "Faces: " << polyhedron.size_of_facets() << "\n";
        // Write faces (unchanged, as indices are unaffected)
        for (auto fit = polyhedron.facets_begin(); fit != polyhedron.facets_end(); ++fit)
        {
            int count = 0;
            auto circ = fit->facet_begin();
            do
            {
                outFile << vertexIndexMap[&(*circ->vertex())] << " ";
                ++circ;
                count++;
            } while (circ != fit->facet_begin());
            outFile << " " << count;
            outFile << "\n";
        }
    }

    void debugClass::write(std::vector<Polyhedron> &polyhedrons)
    {
        const Foam::fileName basePath("/usr/lib/openfoam/openfoam2312/run/debug/CGALpolyhedrons.txt");
        Foam::OFstream outFile(basePath);
        outFile << "Number of Polyhedrons: " << polyhedrons.size() << "\n";
        int polyIndex = 0;
        for (const auto &poly : polyhedrons)
        {
            outFile << "Polyhedron " << polyIndex++ << "\n";
            outFile << "Vertices: " << poly.size_of_vertices() << "\n";
            std::map<const Polyhedron::Vertex *, int> vertexIndexMap;
            int index = 0;
            for (auto vit = poly.vertices_begin(); vit != poly.vertices_end(); ++vit)
            {
                const auto &point = vit->point();
                outFile << CGAL::to_double(point.x()) << " " << CGAL::to_double(point.y()) << " " << CGAL::to_double(point.z()) << "\n";
                vertexIndexMap[&(*vit)] = index++;
            }

            outFile << "Faces: " << poly.size_of_facets() << "\n";
            for (auto fit = poly.facets_begin(); fit != poly.facets_end(); ++fit)
            {
                int count = 0;
                auto circ = fit->facet_begin();
                do
                {
                    outFile << vertexIndexMap[&(*circ->vertex())] << " ";
                    ++circ;
                    count++;
                } while (circ != fit->facet_begin());
                outFile << " " << count;
                outFile << "\n";
            }
        }
        Foam::Info << "Done writing polyhedrons" << Foam::endl;
    }
}
