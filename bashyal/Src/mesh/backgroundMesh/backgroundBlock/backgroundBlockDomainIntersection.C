#include "backgroundBlock.H"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h>

namespace Bashyal
{
    // Define CGAL types with an exact kernel for robustness
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
    typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;

    void backgroundBlock::intersectBoundary(const boundary &domainBoundary, bool keepInside)
    {
        if (dead_)
            return; // Don't process dead blocks

        const CGAL::Nef_polyhedron_3<Kernel> &boundaryNef = domainBoundary.nef_;

        CGAL::Nef_polyhedron_3<Kernel> resultNef;
        std::string opDesc;

        if (keepInside)
        {
            resultNef = nef_.intersection(boundaryNef);
        }
        else // keepOutside (difference)
        {
            resultNef = nef_.difference(boundaryNef);
        }

        // Update the block's Nef polyhedron
        nefBase_ = resultNef;

        intersectedBoundaries_.append(&domainBoundary);

        // Check result and update Foam geometry
        if (nef_.is_empty())
        {
            dead_ = true;
            // Clear Foam geometry explicitly
            points_.clear();
            faces_.clear();
            patches_.clear();
            owners_.clear();
            neighbours_.clear();
            ncells_ = 0;
            return; // Exit early
        }
    }
}