// backgroundBlock.C (or relevant implementation file)

#include "backgroundBlock.H"
#include "aggregate.H"         // Include aggregate header for agg.points(), etc.
#include "foamCGALConverter.H" // Assumed converter header

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h> // Potentially needed by Nef/decomposition

#include <vector>
#include <exception>
#include <cmath>     // For std::abs
#include <algorithm> // For std::sort

namespace Bashyal // Or your appropriate namespace
{
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
    typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
    typedef FoamCGALConverter<Kernel> Converter;

    // Merged function: Subtracts aggregate, decomposes, and rebuilds Foam geometry
    void backgroundBlock::subtractAggregate(
        const aggregate &agg) 
    {
        // --- Pre-checks ---
        if (dead_)
        {
            // Already dead, nothing to do
            return;
        }

        const Nef_polyhedron &aggNef = agg.nef_; 

        // --- Perform Nef Subtraction ---
        Nef_polyhedron resultNef = nef_ - aggNef; // Or use operator-

        nef_ = resultNef; // Update the block's Nef polyhedron

        // --- Handle Empty Result ---
        if (resultNef.number_of_volumes() <= 1)
        {
            Foam::Info << "Block " << blockID_ << " became empty after subtracting aggregate "
                       << agg.identifier_ << "." << Foam::endl; // Use agg.identifier_ if available
            dead_ = true;
            points_.clear();
            faces_.clear();
            patches_.clear();
            owners_.clear();
            neighbours_.clear();
            ncells_ = 0;
            multiple_ = false; // Reset multiple flag
            nef_ = resultNef;  // Update nef_ to the empty state
            edited_ = true;    // Mark as edited
            return;
        }
    }
}
