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

        if (boundaryNef.is_empty() || !boundaryNef.is_simple()) // is_simple is a good check
        {
            WarningIn("backgroundBlock::intersectBoundary()")
                << "Skipping intersection with invalid or empty boundary Nef for patch '"
                << domainBoundary.name() << "' in block " << blockID_ << Foam::endl;
            return;
        }

        try
        {
            CGAL::Nef_polyhedron_3<Kernel> resultNef;
            std::string opDesc;

            if (keepInside)
            {
                resultNef = nef_.intersection(boundaryNef);
                opDesc = "Intersection with Boundary '" + domainBoundary.name().str() + "'";
            }
            else // keepOutside (difference)
            {
                resultNef = nef_.difference(boundaryNef);
                opDesc = "Difference with Boundary '" + domainBoundary.name().str() + "'";
            }

            // Update the block's Nef polyhedron
            nefBase_ = resultNef;

            // Check result and update Foam geometry
            if (nef_.is_empty())
            {
                Foam::Info << "Block " << blockID_ << " became empty after " << opDesc << "." << Foam::endl;
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
            if (!nef_.is_simple())
            {
                WarningIn("backgroundBlock::intersectBoundary()")
                    << "Nef polyhedron became non-simple after " << opDesc
                    << " in block " << blockID_ << ". Result might be unreliable." << Foam::endl;
            }

            // Convert the result back to Foam geometry, passing the boundary object for patch mapping
            bool success = updateFoamGeometryFromNef(opDesc, &domainBoundary);
            if (!success)
            {
                WarningIn("backgroundBlock::intersectBoundary()")
                    << "Failed to convert result back to Foam geometry after " << opDesc
                    << " for block " << blockID_ << ". Marking as dead." << Foam::endl;
                dead_ = true;
                points_.clear();
                faces_.clear();
                patches_.clear();
                owners_.clear();
                neighbours_.clear();
                ncells_ = 0;
            }
            // Note: We don't set 'edited_' flag here, as this boundary intersection
            // is part of the initial setup phase, not subsequent modifications.
        }
        catch (const CGAL::Failure_exception &e)
        {
            WarningIn("backgroundBlock::intersectBoundary()")
                << "CGAL Error during intersection/difference with boundary '" << domainBoundary.name()
                << "' for block " << blockID_ << ": " << e.what()
                << ". Marking block as dead." << Foam::endl;
            nef_.clear(); // Make it empty
            points_.clear();
            faces_.clear();
            patches_.clear();
            owners_.clear();
            neighbours_.clear();
            ncells_ = 0;
            dead_ = true;
        }
        catch (const std::exception &e)
        {
            WarningIn("backgroundBlock::intersectBoundary()")
                << "Standard Error during intersection/difference with boundary '" << domainBoundary.name()
                << "' for block " << blockID_ << ": " << e.what()
                << ". Marking block as dead." << Foam::endl;
            nef_.clear();
            points_.clear();
            faces_.clear();
            patches_.clear();
            owners_.clear();
            neighbours_.clear();
            ncells_ = 0;
            dead_ = true;
        }
        catch (...)
        {
            WarningIn("backgroundBlock::intersectBoundary()")
                << "Unknown Error during intersection/difference with boundary '" << domainBoundary.name()
                << "' for block " << blockID_ << ". Marking block as dead." << Foam::endl;
            nef_.clear();
            points_.clear();
            faces_.clear();
            patches_.clear();
            owners_.clear();
            neighbours_.clear();
            ncells_ = 0;
            dead_ = true;
        }
    }
}
