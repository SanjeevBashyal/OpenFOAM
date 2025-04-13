#include "aggregate.H"
#include "point.H"
#include "scalar.H"
#include "pointField.H"
#include "faceList.H"
#include "List.H"

#include "Time.H"
#include "fileName.H"
#include "Random.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"

// Include the converter header if not already included via quickInclude.H
#include "foamCGALConverter.H"

namespace Bashyal
{

    // Implementation of generateNefPolyhedron
    void aggregate::generateNefPolyhedron()
    {
        using Kernel = CGAL::Epeck;
        using Converter = FoamCGALConverter<Kernel>;

        // Convert globalPoints_ and faces_ to a CGAL Polyhedron_3
        auto polyhedron = Converter::toCGALPolyhedron(globalPoints_, faces_);

        // Construct the Nef_polyhedron_3 from the Polyhedron_3 and assign to nef_
        nef_ = CGAL::Nef_polyhedron_3<Kernel>(polyhedron);
    }

}
