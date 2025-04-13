#include "backgroundBlock.H"
#include "foamCGALConverter.H"

namespace Bashyal
{
    void backgroundBlock::generateNefPolyhedron()
    {
        using Kernel = CGAL::Epeck;                  // Define the exact precision kernel
        using Converter = FoamCGALConverter<Kernel>; // Converter for the specified kernel
        // Convert points and faces from geomObject to CGAL Polyhedron
        auto polyhedron = Converter::toCGALPolyhedron(this->points_, this->faces_);
        // Assign the resulting Nef polyhedron to the nef_ attribute
        this->nefBase_ = CGAL::Nef_polyhedron_3<Kernel>(polyhedron);
    }
}