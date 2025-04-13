#include "boundary.H"
#include "foamCGALConverter.H"

namespace Bashyal
{
    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    // *** Default Constructor Implementation ***
    boundary::boundary()
        : geomObject(), // Call base class default constructor
          nef_(),       // Default construct Nef polyhedron (empty)
          patchTypes_() // Default construct patchTypes_ list (empty) if used
    {
    }

    boundary::boundary(const Foam::pointField &points, const Foam::faceList &faces)
        : geomObject(points, faces)
    {
        // Additional boundary-specific initialization can be added here if needed
    }

    boundary::boundary(const Foam::fileName &stlFile)
        : geomObject(stlFile)
    {
        // Additional boundary-specific initialization can be added here if needed
    }

    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    boundary::~boundary()
    {
    }

    void boundary::generateNefPolyhedron()
    {
        using Kernel = CGAL::Epeck;                  // Define the exact precision kernel
        using Converter = FoamCGALConverter<Kernel>; // Converter for the specified kernel
        // Convert points and faces from geomObject to CGAL Polyhedron
        auto polyhedron = Converter::toCGALPolyhedron(this->points_, this->faces_);
        // Assign the resulting Nef polyhedron to the nef_ attribute
        this->nef_ = CGAL::Nef_polyhedron_3<Kernel>(polyhedron);
    }

}
