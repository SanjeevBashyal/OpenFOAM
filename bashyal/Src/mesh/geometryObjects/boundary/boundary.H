#ifndef boundary_H
#define boundary_H

#include "quickInclude.H"
#include "geomObject.H"
#include "faceOperations.H"

#include <CGAL/Nef_polyhedron_3.h>

namespace Bashyal
{
    class boundary : public geomObject
    {
    private:
        // Add boundary-specific private members here if needed

    public:
        // Existing attribute
        CGAL::Nef_polyhedron_3<CGAL::Epeck> nef_;
        
        // *** Add patch type list attribute ***
        Foam::List<int> patchTypes_;

    public:
        // Constructors
        boundary();
        boundary(const Foam::pointField &points, const Foam::faceList &faces);
        boundary(const Foam::fileName &stlFile);

        // Destructor
        ~boundary();

        // New method to generate the Nef polyhedron
        void generateNefPolyhedron();
    };
}

#endif
