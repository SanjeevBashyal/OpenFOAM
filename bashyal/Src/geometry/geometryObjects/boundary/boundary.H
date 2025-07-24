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

        Foam::word name_; // Name of the boundary

    public:
        // Constructors
        boundary();
        boundary(const Foam::pointField &points, const Foam::faceList &faces);

        // Destructor
        ~boundary();

        // New method to generate the Nef polyhedron
        void generateNefPolyhedron();
        void triangulateFaces(
            const Foam::pointField &points,
            const Foam::faceList &faces,
            Foam::pointField &outPoints,
            Foam::faceList &outFaces);

        Foam::List<int> patchTypes() const { return patchTypes_; }                          // Getter for patch types
        void setPatchTypes(const Foam::List<int> &patchTypes) { patchTypes_ = patchTypes; } // Setter for patch types

        Foam::word name() const { return name_; } // Getter for name
        void setName(const Foam::word &name) { name_ = name; } // Setter for name
    };
}

#endif
