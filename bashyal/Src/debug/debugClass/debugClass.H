#ifndef DEBUGCLASS_H
#define DEBUGCLASS_H

#include "baseClass.H"
#include "quickInclude.H"
#include "quickMesh.H"
#include "OFstream.H"
#include <string>
#include <sstream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polyhedron_3.h>
#include <vector>
#include <map>

namespace Bashyal
{
    // Define CGAL types
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

    class debugClass : public baseClass
    {
    public:
        // Constructor
        debugClass() {}

        // Destructor
        virtual ~debugClass() {}

        // Methods for writing geometry data
        virtual void write(Foam::pointField &points);
        virtual void write(Foam::List<Foam::point>& points);
        virtual void write(Foam::faceList &faces);
        void write(Foam::labelList &labels);

        virtual void write(Polyhedron& polyhedron);
        virtual void write(std::vector<Polyhedron>& polyhedrons);

        void writePolyhedronAsPolyMesh(Polyhedron& polyhedron, const std::string& folderName);
        virtual void writePoly(std::vector<Polyhedron>& polyhedrons);

        void testDisplay();

        // Methods for writing FVM matrix data
    };
}

#endif
