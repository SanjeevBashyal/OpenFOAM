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
        virtual void writePoints(const Foam::pointField &points);
        virtual void writePoints(const Foam::List<Foam::point>& points);
        virtual void writeFaces(const Foam::faceList &faces);
        void writeLabelList(const Foam::labelList &labels);

        virtual void writePolyhedron(const Polyhedron& polyhedron);
        virtual void writePolyhedrons(const std::vector<Polyhedron>& polyhedrons);

        void writePolyhedronAsPolyMesh(const Polyhedron& polyhedron, const std::string& folderName);
        virtual void writePolyhedronsAsPolyMesh(const std::vector<Polyhedron>& polyhedrons);

        void testDisplay();

        // Methods for writing FVM matrix data
    };
}

#endif
