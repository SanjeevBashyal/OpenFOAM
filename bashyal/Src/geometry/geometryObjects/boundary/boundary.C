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

    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    boundary::~boundary()
    {
    }

    void boundary::generateNefPolyhedron()
    {
        using Kernel = CGAL::Epeck;                  // Define the exact precision kernel
        using Converter = FoamCGALConverter<Kernel>; // Converter for the specified kernel
        Foam::pointField triPoints;
        Foam::faceList triFaces;
        this->triangulateFaces(this->vertices(), this->faces(), triPoints, triFaces);

        // Convert points and faces from geomObject to CGAL Polyhedron
        auto polyhedron = Converter::toCGALPolyhedron(triPoints, triFaces);

        // Assign the resulting Nef polyhedron to the nef_ attribute
        this->nef_ = CGAL::Nef_polyhedron_3<Kernel>(polyhedron);
    }

    void boundary::triangulateFaces(
        const Foam::pointField &points, // Input points (not modified)
        const Foam::faceList &faces,    // Input faces to triangulate
        Foam::pointField &outPoints,    // Output points (same as input)
        Foam::faceList &outFaces        // Output triangulated faces
    )
    {
        // Assign input points to output (no modification occurs)
        outPoints = points;

        // Step 1: Calculate total number of triangles
        Foam::label totalTriangles = 0;
        forAll(faces, faceI)
        {
            const Foam::face &f = faces[faceI];
            if (f.size() >= 3) // Only process faces with 3 or more vertices
            {
                totalTriangles += f.size() - 2; // n-2 triangles per n-sided polygon
            }
        }

        // Step 2: Preallocate outFaces with the total number of triangles
        outFaces.setSize(totalTriangles);

        // Step 3: Triangulate each face and store in outFaces
        Foam::label triIndex = 0;
        forAll(faces, faceI)
        {
            const Foam::face &f = faces[faceI];
            if (f.size() < 3)
                continue; // Skip invalid faces
            // Fan triangulation: connect vertex 0 to vertices i and i+1
            for (Foam::label i = 1; i < f.size() - 1; ++i)
            {
                Foam::face tri(3);       // Create a triangle (face with 3 vertices)
                tri[0] = f[0];     // First vertex of the fan
                tri[1] = f[i];     // Current vertex
                tri[2] = f[i + 1]; // Next vertex
                outFaces[triIndex++] = tri;
            }
        }
    }

}
