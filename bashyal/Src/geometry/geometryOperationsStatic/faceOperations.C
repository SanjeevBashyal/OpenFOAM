#include "faceOperations.H"

namespace Bashyal
{
    void faceOperations::triangulateFaces(
        const Foam::pointField& points, // Input points (not modified)
        const Foam::faceList& faces,    // Input faces to triangulate
        Foam::pointField& outPoints,    // Output points (same as input)
        Foam::faceList& outFaces        // Output triangulated faces
    )
    {
        // Assign input points to output (no modification occurs)
        outPoints = points;

        // Step 1: Calculate total number of triangles
        Foam::label totalTriangles = 0;
        forAll(faces, faceI)
        {
            const Foam::face& f = faces[faceI];
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
            const Foam::face& f = faces[faceI];
            if (f.size() < 3)
                continue; // Skip invalid faces
            // Fan triangulation: connect vertex 0 to vertices i and i+1
            for (Foam::label i = 1; i < f.size() - 1; ++i)
            {
                Foam::face tri(3);       // Create a triangle (face with 3 vertices)
                tri[0] = f[0];           // First vertex of the fan
                tri[1] = f[i];           // Current vertex
                tri[2] = f[i + 1];       // Next vertex
                outFaces[triIndex++] = tri;
            }
        }
    }


    void faceOperations::writeToSTL(
        const Foam::pointField& points,
        const Foam::faceList& faces,
        const Foam::fileName& stlFile
    )
    {
        // Create a list of labelled triangles from the faces
        Foam::List<Foam::labelledTri> triangles;

        // Iterate over each face in the provided faces list
        for (const auto& face : faces)
        {
            if (face.size() == 3)
            {
                // Triangular faces are added directly as labelled triangles
                triangles.append(Foam::labelledTri(face[0], face[1], face[2]));
            }
            else if (face.size() == 4)
            {
                // Split quadrilateral faces into two triangles
                triangles.append(Foam::labelledTri(face[0], face[1], face[2]));
                triangles.append(Foam::labelledTri(face[0], face[2], face[3]));
            }
            else
            {
                // Throw an error for unsupported face sizes
                FatalErrorInFunction
                    << "Unsupported face size: " << face.size() << ". Only triangles (3) and quads (4) are supported."
                    << abort(Foam::FatalError);
            }
        }

        // Create a triSurface object from the triangles and points
        Foam::triSurface surface(triangles, points);

        // Write the surface to the specified STL file
        surface.write(stlFile);
    }
}
