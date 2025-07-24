#ifndef faceOperations_H
#define faceOperations_H

#include "quickInclude.H"
#include "triSurface.H"       // For triSurface functionality
#include "labelledTri.H"      // For labelledTri type

namespace Bashyal
{
    class faceOperations
    {
    public:
        // Static method to triangulate faces
        static void triangulateFaces(
            const Foam::pointField &points, // Input points (not modified)
            const Foam::faceList &faces,    // Input faces to triangulate
            Foam::pointField &outPoints,    // Output points (same as input)
            Foam::faceList &outFaces        // Output triangulated faces
        );

        // Static method to write points and faces to an STL file
        static void writeToSTL(
            const Foam::pointField &points,
            const Foam::faceList &faces,
            const Foam::fileName &stlFile);
    };
}

#endif
