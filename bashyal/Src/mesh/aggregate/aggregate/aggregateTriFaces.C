#include "aggregate.H"

namespace Bashyal
{
    // Create triangulated surface for STL output
    void aggregate::createTriSurface()
    {
        Foam::List<Foam::labelledTri> triangles;
        for (const auto& face : this->faces_)
        {
            if (face.size() == 4)
            {
                // Split quads into two triangles
                triangles.append(Foam::labelledTri(face[0], face[1], face[2]));
                triangles.append(Foam::labelledTri(face[0], face[2], face[3]));
            }
            else
            {
                // Triangular faces remain as is
                triangles.append(Foam::labelledTri(face[0], face[1], face[2]));
            }
        }

        Foam::triSurface surface(triangles, this->points_);
        Foam::fileName outputFile("/usr/lib/openfoam/openfoam2312/run/debug/aggregate.stl"); // Second argument is the output file name
        surface.write(outputFile);
        this->surface_ = surface;
    }
}