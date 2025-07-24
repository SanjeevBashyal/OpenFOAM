#include "geomObject.H"
#include <OpenFOAM/OFstream.H>

using namespace Bashyal;
using namespace Foam;

int main(int argc, char *argv[])
{
    // Define cube vertices
    pointField vertices = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}
    };
    // Define cube faces (each face as a list of vertex indices)
    faceList faces = {
        {0, 3, 2, 1}, // bottom
        {4, 5, 6, 7}, // top
        {0, 1, 5, 4}, // front
        {2, 3, 7, 6}, // back
        {1, 2, 6, 5}, // right
        {3, 0, 4, 7}  // left
    };

    // Create geomObject
    geomObject cube(vertices, faces);

    // Write to VTP file
    cube.writeVtp("cube.vtp");

    Info << "Cube geometry written to cube.vtp" << endl;
    return 0;
} 