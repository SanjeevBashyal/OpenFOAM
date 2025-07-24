#include "fvCFD.H"
#include "IndexedFaceSet.H"
#include "HalfEdgeMesh.H"
#include "ImplicitPlanes.H"

using namespace Foam;
using namespace Foam::particleModels;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "--- Particle Representation Conversion Demo ---" << nl << endl;

    // 1. Create a simple cube using the basic IndexedFaceSet
    pointField vertices =
    {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}
    };
    faceList faces =
    {
        {0, 3, 2, 1}, // bottom
        {4, 5, 6, 7}, // top
        {0, 1, 5, 4}, // front
        {2, 3, 7, 6}, // back
        {1, 2, 6, 5}, // right
        {3, 0, 4, 7}  // left
    };
    IndexedFaceSet originalParticle(vertices, faces);
    originalParticle.print(Info);
    Info << endl;

    // 2. Convert from IndexedFaceSet -> HalfEdgeMesh
    Info << "Converting to HalfEdgeMesh..." << endl;
    HalfEdgeMesh heParticle = originalParticle.toHalfEdgeMesh();
    heParticle.print(Info);
    Info << endl;

    // 3. Convert from HalfEdgeMesh -> ImplicitPlanes
    Info << "Converting to ImplicitPlanes..." << endl;
    ImplicitPlanes implicitParticle = heParticle.toImplicitPlanes();
    implicitParticle.print(Info);
    Info << endl;

    // 4. Convert back from ImplicitPlanes -> IndexedFaceSet
    Info << "Converting back to IndexedFaceSet..." << endl;
    IndexedFaceSet finalParticle = implicitParticle.toIndexedFaceSet();
    finalParticle.print(Info);
    Info << endl;

    Info << "Final particle vertices and faces:" << endl;
    Info << "Vertices: " << finalParticle.vertices() << endl;
    Info << "Faces: " << finalParticle.faces() << endl;

    Info<< "\nEnd of application.\n" << endl;
    return 0;
} 