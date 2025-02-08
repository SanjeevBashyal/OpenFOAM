

// Bashyal::backgroundMesh bMesh(&runTime, point(0, 0, 0), point(2.5, 0.5, 0.5), 0.1);
Bashyal::backgroundMesh bMesh(&runTime, point(0, 0, 0), point(0.01, 0.4, 0.01), 0.001);

#include "createPatchTypes.H"

bMesh.setBoundaryPatchType(boundaryDict);

bMesh.developMesh();
Foam::fileName runDir = args.path();
fileName meshDir0 = "constant/polyMesh";
bMesh.writeBackgroundMesh(runDir / meshDir0);
fileName polyMeshThis = "polyMesh";