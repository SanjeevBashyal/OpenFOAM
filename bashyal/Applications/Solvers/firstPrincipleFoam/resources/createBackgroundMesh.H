

Bashyal::backgroundMesh bMesh(&runTime);

#include "createPatchTypes.H"
bMesh.setBoundaryPatchType(boundaryDict);

// bMesh.developMesh();
Foam::fileName runDir = args.path();
fileName meshDir0 = "constant/polyMesh";
// bMesh.writeBackgroundMesh(runDir / meshDir0);