

// Bashyal::backgroundMesh bMesh(&runTime, point(0, 0.1, 0.1), point(0.1, 0.3, 0.2), 0.1);
Bashyal::backgroundMesh bMesh(&runTime, point(0, 0, 0), point(0.3, 0.3, 0.3), 0.1);

#include "createPatchTypes.H"

bMesh.setBoundaryPatchType(boundaryDict);

// bMesh.developMesh();
Foam::fileName runDir = args.path();
fileName meshDir0 = "constant/polyMesh";
// bMesh.writeBackgroundMesh(runDir / meshDir0);