

Bashyal::backgroundMesh bMesh(&runTime, point(0, 0, 0), point(0.2, 0.2, 0.2), 0.1);
bMesh.developMesh();
Foam::fileName runDir = args.path();
fileName meshDir0 = "constant/polyMesh";
bMesh.writeBackgroundMesh(runDir / meshDir0);