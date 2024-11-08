

Bashyal::backgroundMesh bMesh(&runTime, point(0, 0, 0), point(2.5, 0.3, 0.3), 0.1);
bMesh.developMesh();
Foam::fileName runDir = args.path();
fileName meshDir0 = "constant/polyMesh";
bMesh.writeBackgroundMesh(runDir / meshDir0);