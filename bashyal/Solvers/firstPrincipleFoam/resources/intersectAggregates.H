// bMesh.resetBlocks();
// bMesh.intersectCubeCGAL(a);
bMesh.developMesh();

bMesh.writeBackgroundMesh(runDir / meshDir0);