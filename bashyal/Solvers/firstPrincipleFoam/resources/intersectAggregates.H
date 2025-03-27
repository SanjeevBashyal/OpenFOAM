// bMesh.resetBlocks();
// bMesh.intersectCubeCGAL(a);
bMesh.developMeshPlain();

bMesh.writeBackgroundMesh(runDir / meshDir0);