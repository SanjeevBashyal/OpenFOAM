bMesh.resetBlocks();
bMesh.intersectCube(a);
bMesh.developMesh();

bMesh.writeBackgroundMesh(runDir / meshDir0);