// bMesh.resetBlocks();
bMesh.intersect(a);
bMesh.developMesh();

bMesh.writeBackgroundMesh(runDir / meshDir0);