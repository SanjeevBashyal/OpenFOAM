Foam::dictionary boundaryDict;

boundaryDict.add("YZ-Xmin", "wall");
boundaryDict.add("YZ-Xmax", "wall");
boundaryDict.add("XZ-Ymin", "wall");
boundaryDict.add("XZ-Ymax", "patch");
boundaryDict.add("XY-Zmin", "wall");
boundaryDict.add("XY-Zmax", "wall");
boundaryDict.add("aggregate", "wall");
boundaryDict.add("frontAndBack", "empty");

// boundaryDict.add("YZ-Xmin", "patch");
// boundaryDict.add("YZ-Xmax", "patch");
// boundaryDict.add("XZ-Ymin", "patch");
// boundaryDict.add("XZ-Ymax", "patch");
// boundaryDict.add("XY-Zmin", "patch");
// boundaryDict.add("XY-Zmax", "patch");
// boundaryDict.add("aggregate", "patch");
// boundaryDict.add("frontAndBack", "empty");