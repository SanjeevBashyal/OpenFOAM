// patchTypes.H
#ifndef PATCH_TYPES_H
#define PATCH_TYPES_H

namespace Bashyal
{
    enum patchType
    {
        YZ_Xmin = -1, // Left face (minimum X)
        YZ_Xmax = 1,  // Right face (maximum X)
        XZ_Ymin = -2, // Front face (minimum Y)
        XZ_Ymax = 2,   // Back face (maximum Y)
        XY_Zmin = -3, // Bottom face (minimum Z)
        XY_Zmax = 3   // Top face (maximum Z)
    };
}

#endif // PATCH_TYPES_H
