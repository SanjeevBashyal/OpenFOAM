#ifndef FOAM_CGAL_CONVERTER_H
#define FOAM_CGAL_CONVERTER_H

#include <vector>

// Foam headers
#include "quickInclude.H"

// CGAL headers
#include <CGAL/Surface_mesh.h>

template <typename Kernel>
class FoamCGALConverter {
public:
    // Type aliases for CGAL objects
    using CGALPoint = typename Kernel::Point_3;
    using CGALVector = typename Kernel::Vector_3;
    using CGALMesh = CGAL::Surface_mesh<CGALPoint>;

    // --- Point Conversion ---
    /** Convert a Foam point to a CGAL Point_3 */
    static CGALPoint toCGALPoint(const Foam::point& p);

    /** Convert a CGAL Point_3 to a Foam point */
    static Foam::point toFoamPoint(const CGALPoint& p);

    // --- Vector Conversion ---
    /** Convert a Foam vector to a CGAL Vector_3 */
    static CGALVector toCGALVector(const Foam::vector& v);

    /** Convert a CGAL Vector_3 to a Foam vector */
    static Foam::vector toFoamVector(const CGALVector& v);

    // --- Face Conversion ---
    /** 
     * Convert a Foam face to a vector of CGAL points 
     * @param face The Foam face (indices into points)
     * @param points The Foam point field referenced by the face
     */
    static std::vector<CGALPoint> toCGALFace(const Foam::face& face, 
                                             const Foam::pointField& points);

    /** 
     * Convert a Foam faceList to a vector of CGAL faces 
     * @param faces The Foam list of faces
     * @param points The Foam point field referenced by the faces
     */
    static std::vector<std::vector<CGALPoint>> toCGALFaces(const Foam::faceList& faces, 
                                                           const Foam::pointField& points);

    // --- Mesh Conversion ---
    /** 
     * Convert Foam pointField and faceList to a CGAL Surface_mesh 
     * @param points The Foam point field
     * @param faces The Foam list of faces
     */
    static CGALMesh toCGALMesh(const Foam::pointField& points, 
                               const Foam::faceList& faces);

    /** 
     * Convert a CGAL Surface_mesh to Foam pointField and faceList 
     * @param mesh The CGAL surface mesh
     * @param points Output Foam point field
     * @param faces Output Foam list of faces
     */
    static void toFoamMesh(const CGALMesh& mesh, 
                           Foam::pointField& points, 
                           Foam::faceList& faces);
};

#endif // FOAM_CGAL_CONVERTER_H
