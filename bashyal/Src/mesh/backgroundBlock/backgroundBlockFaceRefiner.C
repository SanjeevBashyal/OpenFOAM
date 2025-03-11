#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <map>

#include "backgroundMesh.H"

using namespace Foam;
namespace Bashyal
{
    void backgroundBlock::discretizeFaceToQuads(
        const face &inputFace,
        const pointField &inputPoints,
        pointField &outputPoints,
        faceList &outputFaces)
    {
        label nPoints = inputFace.size();

        // **Handle Invalid Cases**
        if (nPoints < 3)
        {
            Info << "Warning: Face has fewer than 3 points (" << nPoints << "). No quads generated." << endl;
            return;
        }

        // **Special Case for Quadrilaterals**
        if (nPoints == 4)
        {
            labelList newPointIndices(4);
            forAll(inputFace, i)
            {
                newPointIndices[i] = outputPoints.size();
                outputPoints.append(inputPoints[inputFace[i]]);
            }
            face quad(4);
            quad[0] = newPointIndices[0];
            quad[1] = newPointIndices[1];
            quad[2] = newPointIndices[2];
            quad[3] = newPointIndices[3];
            outputFaces.append(quad);
            return;
        }

        // **Compute the Plane**
        // Calculate the face normal and define an origin point
        vector normal = inputFace.areaNormal(inputPoints);
        normal /= mag(normal);                    // Normalize the normal vector
        point origin = inputPoints[inputFace[0]]; // Use the first point as the origin

        // Define basis vectors u and v in the plane
        vector u;
        if (mag(normal & vector(1, 0, 0)) < 0.99) // Avoid parallelism with x-axis
        {
            u = normal ^ vector(1, 0, 0);
        }
        else
        {
            u = normal ^ vector(0, 1, 0);
        }
        u /= mag(u); // Normalize u
        vector v = normal ^ u;
        v /= mag(v); // Normalize v, ensuring u, v, and normal are orthogonal

        // **Project Points to 2D**
        typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
        typedef Kernel::Point_2 Point_2;
        List<Point_2> points2D(nPoints);
        forAll(inputFace, i)
        {
            point p = inputPoints[inputFace[i]];
            scalar x = (p - origin) & u; // Dot product to get x-coordinate
            scalar y = (p - origin) & v; // Dot product to get y-coordinate
            points2D[i] = Point_2(x, y);
        }

        // **Perform Constrained Delaunay Triangulation**
        typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel> CDT;
        CDT cdt;
        std::map<CDT::Vertex_handle, label> vertexToIndex;
        List<CDT::Vertex_handle> vertexHandles(nPoints);

        // Insert points into the triangulation and map vertex handles to original indices
        forAll(points2D, i)
        {
            vertexHandles[i] = cdt.insert(points2D[i]);
            vertexToIndex[vertexHandles[i]] = inputFace[i];
        }

        // Add constraints to preserve the polygon’s boundary edges
        forAll(inputFace, i)
        {
            label iNext = (i + 1) % nPoints; // Wrap around to close the polygon
            cdt.insert_constraint(vertexHandles[i], vertexHandles[iNext]);
        }

        // **Extract Triangles**
        List<face> triangles;
        for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
        {
            face tri(3);
            for (int j = 0; j < 3; ++j)
            {
                CDT::Vertex_handle vh = fit->vertex(j);
                tri[j] = vertexToIndex[vh]; // Map back to original point indices
            }
            triangles.append(tri);
        }

        // **Decompose Each Triangle into Three Quads**
        forAll(triangles, triI)
        {
            const face &tri = triangles[triI];
            // Triangle vertices
            point p0 = inputPoints[tri[0]];
            point p1 = inputPoints[tri[1]];
            point p2 = inputPoints[tri[2]];

            // Compute midpoints of triangle edges
            point m01 = 0.5 * (p0 + p1);
            point m12 = 0.5 * (p1 + p2);
            point m20 = 0.5 * (p2 + p0);

            // Compute centroid
            point centroid = (p0 + p1 + p2) / 3.0;

            // Add points to outputPoints and assign new indices
            label p0Idx = outputPoints.size();
            outputPoints.append(p0);

            label p1Idx = outputPoints.size();
            outputPoints.append(p1);

            label p2Idx = outputPoints.size();
            outputPoints.append(p2);

            label m01Idx = outputPoints.size();
            outputPoints.append(m01);

            label m12Idx = outputPoints.size();
            outputPoints.append(m12);

            label m20Idx = outputPoints.size();
            outputPoints.append(m20);

            label centroidIdx = outputPoints.size();
            outputPoints.append(centroid);

            // Define three quads
            // Quad 1: p0, m01, centroid, m20
            face quad1(4);
            quad1[0] = p0Idx;
            quad1[1] = m01Idx;
            quad1[2] = centroidIdx;
            quad1[3] = m20Idx;
            outputFaces.append(quad1);

            // Quad 2: m01, p1, m12, centroid
            face quad2(4);
            quad2[0] = m01Idx;
            quad2[1] = p1Idx;
            quad2[2] = m12Idx;
            quad2[3] = centroidIdx;
            outputFaces.append(quad2);

            // Quad 3: m20, centroid, m12, p2
            face quad3(4);
            quad3[0] = m20Idx;
            quad3[1] = centroidIdx;
            quad3[2] = m12Idx;
            quad3[3] = p2Idx;
            outputFaces.append(quad3);
        }
    }
}
