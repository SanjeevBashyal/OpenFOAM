#include "backgroundBlock.H"
#include "foamCGALConverter.H" // Assumed converter header
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>

namespace Bashyal
{
    // Define CGAL types with an exact kernel for robustness
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
    typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
    typedef FoamCGALConverter<Kernel> Converter; // Adjust namespace as needed

    void backgroundBlock::intersectClosedSurfaceCGAL(
        const Foam::faceList &faces,
        const Foam::pointField &points,
        Foam::point insidePoint,
        unsigned int identifier)
    {
        // Step 1: Convert block's geometry (A) to CGAL Polyhedron
        Polyhedron A;
        try
        {
            A = Converter::toCGALPolyhedron(points_, faces_);
            if (!A.is_valid())
            {
                Foam::Warning << "Block polyhedron A is invalid." << Foam::endl;
                dead_ = true;
                ncells_ = 0;
                return;
            }
        }
        catch (const std::exception &e)
        {
            Foam::Warning << "Failed to convert block to Polyhedron: " << e.what() << Foam::endl;
            dead_ = true;
            ncells_ = 0;
            return;
        }

        // Step 2: Convert input geometry (B) to CGAL Polyhedron
        Polyhedron B;
        try
        {
            B = Converter::toCGALPolyhedron(points, faces);
            if (!B.is_valid())
            {
                Foam::Warning << "Input polyhedron B is invalid." << Foam::endl;
                dead_ = true;
                ncells_ = 0;
                return;
            }
        }
        catch (const std::exception &e)
        {
            Foam::Warning << "Failed to convert input to Polyhedron: " << e.what() << Foam::endl;
            dead_ = true;
            ncells_ = 0;
            return;
        }

        // Step 3: Convert to Nef_polyhedron_3 for subtraction
        Nef_polyhedron nefA(A);
        Nef_polyhedron nefB(B);

        // Step 4: Compute difference A - B
        Nef_polyhedron result = nefA - nefB;

        // Step 5: Handle the result
        if (result.number_of_volumes() <= 1)
        { // Only the exterior volume, i.e., result is empty
            points_.clear();
            faces_.clear();
            patches_.clear();
            dead_ = true;
            ncells_ = 0;
            edited_ = true;
            Foam::Info << "Subtraction resulted in an empty polyhedron. Block marked as dead." << Foam::endl;
            return;
        }

        // Extract the first connected component
        Polyhedron resultPoly;
        Nef_polyhedron::Volume_const_iterator ci = ++result.volumes_begin(); // Skip exterior
        result.convert_inner_shell_to_polyhedron(ci->shells_begin(), resultPoly);

        // Step 6: Convert back to OpenFOAM format
        Foam::pointField newPoints;
        Foam::faceList newFaces;
        try
        {
            Converter::toFoamPolyhedron(resultPoly, newPoints, newFaces);
        }
        catch (const std::exception &e)
        {
            Foam::Warning << "Failed to convert result back to Foam format: " << e.what() << Foam::endl;
            dead_ = true;
            ncells_ = 0;
            return;
        }

        // Step 7: Update block geometry
        points_ = newPoints;
        faces_ = newFaces;

        // Step 8: Update patches and stringPtrs, mimicking original behavior
        patches_.setSize(newFaces.size());
        forAll(newFaces, i)
        {
            // For simplicity, assign "subtracted" to all faces; adjust as needed
            // patches_[i] = "subtracted";
        }

        // Step 9: Mark as edited
        edited_ = true;

        // Step 10: Check validity
        if (points_.empty() || faces_.empty())
        {
            dead_ = true;
            ncells_ = 0;
            Foam::Info << "Resulting block has no geometry. Marked as dead." << Foam::endl;
        }
        else
        {
            // Optionally update ncells_ based on new geometry; for now, keep as is
            Foam::Info << "Block updated with " << faces_.size() << " faces after subtraction." << Foam::endl;
        }
    }

    void backgroundBlock::mapNewFacesToBoundaries(
        const Foam::pointField &newPoints,
        const Foam::faceList &newFaces,
        const Foam::pointField &pointsA,
        const Foam::faceList &facesA,
        const Foam::List<int> &patchesA,
        const Foam::pointField &pointsB,
        const Foam::faceList &facesB,
        const Foam::List<int> &patchesB,
        const Foam::label indentifierB,
        Foam::List<int> &newPatches)
    {
        // Define CGAL kernel
        typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

        // Define a tolerance for floating-point comparisons
        const double tolerance = 1e-10;

        // Helper function to compute the plane of a face from its points
        auto computePlane = [](const Foam::face &face, const Foam::pointField &points) -> CGAL::Plane_3<Kernel>
        {
            if (face.size() < 3)
            {
                throw std::invalid_argument("Face must have at least 3 points to define a plane.");
            }
            CGAL::Point_3<Kernel> p0 = Converter::toCGALPoint(points[face[0]]);
            CGAL::Point_3<Kernel> p1 = Converter::toCGALPoint(points[face[1]]);
            CGAL::Point_3<Kernel> p2 = Converter::toCGALPoint(points[face[2]]);
            return CGAL::Plane_3<Kernel>(p0, p1, p2); // Plane defined by three points
        };

        // Step 1: Compute planes for all faces of polyhedron A
        std::vector<CGAL::Plane_3<Kernel>> planesA(facesA.size());
        for (Foam::label i = 0; i < facesA.size(); ++i)
        {
            planesA[i] = computePlane(facesA[i], pointsA);
        }

        // Step 2: Compute planes for all faces of polyhedron B
        std::vector<CGAL::Plane_3<Kernel>> planesB(facesB.size());
        for (Foam::label i = 0; i < facesB.size(); ++i)
        {
            planesB[i] = computePlane(facesB[i], pointsB);
        }

        // Step 3: Resize newPatches to match the number of new faces and initialize with a default value
        newPatches.setSize(newFaces.size());
        forAll(newPatches, i)
        {
            newPatches[i] = 0; // Default value for unmatched faces
        }

        // Step 4: Map each new face to a boundary and assign patch names
        for (Foam::label i = 0; i < newFaces.size(); ++i)
        {
            const Foam::face &newFace = newFaces[i];
            CGAL::Plane_3<Kernel> newPlane = computePlane(newFace, newPoints);

            // Helper lambda to find a matching plane and return its index
            auto findMatchingPlane = [&](const std::vector<CGAL::Plane_3<Kernel>> &planes) -> Foam::label
            {
                for (Foam::label j = 0; j < planes.size(); ++j)
                {
                    const CGAL::Plane_3<Kernel> &plane = planes[j];
                    // Check if planes are parallel (cross product of normal vectors near zero)
                    if (CGAL::cross_product(plane.orthogonal_vector(), newPlane.orthogonal_vector()).squared_length() < tolerance)
                    {
                        // Check if a point from the new face lies on the plane
                        CGAL::Point_3<Kernel> p = Converter::toCGALPoint(newPoints[newFace[0]]);
                        // Compute the plane equation value at point p
                        auto plane_eq_value = plane.a() * p.x() + plane.b() * p.y() + plane.c() * p.z() + plane.d();
                        if (std::abs(CGAL::to_double(plane_eq_value)) < tolerance)
                        {
                            return j; // Return index of matching plane
                        }
                    }
                }
                return -1; // No match found
            };

            // Check for a match in polyhedron A
            Foam::label matchA = findMatchingPlane(planesA);
            if (matchA != -1)
            {
                // Assign the patch name from the matching face in A
                newPatches[i] = patchesA[matchA];
            }
            // If no match in A, check for a match in polyhedron B
            else
            {
                Foam::label matchB = findMatchingPlane(planesB);
                if (matchB != -1)
                {
                    // Assign "aggregate" for matches with B
                    newPatches[i] = indentifierB;
                }
            }
        }
    }
}
