#ifndef CYLINDRICALBOUNDARY_H
#define CYLINDRICALBOUNDARY_H

#include "boundary.H" // Include the base class header
#include "scalar.H"   // Include Foam::scalar type
#include "point.H"    // Include Foam::point type
#include "boundBox.H" // Include boundBox header
#include "List.H"     // Include List header
#include "label.H"    // Include Foam::label type

namespace Bashyal
{
    /**
     * @class cylindricalBoundary
     * @brief Represents a cylindrical boundary geometry derived from geomObject.
     * Axis aligned with Z. Defined by radius, center XY, Z range, and segmentation.
     * Uses initialize() to set parameters and develop() to generate geometry.
     * Stores its own bounding box and a single patch type (20) for all faces.
     */
    class cylindricalBoundary : public boundary
    {
    private:
        // --- Member Data ---

        // Geometric parameters defining the cylinder
        Foam::scalar radius_;
        Foam::point center_; // Use point for center (X, Y, reference Z - often Zmin)
        Foam::scalar zMin_;
        Foam::scalar zMax_;
        Foam::label nSegments_; // Number of segments to approximate the circular base

        Foam::List<int> patchTypes_; // List to store patch type for each face

        // --- Private Methods ---

        /**
         * @brief Generates the points, faces, and bounding box based on stored parameters.
         * Called internally by initialize().
         */
        void develop();

        /**
         * @brief Sets the patch type identifiers (all 20) for the generated faces.
         * Called internally by develop() after faces are created.
         * @param nFaces The total number of faces generated for the cylinder.
         */
        void initializePatchTypes();

    public:
        // --- Constructors ---

        /**
         * @brief Default constructor.
         * Initializes parameters to default values.
         * Call initialize() afterwards to set specific dimensions and generate geometry.
         */
        cylindricalBoundary();

        // Default destructor is sufficient
        virtual ~cylindricalBoundary() = default;

        // --- Member Functions ---

        /**
         * @brief Initializes the cylinder parameters and calls develop() to generate geometry.
         * @param radius The radius of the cylinder.
         * @param center The center point (X, Y) of the cylinder base. Z component is ignored.
         * @param zMin The minimum Z coordinate of the cylinder.
         * @param zMax The maximum Z coordinate of the cylinder.
         * @param nSegments The number of line segments to approximate the circular parts (min 3).
         */
        void initialize(
            Foam::scalar radius,
            const Foam::point &center, // Pass center point
            Foam::scalar zMin,
            Foam::scalar zMax,
            Foam::label nSegments = 32 // Default nSegments if not provided
        );

        // --- Accessors ---

        /** @brief Get cylinder radius. */
        Foam::scalar getRadius() const { return radius_; }
        /** @brief Get cylinder center point (X, Y). Z coordinate is reference. */
        const Foam::point &getCenter() const { return center_; }
        /** @brief Get minimum Z coordinate. */
        Foam::scalar getZMin() const { return zMin_; }
        /** @brief Get maximum Z coordinate. */
        Foam::scalar getZMax() const { return zMax_; }
        /** @brief Get number of segments used for approximation. */
        Foam::label getNumSegments() const { return nSegments_; }

        /**
         * @brief Get the list of patch type identifiers for each face (all should be 20).
         * @return Constant reference to the list of patch types.
         */
        const Foam::List<int> &getPatchTypes() const;

    }; // End class cylindricalBoundary

} // End namespace Bashyal

#endif // CYLINDRICALBOUNDARY_H
