#ifndef RECTANGULARCHANNEL_H
#define RECTANGULARCHANNEL_H

#include "boundary.H" // Include the base class header
#include "scalar.H"   // Include Foam::scalar type
#include "boundBox.H" // Include boundBox header
#include "List.H"     // *** Include List header ***

namespace Bashyal
{
    /**
     * @class rectangularChannel
     * @brief Represents a rectangular channel geometry derived from geomObject.
     * Defines the channel based on length, width, height, and slope.
     * Uses initialize() to set parameters and develop() to generate geometry.
     * Stores its own bounding box and patch type identifiers for each face.
     */
    class rectangularChannel : public boundary
    {
    private:
        // --- Member Data ---

        // Geometric parameters defining the channel
        Foam::scalar length_ = 2.0; // Default value
        Foam::scalar width_ = 0.3;  // Default value
        Foam::scalar height_ = 0.3; // Default value
        Foam::scalar slope_ = 0.01; // Default value

        // --- Private Methods ---

        /**
         * @brief Generates the points, faces, and bounding box based on stored parameters.
         * Called internally by initialize().
         */
        void develop();

        /**
         * @brief Sets the default patch type identifiers for the 6 faces.
         * Called internally by the constructor.
         */
        void initializeDefaultPatchTypes();

    public:
        // --- Constructors ---

        /**
         * @brief Default constructor.
         * Initializes parameters to default values and sets default patch types.
         * Call initialize() afterwards to set specific dimensions and generate geometry.
         */
        rectangularChannel();

        // Default destructor is sufficient
        virtual ~rectangularChannel() = default;

        // --- Member Functions ---

        /**
         * @brief Initializes the channel parameters and calls develop() to generate geometry.
         * @param length Channel length (typically along X).
         * @param width Channel width (typically along Y).
         * @param height Channel height (typically along Z at the inlet, x=0).
         * @param slope Change in Z per unit length (negative slope means outlet is lower).
         */
        void initialize(
            Foam::scalar length,
            Foam::scalar width,
            Foam::scalar height,
            Foam::scalar slope);

        // --- Accessors ---

        /** @brief Get channel length. */
        Foam::scalar getLength() const { return length_; }
        /** @brief Get channel width. */
        Foam::scalar getWidth() const { return width_; }
        /** @brief Get channel height at inlet (x=0). */
        Foam::scalar getHeight() const { return height_; }
        /** @brief Get channel slope. */
        Foam::scalar getSlope() const { return slope_; }

        /**
         * @brief Get the list of patch type identifiers for each face.
         * @return Constant reference to the list of patch types.
         */
        const Foam::List<int> &getPatchTypes() const; // *** Add getter for patchTypes_ ***

    }; // End class rectangularChannel

} // End namespace Bashyal

#endif // RECTANGULARCHANNEL_H
