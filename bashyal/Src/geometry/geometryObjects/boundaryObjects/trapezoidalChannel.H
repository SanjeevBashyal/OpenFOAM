#ifndef TRAPEZOIDALCHANNEL_H
#define TRAPEZOIDALCHANNEL_H

#include "boundary.H"
#include "scalar.H"
#include "boundBox.H"
#include "List.H"

namespace Bashyal
{
    /**
     * @class trapezoidalChannel
     * @brief Represents a trapezoidal channel geometry derived from boundary.
     * Defines the channel based on length, bottom width, top width, height, and slope.
     * Uses initialize() to set parameters and develop() to generate geometry.
     */
    class trapezoidalChannel : public boundary
    {
    private:
        // --- Member Data ---

        // Geometric parameters defining the channel
        Foam::scalar length_;    // Channel length (along X)
        Foam::scalar bottomWidth_; // Width at the bottom (along Y at z=0)
        Foam::scalar topWidth_;   // Width at the top (along Y at z=height)
        Foam::scalar height_;    // Channel height at inlet (along Z at x=0)
        Foam::scalar slope_;     // Change in Z per unit length (negative slope means outlet is lower)

        // --- Private Methods ---

        /**
         * @brief Generates the points and faces based on stored parameters.
         * Called internally by initialize() or constructor.
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
         * Generates geometry with default parameters.
         */
        trapezoidalChannel();

        // Default destructor
        virtual ~trapezoidalChannel() = default;

        // --- Member Functions ---

        /**
         * @brief Initializes the channel parameters and calls develop() to generate geometry.
         * @param length Channel length (typically along X).
         * @param bottomWidth Channel bottom width (typically along Y at z=0).
         * @param topWidth Channel top width (typically along Y at z=height).
         * @param height Channel height (typically along Z at the inlet, x=0).
         * @param slope Change in Z per unit length (negative slope means outlet is lower).
         */
        void initialize(
            Foam::scalar length,
            Foam::scalar bottomWidth,
            Foam::scalar topWidth,
            Foam::scalar height,
            Foam::scalar slope
        );

        // --- Accessors ---

        /** @brief Get channel length. */
        Foam::scalar getLength() const { return length_; }
        /** @brief Get channel bottom width. */
        Foam::scalar getBottomWidth() const { return bottomWidth_; }
        /** @brief Get channel top width. */
        Foam::scalar getTopWidth() const { return topWidth_; }
        /** @brief Get channel height at inlet (x=0). */
        Foam::scalar getHeight() const { return height_; }
        /** @brief Get channel slope. */
        Foam::scalar getSlope() const { return slope_; }

        /**
         * @brief Get the list of patch type identifiers for each face.
         * @return Constant reference to the list of patch types.
         */
        const Foam::List<int>& getPatchTypes() const;
    };
}

#endif