// --- START OF FILE particle.H ---

#ifndef particle_H
#define particle_H

#include "boundary.H"
#include "quaternion.H"

namespace Bashyal
{
    /**
     * @class particle
     * @brief Represents a physical particle for Discrete Element Method (DEM) simulations.
     *
     * This class extends the geometric 'boundary' with physical properties
     * such as mass, velocity, and forces. It provides methods to update its
     * state over time based on applied forces and torques, forming the core
     * of a DEM particle.
     */
    class particle : public boundary
    {
    private:
        //- Physical properties
        Foam::scalar mass_;                 // Mass of the particle
        Foam::tensor momentOfInertia_;      // Moment of inertia tensor in the body's local frame
        Foam::tensor invMomentOfInertia_;   // Inverse of the moment of inertia tensor

        //- Kinematic properties (State)
        Foam::point position_;              // Position of the center of mass in the world frame
        Foam::quaternion orientation_;      // Orientation of the particle (rotation from local to world frame)
        Foam::vector velocity_;             // Linear velocity of the center of mass
        Foam::vector angularVelocity_;      // Angular velocity in the world frame

        //- Dynamic properties (Forces and Torques)
        Foam::vector force_;                // Total force accumulated on the particle
        Foam::vector torque_;               // Total torque accumulated on the particle

    public:
        // * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

        /**
         * @brief Default constructor.
         * Initializes a particle with default physical properties (e.g., mass=1, zero velocity).
         */
        particle();

        /**
         * @brief Constructor from a boundary object and physical properties.
         * @param b The boundary object defining the particle's geometry.
         * @param mass The mass of the particle.
         * @param moi The moment of inertia tensor of the particle.
         * @param initialPosition The initial position of the particle's center of mass.
         */
        particle(
            const boundary& b,
            const Foam::scalar mass,
            const Foam::tensor& moi,
            const Foam::point& initialPosition
        );

        // * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

        /**
         * @brief Destructor.
         */
        ~particle();


        // * * * * * * * * * * * * * * * * Methods * * * * * * * * * * * * * * * //

        /**
         * @brief Updates the particle's state over a given time step 'dt'.
         * Uses a semi-implicit Euler integration scheme to update position and velocity.
         * @param dt The time step for the integration.
         */
        void update(Foam::scalar dt);

        /**
         * @brief Applies a force at a specific point on the particle.
         * Automatically calculates and adds the resulting torque.
         * @param force The force vector to apply.
         * @param applicationPoint The point in world coordinates where the force is applied.
         */
        void applyForce(const Foam::vector& force, const Foam::point& applicationPoint);

        /**
         * @brief Clears the accumulated forces and torques.
         * Should be called at the beginning of each time step before calculating new interactions.
         */
        void clearForceAndTorque();

        /**
         * @brief Returns the points of the geometry transformed into the world frame.
         * @return A pointField containing the vertices in their current world positions.
         */
        Foam::pointField getGlobalPoints() const;

        /**
         * @brief Writes the particle geometry to a VTP file using global coordinates.
         * @param filename The output VTP file name.
         */
        void writeVtp(const Foam::fileName& filename) const;


        // * * * * * * * * * * * * * * Accessors (Getters) * * * * * * * * * * * * * //

        Foam::scalar mass() const { return mass_; }
        const Foam::point& position() const { return position_; }
        const Foam::quaternion& orientation() const { return orientation_; }
        const Foam::vector& velocity() const { return velocity_; }
        const Foam::vector& angularVelocity() const { return angularVelocity_; }
        const Foam::vector& force() const { return force_; }
        const Foam::vector& torque() const { return torque_; }
        const Foam::tensor& momentOfInertia() const { return momentOfInertia_; }


        // * * * * * * * * * * * * * * Modifiers (Setters) * * * * * * * * * * * * * * //

        void setPosition(const Foam::point& p) { position_ = p; }
        void setVelocity(const Foam::vector& v) { velocity_ = v; }
        void setAngularVelocity(const Foam::vector& w) { angularVelocity_ = w; }
        void setOrientation(const Foam::quaternion& q) { orientation_ = q; }
    };
}

#endif

// --- END OF FILE particle.H ---