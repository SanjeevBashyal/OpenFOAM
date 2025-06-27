// --- START OF FILE particle.C ---

#include "particle.H"

namespace Bashyal
{
    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    particle::particle()
    :   boundary(), // Call base class default constructor
        mass_(1.0),
        momentOfInertia_(Foam::tensor::I),
        invMomentOfInertia_(Foam::tensor::I),
        position_(Foam::point::zero),
        orientation_(Foam::quaternion::I),
        velocity_(Foam::vector::zero),
        angularVelocity_(Foam::vector::zero),
        force_(Foam::vector::zero),
        torque_(Foam::vector::zero)
    {
    }

    particle::particle(
        const boundary& b,
        const Foam::scalar mass,
        const Foam::tensor& moi,
        const Foam::point& initialPosition
    )
    :   boundary(b), // Use boundary's copy constructor for geometry
        mass_(mass),
        momentOfInertia_(moi),
        position_(initialPosition),
        orientation_(Foam::quaternion::I), // Start with no rotation
        velocity_(Foam::vector::zero),
        angularVelocity_(Foam::vector::zero),
        force_(Foam::vector::zero),
        torque_(Foam::vector::zero)
    {
        // Pre-calculate the inverse of the moment of inertia
        if (mass_ <= 0)
        {
            FatalErrorInFunction
                << "Particle mass must be positive."
                << abort(Foam::FatalError);
        }
        invMomentOfInertia_ = inv(momentOfInertia_);
    }


    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    particle::~particle()
    {
    }


    // * * * * * * * * * * * * * * * * Methods  * * * * * * * * * * * * * * * //

    void particle::update(Foam::scalar dt)
    {
        // --- Translational Update (Semi-implicit Euler) ---
        // F = m*a  =>  a = F/m
        Foam::vector acceleration = force_ / mass_;
        velocity_ += acceleration * dt;
        position_ += velocity_ * dt;

        // --- Rotational Update ---
        // T = I*alpha => alpha = I^-1 * T
        // Transform inverse MOI tensor from local to world frame
        Foam::tensor R = Foam::tensor(orientation_); // Rotation matrix from orientation
        Foam::tensor invI_world = R * invMomentOfInertia_ * R.T();

        // Calculate angular acceleration
        Foam::vector angularAcceleration = invI_world & torque_;

        // Update angular velocity
        angularVelocity_ += angularAcceleration * dt;

        // Update orientation using the new angular velocity
        if (mag(angularVelocity_) > SMALL)
        {
            Foam::quaternion deltaQ(angularVelocity_ * dt); // Create quaternion from rotation vector
            orientation_ = deltaQ * orientation_;
            orientation_.normalise(); // Normalize to prevent numerical drift
        }
    }

    void particle::applyForce(const Foam::vector& force, const Foam::point& applicationPoint)
    {
        // Add to total force
        force_ += force;

        // Add to total torque: T = r x F
        // where r is the vector from the center of mass to the application point
        Foam::vector r = applicationPoint - position_;
        torque_ += (r ^ force);
    }

    void particle::clearForceAndTorque()
    {
        force_ = Foam::vector::zero;
        torque_ = Foam::vector::zero;
    }

    Foam::pointField particle::getGlobalPoints() const
    {
        // Get the local points from the base class
        const Foam::pointField& localPoints = this->points();
        Foam::pointField globalPoints(localPoints.size());

        // Create the transformation matrix from the orientation quaternion
        Foam::tensor R(orientation_);

        // Transform each local point to its global position
        forAll(localPoints, i)
        {
            // 1. Rotate the point from local to world orientation
            // 2. Translate the point by the particle's position
            globalPoints[i] = (R & localPoints[i]) + position_;
        }

        return globalPoints;
    }

} // End namespace Bashyal

// --- END OF FILE particle.C ---