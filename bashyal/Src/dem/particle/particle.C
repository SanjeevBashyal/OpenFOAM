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
        Foam::tensor R = orientation_.R(); // Get rotation matrix from quaternion
        Foam::tensor invI_world = R & invMomentOfInertia_ & R.T();

        // Calculate angular acceleration
        Foam::vector angularAcceleration = invI_world & torque_;

        // Update angular velocity
        angularVelocity_ += angularAcceleration * dt;

        Foam::scalar omegaMag = Foam::mag(angularVelocity_);

        // Update orientation using the new angular velocity
        if (omegaMag > Foam::SMALL)
        {
            Foam::vector rotationAxis = angularVelocity_ / omegaMag; // Normalize to get the axis
            Foam::scalar rotationAngle = omegaMag * dt;         // angle = speed * time
            Foam::quaternion deltaQ(rotationAxis, rotationAngle); // Create quaternion from rotation vector
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
        const Foam::pointField& localPoints = this->vertices();
        Foam::pointField globalPoints(localPoints.size());

        // Create the transformation matrix from the orientation quaternion
        Foam::tensor R = orientation_.R();

        // Transform each local point to its global position
        forAll(localPoints, i)
        {
            // 1. Rotate the point from local to world orientation
            // 2. Translate the point by the particle's position
            globalPoints[i] = (R & localPoints[i]) + position_;
        }

        return globalPoints;
    }

    void particle::writeVtp(const Foam::fileName& filename) const
    {
        Foam::OFstream vtpFile(filename);
        if (!vtpFile.good())
        {
            FatalErrorIn("writeVtp") << "Cannot open file " << filename << Foam::exit(Foam::FatalError);
        }

        // Get the global points (transformed to world coordinates)
        const Foam::pointField& globalVertices = getGlobalPoints();
        const Foam::faceList& faces = this->faces();

        // VTK XML Header
        vtpFile << "<?xml version=\"1.0\"?>" << Foam::endl;
        vtpFile << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">" << Foam::endl;
        vtpFile << "  <PolyData>" << Foam::endl;

        // Piece defines the geometry. We have one object, so one piece.
        vtpFile << "    <Piece NumberOfPoints=\"" << globalVertices.size()
                << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\""
                << faces.size() << "\">" << Foam::endl;

        // 1. Write Points (Vertices) - using global coordinates
        vtpFile << "      <Points>" << Foam::endl;
        vtpFile << "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << Foam::endl;
        for (const Foam::point& pt : globalVertices)
        {
            vtpFile << "          " << pt.x() << " " << pt.y() << " " << pt.z() << Foam::endl;
        }
        vtpFile << "        </DataArray>" << Foam::endl;
        vtpFile << "      </Points>" << Foam::endl;

        // 2. Write Polygons (Faces)
        vtpFile << "      <Polys>" << Foam::endl;
        // a) Connectivity: a flat list of all vertex indices for all faces
        vtpFile << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << Foam::endl;
        for (const Foam::face& f : faces)
        {
            vtpFile << "          ";
            for (const Foam::label vIdx : f)
            {
                vtpFile << vIdx << " ";
            }
            vtpFile << Foam::endl;
        }
        vtpFile << "        </DataArray>" << Foam::endl;

        // b) Offsets: the cumulative count of vertices per face
        vtpFile << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << Foam::endl;
        vtpFile << "          ";
        Foam::label offset = 0;
        for (const Foam::face& f : faces)
        {
            offset += f.size();
            vtpFile << offset << " ";
        }
        vtpFile << Foam::endl;
        vtpFile << "        </DataArray>" << Foam::endl;
        vtpFile << "      </Polys>" << Foam::endl;

        // VTK XML Footer
        vtpFile << "    </Piece>" << Foam::endl;
        vtpFile << "  </PolyData>" << Foam::endl;
        vtpFile << "</VTKFile>" << Foam::endl;
    }

} // End namespace Bashyal

// --- END OF FILE particle.C ---