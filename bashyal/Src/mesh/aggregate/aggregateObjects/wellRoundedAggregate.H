#ifndef wellRoundedAggregate_H
#define wellRoundedAggregate_H

#include "quickInclude.H"
#include "aggregate.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"

namespace Bashyal
{
    class wellRoundedAggregate : public aggregate
    {
    public:
        // Member variables
        float s1_;              // Minimum size
        float s2_;              // Maximum size
        float roundnessFactor_; // Controls deviation from a sphere (0 to 1)
        int resolution_;        // Number of divisions for point generation
        Foam::label identifier_; // Unique identifier

    public:
        // Constructor
        wellRoundedAggregate(float s1, float s2, float roundnessFactor, int resolution, Foam::label identifier);

        // Destructor
        ~wellRoundedAggregate() = default;

        // Member functions
        void createFaces();                       // Generate faces
        Foam::List<Foam::face> createQuadFaces(); // Generate quadrilateral faces
        void createTriSurface();                  // Create triangulated surface

    private:
        // Helper functions
        Foam::pointField generatePoints(Foam::scalar s, Foam::scalar sx, Foam::scalar sy, Foam::scalar sz, int resolution);
        Foam::pointField createSpherePoints(Foam::scalar radius, int resolution);
    };
}

#endif
