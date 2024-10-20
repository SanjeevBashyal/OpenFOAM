#include "quickMesh.H"
#include "backgroundMesh.H"

using namespace Foam;
namespace Bashyal
{
    class backgroundBlock
    {
        public:
            backgroundMesh* ref_;  // Reference to another backgroundBlock (if needed)
            label identity_;        // Identifier for the block
            boundBox bounds_;       // Bounding box for the block
            pointField points_;    // Stores the vertices of the cube
            faceList faces_;       // Stores the six faces of the cube
            labelList owners_;     // Stores owner cells
            labelList neighbours_; // Stores neighbor cells
            
        private:
            /* data */
        public:
            backgroundBlock(backgroundMesh* ref, const label identity, const boundBox& bounds);
            void generateCubeGeometry();
            bool contains(const point &pt) const;
            void write(const std::string &meshDir) const;

            const pointField& getPoints() const { return points_; }
            const faceList& getFaces() const { return faces_; }
            const labelList& getOwners() const { return owners_; }
            const labelList& getNeighbours() const { return neighbours_; }

            ~backgroundBlock();
    };
}




