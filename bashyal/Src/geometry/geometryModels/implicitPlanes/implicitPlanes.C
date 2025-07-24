#include "halfEdgeMesh.H"
#include "implicitPlanes.H"
#include "indexedFaceSet.H"
#include "foamCGALConverter.H"
#include <vector>

namespace Foam
{
  namespace particleModels
  {
    implicitPlanes::implicitPlanes() = default;
    implicitPlanes::implicitPlanes(const List<Plane> &planes, const point &centroid)
        : planes_(planes), centroid_(centroid) {}

    bool implicitPlanes::isInside(const point &p) const
    {
      for (const auto &plane : planes_)
      {
        if ((plane.normal & (p - centroid_)) - plane.distance > SMALL)
        {
          return false;
        }
      }
      return true;
    }

    indexedFaceSet implicitPlanes::toIndexedFaceSet() const
    {
      std::vector<CgalPlane_3> cgalPlanes;
      for (const auto &plane : planes_)
      {
        const vector &n = plane.normal;
        scalar plane_d = -((n & centroid_) + plane.distance);
        cgalPlanes.emplace_back(n.x(), n.y(), n.z(), plane_d);
      }
      CgalPolyhedron cgalPoly;
      CGAL::halfspace_intersection_3(cgalPlanes.begin(), cgalPlanes.end(),
                                     cgalPoly);
      HalfEdgeMesh heMesh(cgalPoly);
      return heMesh.toIndexedFaceSet();
    }

    void implicitPlanes::print(Ostream &os) const
    {
      os << "--- implicitPlanes ---" << endl
         << "Planes: " << planes_.size() << nl << "Centroid: " << centroid_ << endl;
    }
  } // namespace particleModels
} // namespace Foam