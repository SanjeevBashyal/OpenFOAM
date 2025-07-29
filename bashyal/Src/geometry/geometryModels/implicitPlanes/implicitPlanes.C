#include "halfEdgeMesh.H"
#include "implicitPlanes.H"
#include "indexedFaceSet.H"
#include "foamCGALConverter.H"
#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using CGALPlane_3 = Kernel::Plane_3;
using CGALPolyhedron = CGAL::Polyhedron_3<Kernel>;

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
      std::vector<CGALPlane_3> cgalPlanes;
      for (const auto &plane : planes_)
      {
        const vector &n = plane.normal;
        scalar plane_d = -((n & centroid_) + plane.distance);
        cgalPlanes.emplace_back(n.x(), n.y(), n.z(), plane_d);
      }
      CGALPolyhedron cgalPoly;
      CGAL::halfspace_intersection_3(cgalPlanes.begin(), cgalPlanes.end(), cgalPoly);
      Foam::particleModels::halfEdgeMesh heMesh(cgalPoly);
      return heMesh.toIndexedFaceSet();
    }

    void implicitPlanes::print(Ostream &os) const
    {
      os << "--- implicitPlanes ---" << endl
         << "Planes: " << planes_.size() << nl << "Centroid: " << centroid_ << endl;
    }
  } // namespace particleModels
} // namespace Foam
