// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitBox.hxx>
#include <BRepClass3d_ImplicitWork.hxx>
#include <BRep_TImplicitFace.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <typeinfo>

namespace
{
using Witness             = BRepClass3d_ImplicitBox::TopologyWitness;
using Triangle            = std::array<int, 3>;
using Edge                = std::pair<int, int>;
constexpr double Roundoff = std::numeric_limits<double>::epsilon();

bool exactIdentity(const TopLoc_Location& location)
{
  const auto matrix = location.Transformation();
  for (int row = 1; row <= 3; ++row)
    for (int column = 1; column <= 4; ++column)
      if (matrix.Value(row, column) != (row == column ? 1. : 0.))
        return false;
  return true;
}

struct Simplex
{
  int                   Count = 0;
  std::array<gp_Pnt, 4> Points;
  std::array<double, 4> Linear{};
};

// Tangential gradient of the nodal-linear function. Gram-Schmidt followed by
// a small triangular solve avoids a dimension-dependent projection contract.
bool direction(const Simplex& cell, gp_Vec& unit, double& norm, double& condition)
{
  std::array<gp_Vec, 3> basis;
  double                derivative[3]{};
  double                minimumAltitude = std::numeric_limits<double>::infinity();
  double                maximumEdge     = 0;
  for (int i = 0; i + 1 < cell.Count; ++i)
  {
    const gp_Vec edge(cell.Points[0], cell.Points[size_t(i + 1)]);
    maximumEdge = std::max(maximumEdge, edge.Magnitude());
    auto   q    = edge;
    double rhs  = cell.Linear[size_t(i + 1)] - cell.Linear[0];
    for (int j = 0; j < i; ++j)
    {
      const double coefficient = edge.Dot(basis[size_t(j)]);
      q -= basis[size_t(j)] * coefficient;
      rhs -= derivative[j] * coefficient;
    }
    const double altitude = q.Magnitude();
    if (!(altitude > 4096 * Roundoff * edge.Magnitude()) || !std::isfinite(altitude))
      return false;
    minimumAltitude  = std::min(minimumAltitude, altitude);
    basis[size_t(i)] = q / altitude;
    derivative[i]    = rhs / altitude;
  }
  gp_Vec gradient(0, 0, 0);
  for (int i = 0; i + 1 < cell.Count; ++i)
    gradient += basis[size_t(i)] * derivative[i];
  norm = gradient.Magnitude();
  // Allow for the triangular solve as well as anisotropic bases. This is an
  // engineering rounding reserve; ill-conditioned settings fail closed.
  condition = std::pow(1 + maximumEdge / minimumAltitude, 3) / minimumAltitude;
  if (!(norm > 0) || !std::isfinite(norm) || !std::isfinite(condition))
    return false;
  unit = gradient / norm;
  return true;
}

struct Checker
{
  const Geom_ImplicitGyroid& Field;
  Witness&                   Result;
  double                     L, H, Reserve, FieldReserve, GradientReserve;
  int                        MaximumDepth;
  size_t                     MaximumCells, Work = 0;

  bool check(const Simplex& cell, const gp_Vec& unit, double linearLower, int depth)
  {
    BRepClass3d_ImplicitWork::Checkpoint();
    const int dimension        = cell.Count - 1;
    Result.UnresolvedDimension = dimension;
    gp_Pnt low = cell.Points[0], high = low;
    gp_XYZ sum(0, 0, 0);
    double minLinear = cell.Linear[0], maxLinear = minLinear;
    for (int i = 0; i < cell.Count; ++i)
    {
      const auto& p = cell.Points[size_t(i)];
      sum += p.XYZ();
      minLinear = std::min(minLinear, cell.Linear[size_t(i)]);
      maxLinear = std::max(maxLinear, cell.Linear[size_t(i)]);
      for (int axis = 1; axis <= 3; ++axis)
      {
        low.SetCoord(axis, std::min(low.Coord(axis), p.Coord(axis)));
        high.SetCoord(axis, std::max(high.Coord(axis), p.Coord(axis)));
      }
    }
    Result.UnresolvedMinimum = low;
    Result.UnresolvedMaximum = high;
    if (Work >= MaximumCells)
    {
      Result.Reason = "homotopy_cell_budget_exhausted";
      return false;
    }
    ++Work;
    ++Result.Visited[size_t(dimension)];
    const gp_Pnt center(sum / cell.Count);
    double       radius = 0;
    for (int i = 0; i < cell.Count; ++i)
      radius = std::max(radius, center.Distance(cell.Points[size_t(i)]));
    radius += 8 * Reserve;
    const auto   g     = Field.Gradient(center);
    const double value = Field.Value(center);
    const double variation =
      std::min(L * radius, (g.Magnitude() + GradientReserve) * radius + .5 * H * radius * radius)
        * (1 + 128 * Roundoff)
      + FieldReserve;
    const double linearReserve =
      FieldReserve
      + 128 * Roundoff * (depth + 1) * (1 + std::max(std::abs(minLinear), std::abs(maxLinear)));
    // Same strict sign for F and L excludes every convex homotopy value.
    if ((value > variation && minLinear > linearReserve)
        || (value < -variation && maxLinear < -linearReserve))
    {
      ++Result.Excluded[size_t(dimension)];
      return true;
    }
    const double lower = std::min(linearLower, g.Dot(unit) - H * radius - GradientReserve);
    if (dimension > 0 && lower > GradientReserve && std::isfinite(lower))
    {
      ++Result.Transverse[size_t(dimension)];
      auto& minimum = Result.MinimumHomotopyDerivative;
      minimum       = minimum == 0 ? lower : std::min(minimum, lower);
      return true;
    }
    if (depth >= MaximumDepth || dimension == 0)
    {
      Result.Reason = "homotopy_transversality_unresolved";
      return false;
    }
    int    a = 0, b = 1;
    double longest = 0;
    for (int i = 0; i < cell.Count; ++i)
      for (int j = i + 1; j < cell.Count; ++j)
        if (cell.Points[size_t(i)].SquareDistance(cell.Points[size_t(j)]) > longest)
        {
          a       = i;
          b       = j;
          longest = cell.Points[size_t(i)].SquareDistance(cell.Points[size_t(j)]);
        }
    const gp_Pnt midpoint((cell.Points[size_t(a)].XYZ() + cell.Points[size_t(b)].XYZ()) * .5);
    if (midpoint.Distance(cell.Points[size_t(a)]) <= Reserve
        || midpoint.Distance(cell.Points[size_t(b)]) <= Reserve)
    {
      Result.Reason = "homotopy_subdivision_precision_exhausted";
      return false;
    }
    auto left = cell, right = cell;
    left.Points[size_t(a)] = right.Points[size_t(b)] = midpoint;
    left.Linear[size_t(a)]                           = right.Linear[size_t(b)] =
      (cell.Linear[size_t(a)] + cell.Linear[size_t(b)]) * .5;
    return check(left, unit, linearLower, depth + 1) && check(right, unit, linearLower, depth + 1);
  }

  bool simplex(const Simplex& cell)
  {
    gp_Vec unit(0, 0, 0);
    double norm = 0, condition = 0, magnitude = 1;
    for (int i = 0; i < cell.Count; ++i)
      magnitude = std::max(magnitude, std::abs(cell.Linear[size_t(i)]));
    const bool   valid = direction(cell, unit, norm, condition);
    const double lower = valid ? norm - 16384 * Roundoff * magnitude * condition : 0;
    return check(cell, unit, lower, 0);
  }
};

bool referenceTopology(Witness& result, std::vector<int>* triangleComponents = nullptr)
{
  struct Use
  {
    int Count = 0, Balance = 0, From = -1, To = -1;
  };

  std::map<Edge, Use>                       edges;
  std::vector<std::map<int, std::set<int>>> links(result.Nodes.size());
  std::vector<int>                          parent(result.Nodes.size());
  for (size_t i = 0; i < parent.size(); ++i)
    parent[i] = int(i);
  auto root = [&](int a) {
    while (parent[size_t(a)] != a)
    {
      parent[size_t(a)] = parent[size_t(parent[size_t(a)])];
      a                 = parent[size_t(a)];
    }
    return a;
  };
  for (const auto& triangle : result.Triangles)
    for (int i = 0; i < 3; ++i)
    {
      const int a = triangle[size_t(i)], b = triangle[size_t((i + 1) % 3)],
                c = triangle[size_t((i + 2) % 3)];
      if (a == b || b == c || c == a)
        return false;
      auto& use = edges[std::minmax(a, b)];
      ++use.Count;
      use.Balance += a < b ? 1 : -1;
      use.From                = a;
      use.To                  = b;
      parent[size_t(root(a))] = root(b);
      links[size_t(a)][b].insert(c);
      links[size_t(a)][c].insert(b);
    }
  std::map<int, int> next;
  std::set<int>      incoming;
  for (const auto& entry : edges)
  {
    const auto& e = entry.second;
    if (e.Count != 1 && !(e.Count == 2 && e.Balance == 0))
      return false;
    if (e.Count == 1 && (!next.emplace(e.From, e.To).second || !incoming.insert(e.To).second))
      return false;
  }
  if (next.size() != incoming.size())
    return false;
  // A connected path/cycle link rules out pinched vertices that paired mesh
  // edges alone would miss. Boundary vertices must have path links.
  for (size_t v = 0; v < links.size(); ++v)
  {
    const auto& link = links[v];
    if (link.empty())
      return false;
    int ends = 0;
    for (const auto& entry : link)
      if (entry.second.size() == 1)
        ++ends;
      else if (entry.second.size() != 2)
        return false;
    if (ends != (next.count(int(v)) ? 2 : 0))
      return false;
    std::set<int>    seen;
    std::vector<int> todo{link.begin()->first};
    while (!todo.empty())
    {
      const int n = todo.back();
      todo.pop_back();
      if (!seen.insert(n).second)
        continue;
      for (int neighbor : link.at(n))
        todo.push_back(neighbor);
    }
    if (seen.size() != link.size())
      return false;
  }
  while (!next.empty())
  {
    const int        start = next.begin()->first;
    int              node  = start;
    std::vector<int> loop{start};
    do
    {
      auto found = next.find(node);
      if (found == next.end())
        return false;
      node = found->second;
      next.erase(found);
      loop.push_back(node);
    } while (node != start);
    result.BoundaryLoops.push_back(std::move(loop));
  }
  std::map<int, int> euler;
  for (size_t v = 0; v < parent.size(); ++v)
    ++euler[root(int(v))];
  for (const auto& edge : edges)
    --euler[root(edge.first.first)];
  for (const auto& triangle : result.Triangles)
    ++euler[root(triangle[0])];
  for (const auto& component : euler)
    result.ComponentEuler.push_back(component.second);
  std::sort(result.ComponentEuler.begin(), result.ComponentEuler.end());
  if (triangleComponents)
  {
    std::vector<std::pair<int, int>> ordered;
    for (const auto& component : euler)
      ordered.emplace_back(component.second, component.first);
    std::sort(ordered.begin(), ordered.end());
    std::map<int, int> labels;
    for (size_t i = 0; i < ordered.size(); ++i)
      labels.emplace(ordered[i].second, int(i));
    for (const auto& triangle : result.Triangles)
      triangleComponents->push_back(labels.at(root(triangle[0])));
  }
  return true;
}
} // namespace

BRepClass3d_ImplicitBox::TopologyWitness BRepClass3d_ImplicitBox::CheckTopology(
  const TopoDS_Face& face,
  int                resolution,
  int                maximumDepth,
  size_t             maximumCells)
{
  return CheckTopology(face, resolution, maximumDepth, maximumCells, CornerPolicy::Reject);
}

BRepClass3d_ImplicitBox::TopologyWitness BRepClass3d_ImplicitBox::CheckTopology(
  const TopoDS_Face& face,
  int                resolution,
  int                maximumDepth,
  size_t             maximumCells,
  CornerPolicy       policy)
{
  BRepClass3d_ImplicitWork::EnterPhase(BRepClass3d_ImplicitWork::Phase::Topology);
  if (resolution < 4 || resolution > 64 || (resolution & (resolution - 1)) || maximumDepth < 0
      || maximumDepth > 16 || !maximumCells || maximumCells > 10000000)
    throw Standard_DomainError("Invalid native box topology settings");
  if (policy != CornerPolicy::Reject && policy != CornerPolicy::Regularized)
    throw Standard_DomainError("Unknown box topology corner policy");
  const bool regularized = policy == CornerPolicy::Regularized;
  // Odd subdivisions avoid systematic quarter/eighth-period interior samples;
  // only checked BOX corners may be exact zeros in this reference policy.
  if (regularized)
    --resolution;
  const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  if (data.IsNull() || typeid(*data) != typeid(BRep_TImplicitFace)
      || typeid(*data->Field()) != typeid(Geom_ImplicitGyroid)
      || (regularized ? !exactIdentity(face.Location() * data->Location())
                      : !face.Location().IsIdentity() || !data->Location().IsIdentity())
      || (face.Orientation() != TopAbs_FORWARD && face.Orientation() != TopAbs_REVERSED))
    throw Standard_NotImplemented("Topology witness requires an identity-placed native gyroid box");
  const auto&  field = *static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  const gp_Pnt low = data->Minimum(), high = data->Maximum();
  const double k = 2 * std::acos(-1.) / field.Period(), h = 4 * k * k, l = field.LipschitzBound();
  double       magnitude = 1;
  for (int axis = 1; axis <= 3; ++axis)
  {
    magnitude = std::max({magnitude, std::abs(low.Coord(axis)), std::abs(high.Coord(axis))});
    for (double value : {low.Coord(axis), high.Coord(axis)})
      if (!std::isfinite(k * value) || std::abs(k * value) >= 1.e6)
        throw Standard_DomainError("Topology trigonometric precision exhausted");
  }
  if (magnitude > 1.e50 || !std::isfinite(h) || h <= 0 || !std::isfinite(l))
    throw Standard_DomainError("Topology derivative/coordinate precision exhausted");
  const double reserve = 4096 * Roundoff * magnitude;
  for (int axis = 1; axis <= 3; ++axis)
    if (!((high.Coord(axis) - low.Coord(axis)) / resolution > 32 * reserve))
      throw Standard_DomainError("Topology grid coordinate precision exhausted");
  const double fieldReserve = 2048 * Roundoff * (1 + std::abs(field.IsoValue())) + l * reserve;
  Witness      result;
  result.Policy = policy;
  CornerModel cornerModel;
  if (regularized)
  {
    cornerModel = CheckCornerGerms(face, data->Tolerance());
    if (!cornerModel.Resolved)
    {
      result.Reason = "reference_corner_germs_unresolved";
      return result;
    }
    for (const auto& corner : cornerModel.Corners)
      if (corner.ExactZero && corner.Kind != CornerKind::SplitCone)
        result.PointContacts.push_back(corner.Anchor);
  }
  result.Resolution             = resolution;
  result.MaximumSimplexDiameter = low.Distance(high) / resolution;
  Checker             checker{field,
                              result,
                              l,
                              h,
                              reserve,
                              fieldReserve,
                              h * reserve + 2048 * Roundoff * l,
                              maximumDepth,
                              maximumCells};
  const int           side = resolution + 1;
  auto                id   = [side](int x, int y, int z) { return (z * side + y) * side + x; };
  std::vector<gp_Pnt> points(size_t(side) * side * side);
  std::vector<double> values(points.size());
  std::vector<int> cornerIndex(points.size(), -1);
  for (int z = 0; z <= resolution; ++z)
    for (int y = 0; y <= resolution; ++y)
      for (int x = 0; x <= resolution; ++x)
      {
        BRepClass3d_ImplicitWork::Checkpoint();
        const int index    = id(x, y, z);
        const int coord[3] = {x, y, z};
        gp_Pnt    p;
        for (int axis = 1; axis <= 3; ++axis)
          p.SetCoord(axis,
                     coord[axis - 1] == resolution
                       ? high.Coord(axis)
                       : low.Coord(axis)
                           + (high.Coord(axis) - low.Coord(axis)) * coord[axis - 1] / resolution);
        points[size_t(index)] = p;
        values[size_t(index)] = field.Value(p);
        bool fixedCorner = false;
        if (regularized && (x == 0 || x == resolution) && (y == 0 || y == resolution)
            && (z == 0 || z == resolution))
        {
          const int corner = int(x != 0) + 2 * int(y != 0) + 4 * int(z != 0);
          if (cornerModel.Corners[size_t(corner)].ExactZero)
          {
            values[size_t(index)] = 0; // Exact mathematical identity, not near-zero snapping.
            cornerIndex[size_t(index)] = corner;
            fixedCorner = true;
          }
        }
        if (checker.Work >= maximumCells
            || (!fixedCorner && !(std::abs(values[size_t(index)]) > fieldReserve)))
        {
          result.Reason = checker.Work >= maximumCells ? "homotopy_cell_budget_exhausted"
                                                       : "reference_vertex_sign_unresolved";
          result.UnresolvedDimension = 0;
          result.UnresolvedMinimum = result.UnresolvedMaximum = p;
          return result;
        }
        ++checker.Work;
        ++result.Visited[0];
        if (fixedCorner)
          ++result.FixedCornerCount;
        else
          ++result.Excluded[0];
      }
  std::map<Edge, int> roots;
  auto                root = [&](int a, int b) {
    if (values[size_t(a)] == 0)
      b = a;
    else if (values[size_t(b)] == 0)
      a = b;
    const Edge edge  = std::minmax(a, b);
    auto       found = roots.find(edge);
    if (found != roots.end())
      return found->second;
    const double t     = a == b ? 0 : values[size_t(a)] / (values[size_t(a)] - values[size_t(b)]);
    const int    index = int(result.Nodes.size());
    result.Nodes.emplace_back(points[size_t(a)].XYZ() * (1 - t) + points[size_t(b)].XYZ() * t);
    result.GridEdges.push_back({edge.first, edge.second});
    roots.emplace(edge, index);
    return index;
  };
  constexpr int offsets[8][3] =
    {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
  constexpr int tets[6][4] =
    {{0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7}, {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}};
  for (int z = 0; z < resolution; ++z)
    for (int y = 0; y < resolution; ++y)
      for (int x = 0; x < resolution; ++x)
      {
        BRepClass3d_ImplicitWork::Checkpoint();
        int nodes[8];
        for (int i = 0; i < 8; ++i)
          nodes[i] = id(x + offsets[i][0], y + offsets[i][1], z + offsets[i][2]);
        for (const auto& tet : tets)
        {
          if (regularized)
          {
            Simplex incident;
            incident.Count = 4;
            int fixed = -1;
            double valueMagnitude = 1;
            for (int i = 0; i < 4; ++i)
            {
              const int index = nodes[tet[i]];
              incident.Points[size_t(i)] = points[size_t(index)];
              incident.Linear[size_t(i)] = values[size_t(index)];
              valueMagnitude = std::max(valueMagnitude, std::abs(values[size_t(index)]));
              if (cornerIndex[size_t(index)] >= 0)
                fixed = cornerIndex[size_t(index)];
            }
            if (fixed >= 0)
            {
              gp_Vec unit;
              double norm = 0, condition = 0;
              if (!direction(incident, unit, norm, condition))
              {
                result.Reason = "reference_corner_linear_gradient_unresolved";
                return result;
              }
              const double error = 16384 * Roundoff * valueMagnitude * condition
                                   + checker.GradientReserve;
              const auto& germ = cornerModel.Corners[size_t(fixed)];
              for (int axis = 0; axis < 3; ++axis)
              {
                const double inward = unit.Coord(axis + 1) * norm
                                      * (fixed & (1 << axis) ? -1 : 1);
                if (!(germ.InwardDerivativeLower[size_t(axis)] > 0 ? inward > error
                                                                   : inward < -error))
                {
                  result.Reason = "reference_corner_homotopy_cone_unresolved";
                  return result;
                }
              }
            }
          }
          // Interior grid faces/edges are NOT boundary strata. Requiring their
          // restrictions to be transverse would reject harmless tangencies to
          // an artificial grid. Each full tetrahedron instead proves
          // grad(L).grad(F)>0 near possible homotopy zeros. Thus the continuous
          // field-gradient direction increases every incident linear piece,
          // including at internal simplex boundaries. On real box faces/edges
          // use the respective projected gradient and restricted L instead.
          for (int mask = 1; mask < 16; ++mask)
          {
            Simplex cell;
            for (int i = 0; i < 4; ++i)
              if (mask & (1 << i))
              {
                const int index                   = nodes[tet[i]];
                cell.Points[size_t(cell.Count)]   = points[size_t(index)];
                cell.Linear[size_t(cell.Count++)] = values[size_t(index)];
              }
            int fixedBoundaryAxes = 0;
            for (int axis = 1; axis <= 3; ++axis)
              for (double boundary : {low.Coord(axis), high.Coord(axis)})
              {
                bool fixed = true;
                for (int i = 0; i < cell.Count; ++i)
                  fixed = fixed && cell.Points[size_t(i)].Coord(axis) == boundary;
                fixedBoundaryAxes += fixed ? 1 : 0;
              }
            const bool boxStratum = cell.Count == 4 || (cell.Count == 3 && fixedBoundaryAxes == 1)
                                    || (cell.Count == 2 && fixedBoundaryAxes == 2);
            if (boxStratum && !checker.simplex(cell))
              return result;
          }
          Simplex cell;
          cell.Count = 4;
          int neg[4], pos[4], n = 0, p = 0, zero = -1;
          for (int i = 0; i < 4; ++i)
          {
            const int index        = nodes[tet[i]];
            cell.Points[size_t(i)] = points[size_t(index)];
            cell.Linear[size_t(i)] = values[size_t(index)];
            if (values[size_t(index)] == 0)
              zero = index;
            else if (values[size_t(index)] < 0)
              neg[n++] = index;
            else
              pos[p++] = index;
          }
          if (!n || !p)
            continue;
          if (result.Nodes.size() + 4 > 200000 || result.Triangles.size() + 2 > 400000)
          {
            result.Reason = "reference_mesh_budget_exhausted";
            return result;
          }
          gp_Vec unit;
          double norm, condition;
          if (!direction(cell, unit, norm, condition))
          {
            result.Reason = "reference_linear_gradient_unresolved";
            return result;
          }
          auto triangle = [&](int a, int b, int c) {
            const gp_Vec area =
              gp_Vec(result.Nodes[size_t(a)], result.Nodes[size_t(b)])
                .Crossed(gp_Vec(result.Nodes[size_t(a)], result.Nodes[size_t(c)]));
            if (area.Dot(unit) * (face.Orientation() == TopAbs_FORWARD ? 1 : -1) < 0)
              std::swap(b, c);
            result.Triangles.push_back({a, b, c});
          };
          if (zero >= 0)
          {
            // One exact box corner per tetrahedron. Strict signs are required
            // on both sides for a 2-D cut; isolated zeros create no fake face.
            const int a = root(zero, zero), b = root(neg[0], pos[0]);
            triangle(a, b, n == 2 ? root(neg[1], pos[0]) : root(neg[0], pos[1]));
          }
          else if (n == 1 || p == 1)
          {
            const int  a     = n == 1 ? neg[0] : pos[0];
            const int* other = n == 1 ? pos : neg;
            triangle(root(a, other[0]), root(a, other[1]), root(a, other[2]));
          }
          else
          {
            const int a = root(neg[0], pos[0]), b = root(neg[0], pos[1]), c = root(neg[1], pos[0]),
                      d = root(neg[1], pos[1]);
            triangle(a, b, c);
            triangle(b, d, c);
          }
        }
      }
  result.HomotopyResolved    = true;
  result.ReferenceManifold   = referenceTopology(result);
  result.Reason              = result.ReferenceManifold ? "homotopy_and_reference_manifold_resolved"
                                                        : "reference_manifold_unresolved";
  result.UnresolvedDimension = -1;
  return result;
}

BRepClass3d_ImplicitBox::RegionTopologyWitness BRepClass3d_ImplicitBox::CheckRegionTopology(
  const TopoDS_Face& face,
  int                resolution,
  int                maximumDepth,
  size_t             maximumCells)
{
  return CheckRegionTopology(face, resolution, maximumDepth, maximumCells, CornerPolicy::Reject);
}

BRepClass3d_ImplicitBox::RegionTopologyWitness BRepClass3d_ImplicitBox::CheckRegionTopology(
  const TopoDS_Face& face,
  int                resolution,
  int                maximumDepth,
  size_t             maximumCells,
  CornerPolicy       policy)
{
  RegionTopologyWitness result;
  result.Surface = CheckTopology(face, resolution, maximumDepth, maximumCells, policy);
  if (!result.Surface.HomotopyResolved || !result.Surface.ReferenceManifold)
  {
    result.Reason = result.Surface.Reason;
    return result;
  }
  resolution = result.Surface.Resolution;
  const bool regularized = policy == CornerPolicy::Regularized;
  const auto   data  = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  const auto&  field = *static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  const gp_Pnt low = data->Minimum(), high = data->Maximum();
  const int    side  = resolution + 1;
  auto         id    = [side](int x, int y, int z) { return (z * side + y) * side + x; };
  auto         point = [&](int index) {
    const int coordinates[3] = {index % side, (index / side) % side, index / (side * side)};
    gp_Pnt    p;
    for (int axis = 1; axis <= 3; ++axis)
      p.SetCoord(axis,
                 coordinates[axis - 1] == resolution
                   ? high.Coord(axis)
                   : low.Coord(axis)
                       + (high.Coord(axis) - low.Coord(axis)) * coordinates[axis - 1] / resolution);
    return p;
  };
  result.Nodes     = result.Surface.Nodes;
  result.Triangles = result.Surface.Triangles;
  result.TriangleSides.resize(result.Triangles.size(), -1);
  std::map<Edge, int> roots;
  for (size_t i = 0; i < result.Surface.GridEdges.size(); ++i)
  {
    const auto& edge = result.Surface.GridEdges[i];
    roots.emplace(Edge(edge[0], edge[1]), int(i));
  }
  std::map<int, int>    vertices;
  std::vector<double> values(size_t(side) * side * side,
                             std::numeric_limits<double>::quiet_NaN());
  auto                  value = [&](int index) {
    if (std::isfinite(values[size_t(index)]))
      return values[size_t(index)];
    const auto p = point(index);
    double f = field.Value(p);
    const int coordinates[3] = {index % side, (index / side) % side, index / (side * side)};
    if (regularized && (coordinates[0] == 0 || coordinates[0] == resolution)
        && (coordinates[1] == 0 || coordinates[1] == resolution)
        && (coordinates[2] == 0 || coordinates[2] == resolution))
    {
      double exact = 0;
      gp_Vec gradient;
      if (field.LatticeJet(p, exact, gradient) && exact == 0)
        f = 0; // Same exact box-corner identity as the surface witness.
    }
    values[size_t(index)] = f;
    return f;
  };
  auto vertex = [&](int index) {
    auto found = vertices.find(index);
    if (found != vertices.end())
      return found->second;
    const auto shared = roots.find({index, index});
    const int n = shared != roots.end() ? shared->second : int(result.Nodes.size());
    if (shared == roots.end())
      result.Nodes.push_back(point(index));
    vertices.emplace(index, n);
    return n;
  };
  const bool negative = face.Orientation() == TopAbs_FORWARD;
  // Both sides of every exterior grid square use its low-low / high-high
  // diagonal, exactly the restriction of the common Freudenthal complex.
  for (int axis = 0; axis < 3; ++axis)
    for (int highSide = 0; highSide < 2; ++highSide)
      for (int v = 0; v < resolution; ++v)
        for (int u = 0; u < resolution; ++u)
        {
          int       corners[4];
          const int du[4] = {0, 1, 1, 0}, dv[4] = {0, 0, 1, 1};
          for (int corner = 0; corner < 4; ++corner)
          {
            int coordinates[3];
            coordinates[axis]           = highSide ? resolution : 0;
            coordinates[(axis + 1) % 3] = u + du[corner];
            coordinates[(axis + 2) % 3] = v + dv[corner];
            corners[corner]             = id(coordinates[0], coordinates[1], coordinates[2]);
          }
          for (const auto& local : {std::array<int, 3>{0, 1, 2}, std::array<int, 3>{0, 2, 3}})
          {
            int triangle[3] = {corners[size_t(local[0])],
                               corners[size_t(local[1])],
                               corners[size_t(local[2])]};
            if (!highSide)
              std::swap(triangle[1], triangle[2]);
            bool material = false;
            for (int n : triangle)
              material |= value(n) != 0 && ((value(n) < 0) == negative);
            if (!material)
              continue; // A zero contact alone is not a retained 2-D cap.
            if (result.Nodes.size() + 3 > 200000 || result.Triangles.size() + 2 > 400000)
            {
              result.Reason = "closed_reference_mesh_budget_exhausted";
              return result;
            }
            std::vector<int> polygon;
            for (int edge = 0; edge < 3; ++edge)
            {
              const int  a = triangle[edge], b = triangle[(edge + 1) % 3];
              const double va = value(a), vb = value(b);
              const bool insideA = (va < 0) == negative, insideB = (vb < 0) == negative;
              if (va == 0 || insideA)
                polygon.push_back(vertex(a));
              if (va != 0 && vb != 0 && insideA != insideB)
              {
                const auto found = roots.find(std::minmax(a, b));
                if (found == roots.end())
                {
                  result.Reason = "closed_reference_missing_shared_root";
                  return result;
                }
                polygon.push_back(found->second);
              }
            }
            for (size_t i = 1; i + 1 < polygon.size(); ++i)
            {
              result.Triangles.push_back({polygon[0], polygon[i], polygon[i + 1]});
              result.TriangleSides.push_back(2 * axis + highSide);
            }
          }
        }
  Witness closed;
  closed.Nodes     = result.Nodes;
  closed.Triangles = result.Triangles;
  if (!referenceTopology(closed, &result.TriangleComponents) || !closed.BoundaryLoops.empty())
  {
    result.Reason = "closed_reference_manifold_unresolved";
    return result;
  }
  result.ComponentEuler = std::move(closed.ComponentEuler);
  // Translate before summing scalar triple products. These are diagnostics of
  // the linear reference, not volume/error bounds for the authoritative field.
  const gp_Pnt             origin((low.XYZ() + high.XYZ()) * .5);
  std::vector<long double> volumes(result.ComponentEuler.size(), 0);
  for (size_t i = 0; i < result.Triangles.size(); ++i)
  {
    const auto& t = result.Triangles[i];
    long double a[3], b[3], c[3];
    for (int axis = 1; axis <= 3; ++axis)
    {
      a[axis - 1] =
        static_cast<long double>(result.Nodes[size_t(t[0])].Coord(axis)) - origin.Coord(axis);
      b[axis - 1] =
        static_cast<long double>(result.Nodes[size_t(t[1])].Coord(axis)) - origin.Coord(axis);
      c[axis - 1] =
        static_cast<long double>(result.Nodes[size_t(t[2])].Coord(axis)) - origin.Coord(axis);
    }
    volumes[size_t(result.TriangleComponents[i])] +=
      (a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2])
       + a[2] * (b[0] * c[1] - b[1] * c[0]))
      / 6;
  }
  for (long double volume : volumes)
    result.ComponentSignedVolumes.push_back(double(volume));
  // Within a tetrahedron a retained affine half-space is convex. Every
  // nonempty retained VOLUME cell contains a strictly retained grid vertex.
  // Zero vertices must NOT join components through a mere contact. Joining
  // strict same-sign vertices along the common
  // tetrahedral 1-skeleton therefore gives VOLUME connectivity, not merely
  // connectivity of its possibly disconnected boundary shells.
  const size_t        gridSize = size_t(side) * side * side;
  std::vector<double> nodal(gridSize);
  std::vector<int>    parent(gridSize, -1);
  for (size_t i = 0; i < gridSize; ++i)
  {
    nodal[i] = value(int(i));
    if (nodal[i] != 0 && ((nodal[i] < 0) == negative))
      parent[i] = int(i);
  }
  auto representative = [&](int n) {
    while (parent[size_t(n)] != n)
    {
      parent[size_t(n)] = parent[size_t(parent[size_t(n)])];
      n                 = parent[size_t(n)];
    }
    return n;
  };
  constexpr int offsets[8][3] =
    {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
  constexpr int tets[6][4] =
    {{0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7}, {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}};
  std::map<int, std::vector<int>> zeroNeighbors;
  for (int z = 0; z < resolution; ++z)
    for (int y = 0; y < resolution; ++y)
      for (int x = 0; x < resolution; ++x)
      {
        int cube[8];
        for (int i = 0; i < 8; ++i)
          cube[i] = id(x + offsets[i][0], y + offsets[i][1], z + offsets[i][2]);
        for (const auto& tet : tets)
        {
          int anchor = -1;
          for (int local : tet)
          {
            const int n = cube[local];
            if (parent[size_t(n)] < 0)
              continue;
            if (anchor < 0)
              anchor = n;
            else
              parent[size_t(representative(n))] = representative(anchor);
          }
          if (anchor >= 0)
            for (int local : tet)
              if (nodal[size_t(cube[local])] == 0)
                zeroNeighbors[cube[local]].push_back(anchor);
        }
      }
  std::map<int, int> regions;
  for (size_t i = 0; i < gridSize; ++i)
    if (parent[i] >= 0)
      regions.emplace(representative(int(i)), 0);
  for (auto& entry : regions)
    entry.second = result.RegionCount++;
  std::map<int, int> zeroRegions;
  for (const auto& entry : zeroNeighbors)
  {
    int region = -1;
    for (int neighbor : entry.second)
    {
      const int adjacent = regions.at(representative(neighbor));
      if (region >= 0 && region != adjacent)
      {
        result.Reason = "closed_reference_corner_volume_link_unresolved";
        return result;
      }
      region = adjacent;
    }
    zeroRegions.emplace(entry.first, region);
  }
  std::vector<int> nodeRegion(result.Nodes.size(), -1);
  for (size_t i = 0; i < result.Surface.GridEdges.size(); ++i)
  {
    const auto& edge     = result.Surface.GridEdges[i];
    if (edge[0] == edge[1])
    {
      const auto found = zeroRegions.find(edge[0]);
      if (found == zeroRegions.end())
      {
        result.Reason = "closed_reference_corner_region_unresolved";
        return result;
      }
      nodeRegion[i] = found->second;
      continue;
    }
    const int   retained = parent[size_t(edge[0])] >= 0 ? edge[0] : edge[1];
    if (parent[size_t(retained)] < 0)
    {
      result.Reason = "closed_reference_root_region_unresolved";
      return result;
    }
    nodeRegion[i] = regions.at(representative(retained));
  }
  for (const auto& entry : vertices)
  {
    if (nodal[size_t(entry.first)] == 0)
    {
      const auto found = zeroRegions.find(entry.first);
      if (found == zeroRegions.end())
      {
        result.Reason = "closed_reference_cap_corner_region_unresolved";
        return result;
      }
      nodeRegion[size_t(entry.second)] = found->second;
    }
    else
      nodeRegion[size_t(entry.second)] = regions.at(representative(entry.first));
  }
  result.ComponentRegions.assign(result.ComponentEuler.size(), -1);
  for (size_t i = 0; i < result.Triangles.size(); ++i)
  {
    const auto& triangle        = result.Triangles[i];
    const int   region          = nodeRegion[size_t(triangle[0])];
    auto&       componentRegion = result.ComponentRegions[size_t(result.TriangleComponents[i])];
    if (region < 0 || nodeRegion[size_t(triangle[1])] != region
        || nodeRegion[size_t(triangle[2])] != region
        || (componentRegion >= 0 && componentRegion != region))
    {
      result.Reason = "closed_reference_shell_region_grouping_unresolved";
      return result;
    }
    componentRegion = region;
  }
  std::vector<long double> regionVolumes(size_t(result.RegionCount), 0);
  std::vector<bool>        regionUsed(size_t(result.RegionCount), false);
  for (size_t component = 0; component < volumes.size(); ++component)
  {
    const int region = result.ComponentRegions[component];
    if (region < 0)
    {
      result.Reason = "closed_reference_unassigned_shell";
      return result;
    }
    regionVolumes[size_t(region)] += volumes[component];
    regionUsed[size_t(region)] = true;
  }
  for (bool used : regionUsed)
    if (!used)
    {
      result.Reason = "closed_reference_region_without_boundary";
      return result;
    }
  for (long double volume : regionVolumes)
    result.RegionVolumes.push_back(double(volume));
  result.ClosedReferenceManifold = true;
  result.Reason                  = "closed_region_reference_resolved";
  return result;
}
