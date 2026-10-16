// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitBox.hxx>
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

bool referenceTopology(Witness& result)
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
  return true;
}
} // namespace

BRepClass3d_ImplicitBox::TopologyWitness BRepClass3d_ImplicitBox::CheckTopology(
  const TopoDS_Face& face,
  int                resolution,
  int                maximumDepth,
  size_t             maximumCells)
{
  if (resolution < 4 || resolution > 64 || (resolution & (resolution - 1)) || maximumDepth < 0
      || maximumDepth > 16 || !maximumCells || maximumCells > 10000000)
    throw Standard_DomainError("Invalid native box topology settings");
  const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  if (data.IsNull() || typeid(*data) != typeid(BRep_TImplicitFace)
      || typeid(*data->Field()) != typeid(Geom_ImplicitGyroid) || !face.Location().IsIdentity()
      || !data->Location().IsIdentity()
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
  for (int z = 0; z <= resolution; ++z)
    for (int y = 0; y <= resolution; ++y)
      for (int x = 0; x <= resolution; ++x)
      {
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
        if (checker.Work >= maximumCells || !(std::abs(values[size_t(index)]) > fieldReserve))
        {
          result.Reason = checker.Work >= maximumCells ? "homotopy_cell_budget_exhausted"
                                                       : "reference_vertex_sign_unresolved";
          result.UnresolvedDimension = 0;
          result.UnresolvedMinimum = result.UnresolvedMaximum = p;
          return result;
        }
        ++checker.Work;
        ++result.Visited[0];
        ++result.Excluded[0];
      }
  std::map<Edge, int> roots;
  auto                root = [&](int a, int b) {
    const Edge edge  = std::minmax(a, b);
    auto       found = roots.find(edge);
    if (found != roots.end())
      return found->second;
    const double t     = values[size_t(a)] / (values[size_t(a)] - values[size_t(b)]);
    const int    index = int(result.Nodes.size());
    result.Nodes.emplace_back(points[size_t(a)].XYZ() * (1 - t) + points[size_t(b)].XYZ() * t);
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
        int nodes[8];
        for (int i = 0; i < 8; ++i)
          nodes[i] = id(x + offsets[i][0], y + offsets[i][1], z + offsets[i][2]);
        for (const auto& tet : tets)
        {
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
          int neg[4], pos[4], n = 0, p = 0;
          for (int i = 0; i < 4; ++i)
          {
            const int index        = nodes[tet[i]];
            cell.Points[size_t(i)] = points[size_t(index)];
            cell.Linear[size_t(i)] = values[size_t(index)];
            if (values[size_t(index)] < 0)
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
          if (n == 1 || p == 1)
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
