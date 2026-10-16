// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepMesh_ImplicitBoundary.hxx>
#include <BRep_TImplicitFace.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <unordered_map>

namespace
{
using Result   = BRepMesh_ImplicitBoundary::Result;
using Triangle = std::array<int, 3>;
using Key      = std::pair<int, int>;

void require(bool test, const char* message)
{
  if (!test)
    throw Standard_NotImplemented(message);
}

Key key(int a, int b)
{
  return std::minmax(a, b);
}

struct Use
{
  int Count = 0, Balance = 0, From = -1, To = -1;
};

void add(std::map<Key, Use>& edges, int a, int b)
{
  require(a != b, "Degenerate candidate mesh edge");
  auto& e = edges[key(a, b)];
  ++e.Count;
  e.Balance += a < b ? 1 : -1;
  e.From = a;
  e.To   = b;
}

std::vector<std::vector<int>> boundaryLoops(const std::vector<Triangle>& triangles)
{
  std::map<Key, Use> edges;
  for (const auto& t : triangles)
    for (int i = 0; i < 3; ++i)
      add(edges, t[i], t[(i + 1) % 3]);
  std::map<int, int> next;
  std::set<int>      incoming;
  for (const auto& entry : edges)
  {
    const auto& e = entry.second;
    require(e.Count == 1 || (e.Count == 2 && e.Balance == 0),
            "Nonmanifold interior candidate mesh");
    if (e.Count == 1)
    {
      require(next.emplace(e.From, e.To).second && incoming.insert(e.To).second,
              "Ambiguous candidate mesh boundary vertex");
    }
  }
  require(next.size() == incoming.size(), "Unbalanced candidate mesh boundary");
  std::vector<std::vector<int>> loops;
  while (!next.empty())
  {
    const int        start = next.begin()->first;
    std::vector<int> loop{start};
    int              node = start;
    do
    {
      auto found = next.find(node);
      require(found != next.end(), "Open candidate mesh boundary");
      node = found->second;
      next.erase(found);
      loop.push_back(node);
    } while (node != start);
    require(loop.size() >= 4, "Degenerate candidate mesh loop");
    loops.push_back(std::move(loop));
  }
  return loops;
}

// Triangle-wide interpolation remainder and a monotone normal line witness.
// This bounds distance to the untrimmed field, NOT coverage or clipped roots.
double bound(const Result&              mesh,
             const Triangle&            t,
             const Geom_ImplicitGyroid& field,
             double                     hessian,
             double                     tolerance,
             double                     reserve,
             double                     orientation)
{
  const auto&  a = mesh.Nodes[size_t(t[0])];
  const auto&  b = mesh.Nodes[size_t(t[1])];
  const auto&  c = mesh.Nodes[size_t(t[2])];
  const gp_Pnt center((a.XYZ() + b.XYZ() + c.XYZ()) / 3);
  const gp_Vec area     = gp_Vec(a, b).Crossed(gp_Vec(a, c));
  const auto   gradient = field.Gradient(center);
  if (!(area.SquareMagnitude() > reserve * reserve * reserve * reserve)
      || !(orientation * area.Dot(gradient) > 0))
    return std::numeric_limits<double>::infinity();
  const double diameter2 =
    std::max({a.SquareDistance(b), b.SquareDistance(c), c.SquareDistance(a)});
  const double radius = std::sqrt(
    std::max({center.SquareDistance(a), center.SquareDistance(b), center.SquareDistance(c)}));
  const double residual =
    std::max({std::abs(field.Value(a)), std::abs(field.Value(b)), std::abs(field.Value(c))});
  const double fieldReserve =
    field.LipschitzBound() * reserve + 1024 * std::numeric_limits<double>::epsilon();
  const double lower = gradient.Magnitude() - hessian * radius - hessian * reserve;
  if (!(lower > 0))
    return std::numeric_limits<double>::infinity();
  const double error        = residual + hessian * diameter2 / 6 + fieldReserve;
  const double discriminant = lower * lower - 2 * hessian * error;
  if (!(discriminant > 0))
    return std::numeric_limits<double>::infinity();
  const double distance = 2 * error / (lower + std::sqrt(discriminant)) + reserve;
  return distance <= tolerance ? distance : std::numeric_limits<double>::infinity();
}

// Symmetric sample-to-polyline distance is used ONLY for loop correspondence;
// no nodes are welded. All candidate pairings must satisfy a separation gate.
double pointPath(const gp_Pnt& p, const std::vector<int>& path, const Result& mesh)
{
  double best = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i < path.size(); ++i)
  {
    const auto&  a = mesh.Nodes[size_t(path[i - 1])];
    const auto&  b = mesh.Nodes[size_t(path[i])];
    const gp_Vec ab(a, b), ap(a, p);
    const double u = std::clamp(ap.Dot(ab) / ab.SquareMagnitude(), 0., 1.);
    best           = std::min(best, p.SquareDistance(a.Translated(ab * u)));
  }
  return std::sqrt(best);
}

double pathDistance(const std::vector<int>& a, const std::vector<int>& b, const Result& mesh)
{
  double result = 0;
  for (int node : a)
    result = std::max(result, pointPath(mesh.Nodes[size_t(node)], b, mesh));
  for (int node : b)
    result = std::max(result, pointPath(mesh.Nodes[size_t(node)], a, mesh));
  return result;
}

void collar(Result&                    mesh,
            const std::vector<int>&    native,
            const std::vector<int>&    inner,
            const Geom_ImplicitGyroid& field,
            double                     hessian,
            double                     tolerance,
            double                     reserve,
            double                     orientation,
            int                        maximumTriangles)
{
  const size_t n = native.size() - 1, m = inner.size() - 1;
  require(n * m <= 8000000, "Native mesh collar correspondence budget exhausted");
  size_t startN = 0, startM = 0;
  double closest = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
    {
      const double distance =
        mesh.Nodes[size_t(native[i])].SquareDistance(mesh.Nodes[size_t(inner[j])]);
      if (distance < closest)
      {
        closest = distance;
        startN  = i;
        startM  = j;
      }
    }
  const auto                 a     = [&](size_t i) { return native[(i + startN) % n]; };
  const auto                 b     = [&](size_t j) { return inner[(j + startM) % m]; };
  const size_t               width = m + 1;
  std::vector<double>        costs((n + 1) * width, std::numeric_limits<double>::infinity());
  std::vector<unsigned char> moves(costs.size(), 0);
  costs[0]           = 0;
  const auto quality = [&](const Triangle& triangle) {
    const auto&  p = mesh.Nodes[size_t(triangle[0])];
    const auto&  q = mesh.Nodes[size_t(triangle[1])];
    const auto&  r = mesh.Nodes[size_t(triangle[2])];
    const gp_Pnt center((p.XYZ() + q.XYZ() + r.XYZ()) / 3);
    if (!(orientation * gp_Vec(p, q).Crossed(gp_Vec(p, r)).Dot(field.Gradient(center)) > 0))
      return std::numeric_limits<double>::infinity();
    const double distance = bound(mesh, triangle, field, hessian, tolerance, reserve, orientation);
    const double penalty  = std::isfinite(distance) ? 1 + distance / tolerance : 10;
    return penalty * (p.SquareDistance(q) + q.SquareDistance(r) + r.SquareDistance(p));
  };
  for (size_t i = 0; i <= n; ++i)
    for (size_t j = 0; j <= m; ++j)
    {
      const size_t index = i * width + j;
      if (i > 0)
      {
        const double value = costs[index - width] + quality({a(i - 1), a(i), b(j)});
        if (value < costs[index])
        {
          costs[index] = value;
          moves[index] = 1;
        }
      }
      if (j > 0)
      {
        const double value = costs[index - 1] + quality({a(i), b(j), b(j - 1)});
        if (value < costs[index])
        {
          costs[index] = value;
          moves[index] = 2;
        }
      }
    }
  require(std::isfinite(costs.back()), "No oriented native mesh collar");
  require(mesh.ImplicitTriangles.size() + n + m
            <= size_t(maximumTriangles) - mesh.CapTriangles.size(),
          "Native mesh triangle budget exhausted");
  size_t i = n, j = m;
  while (i || j)
  {
    const auto move = moves[i * width + j];
    require(move != 0, "Incomplete native mesh collar");
    if (move == 1)
    {
      mesh.ImplicitTriangles.push_back({a(i - 1), a(i), b(j)});
      --i;
    }
    else
    {
      mesh.ImplicitTriangles.push_back({a(i), b(j), b(j - 1)});
      --j;
    }
  }
}

void refine(Result&                    mesh,
            const Geom_ImplicitGyroid& field,
            const gp_Pnt&              minimum,
            const gp_Pnt&              maximum,
            double                     hessian,
            double                     tolerance,
            double                     reserve,
            double                     orientation,
            int                        maximumNodes,
            int                        maximumTriangles)
{
  std::set<Key> fixed;
  for (const auto& wire : mesh.Wires)
    for (size_t i = 1; i < wire.size(); ++i)
      fixed.insert(key(wire[i - 1], wire[i]));
  for (int level = 0; level <= 8; ++level)
  {
    std::map<Key, int> split;
    for (const auto& t : mesh.ImplicitTriangles)
      if (!(bound(mesh, t, field, hessian, tolerance, reserve, orientation) <= tolerance))
        for (int i = 0; i < 3; ++i)
        {
          const auto e = key(t[i], t[(i + 1) % 3]);
          if (!fixed.count(e))
            split.emplace(e, -1);
        }
    if (split.empty())
      return; // final independent bound check still applies
    require(level < 8, "Native mesh projected refinement limit exhausted");
    require(mesh.Nodes.size() + split.size() <= size_t(maximumNodes),
            "Native mesh refinement node budget exhausted");
    for (auto& entry : split)
    {
      const gp_Pnt origin(
        (mesh.Nodes[size_t(entry.first.first)].XYZ() + mesh.Nodes[size_t(entry.first.second)].XYZ())
        * .5);
      gp_Pnt p = origin;
      for (int iteration = 0; iteration < 24; ++iteration)
      {
        const double value    = field.Value(p);
        const auto   gradient = field.Gradient(p);
        require(gradient.Magnitude() > field.LipschitzBound() * 1.e-12,
                "Critical native mesh refinement projection");
        if (std::abs(value) <= gradient.Magnitude() * tolerance / 8192)
          break;
        const auto step     = gradient * (-value / gradient.SquareMagnitude());
        bool       accepted = false;
        double     factor   = 1;
        for (int backtrack = 0; backtrack < 24; ++backtrack, factor *= .5)
        {
          const gp_Pnt q      = p.Translated(step * factor);
          bool         inside = true;
          for (int axis = 1; axis <= 3; ++axis)
            inside &= q.Coord(axis) >= minimum.Coord(axis) && q.Coord(axis) <= maximum.Coord(axis);
          if (inside && q.Distance(origin) <= 2 * tolerance
              && std::abs(field.Value(q)) < std::abs(value))
          {
            p        = q;
            accepted = true;
            break;
          }
        }
        require(accepted, "Native mesh refinement projection unresolved");
      }
      entry.second = int(mesh.Nodes.size());
      mesh.Nodes.push_back(p);
    }
    std::vector<Triangle> triangles;
    const auto            append = [&](int a, int b, int c) {
      require(triangles.size() + mesh.CapTriangles.size() < size_t(maximumTriangles),
              "Native mesh refinement triangle budget exhausted");
      triangles.push_back({a, b, c});
    };
    for (const auto& t : mesh.ImplicitTriangles)
    {
      int middle[3], count = 0;
      for (int i = 0; i < 3; ++i)
      {
        const auto found = split.find(key(t[i], t[(i + 1) % 3]));
        middle[i]        = found == split.end() ? -1 : found->second;
        count += middle[i] >= 0;
      }
      if (count == 0)
        append(t[0], t[1], t[2]);
      else if (count == 3)
      {
        append(t[0], middle[0], middle[2]);
        append(middle[0], t[1], middle[1]);
        append(middle[2], middle[1], t[2]);
        append(middle[0], middle[1], middle[2]);
      }
      else
      {
        int i = 0;
        while (!(middle[i] >= 0 && (count == 1 || middle[(i + 1) % 3] >= 0)))
          ++i;
        const int a = t[i], b = t[(i + 1) % 3], c = t[(i + 2) % 3], ab = middle[i];
        if (count == 1)
        {
          append(a, ab, c);
          append(ab, b, c);
        }
        else
        {
          const int bc = middle[(i + 1) % 3];
          append(b, bc, ab);
          append(a, ab, c);
          append(ab, bc, c);
        }
      }
    }
    mesh.ImplicitTriangles = std::move(triangles);
  }
}

Result extract(const Result&              prepared,
               const BRep_TImplicitFace&  face,
               const Geom_ImplicitGyroid& field,
               double                     orientation,
               double                     tolerance,
               int                        cells,
               int                        maximumNodes,
               int                        maximumTriangles)
{
  Result       result  = prepared;
  const auto   minimum = face.Minimum(), maximum = face.Maximum();
  const gp_XYZ extent  = maximum.XYZ() - minimum.XYZ();
  const double hessian = 4 * std::pow(2 * std::acos(-1.) / field.Period(), 2);
  require(std::isfinite(hessian) && hessian > 0 && std::isfinite(field.LipschitzBound()),
          "Native surface mesh derivative bound exhausted");
  double       magnitude = std::max({1., minimum.XYZ().Modulus(), maximum.XYZ().Modulus()});
  const double reserve   = 16384 * std::numeric_limits<double>::epsilon() * magnitude;
  require(tolerance > 4096 * reserve, "Native surface mesh coordinate precision exhausted");
  const double longest = std::max({extent.X(), extent.Y(), extent.Z()});
  const int    nx      = std::max(1, int(std::ceil(cells * (extent.X() / longest))));
  const int    ny      = std::max(1, int(std::ceil(cells * (extent.Y() / longest))));
  const int    nz      = std::max(1, int(std::ceil(cells * (extent.Z() / longest))));
  const int    sx = nx + 1, sy = ny + 1;
  const size_t gridSize = size_t(sx) * sy * (nz + 1);
  require(gridSize <= 2200000, "Native mesh grid sample budget exhausted");
  std::vector<double> values(gridSize);
  const auto          id    = [&](int x, int y, int z) { return (z * sy + y) * sx + x; };
  const auto          point = [&](int node) {
    const int x = node % sx, y = (node / sx) % sy, z = node / (sx * sy);
    return gp_Pnt(minimum.X() + extent.X() * (double(x) / nx),
                  minimum.Y() + extent.Y() * (double(y) / ny),
                  minimum.Z() + extent.Z() * (double(z) / nz));
  };
  for (size_t node = 0; node < gridSize; ++node)
  {
    values[node] = field.Value(point(int(node)));
    require(std::isfinite(values[node]), "Nonfinite native mesh field sample");
  }
  std::unordered_map<unsigned long long, int> roots;
  std::vector<bool>                           boundary(result.Nodes.size(), false);
  const auto                                  root = [&](int a, int b) {
    if (a > b)
      std::swap(a, b);
    if (values[size_t(a)] == 0)
      b = a;
    else if (values[size_t(b)] == 0)
      a = b;
    const auto edge  = (static_cast<unsigned long long>(a) << 32) | unsigned(b);
    const auto found = roots.find(edge);
    if (found != roots.end())
      return found->second;
    gp_Pnt     p = point(a), q = point(b);
    double     value = values[size_t(a)];
    const bool onBox =
      (a % sx == b % sx && (a % sx == 0 || a % sx == nx))
      || ((a / sx) % sy == (b / sx) % sy && ((a / sx) % sy == 0 || (a / sx) % sy == ny))
      || (a / (sx * sy) == b / (sx * sy) && (a / (sx * sy) == 0 || a / (sx * sy) == nz));
    for (int step = 0; step < 64 && p.Distance(q) > tolerance / 4096; ++step)
    {
      const gp_Pnt midpoint((p.XYZ() + q.XYZ()) * .5);
      const double v = field.Value(midpoint);
      if (v == 0)
      {
        p = q = midpoint;
        break;
      }
      if ((v < 0) == (value < 0))
      {
        p     = midpoint;
        value = v;
      }
      else
        q = midpoint;
    }
    require(result.Nodes.size() < size_t(maximumNodes), "Native mesh node budget exhausted");
    const int index = int(result.Nodes.size());
    result.Nodes.emplace_back((p.XYZ() + q.XYZ()) * .5);
    boundary.push_back(onBox);
    roots.emplace(edge, index);
    return index;
  };
  const auto triangle = [&](int a, int b, int c) {
    if (a == b || a == c || b == c)
      return;
    if (boundary[size_t(a)] || boundary[size_t(b)] || boundary[size_t(c)])
      return;
    const gp_Pnt center((result.Nodes[size_t(a)].XYZ() + result.Nodes[size_t(b)].XYZ()
                         + result.Nodes[size_t(c)].XYZ())
                        / 3);
    const auto   area = gp_Vec(result.Nodes[size_t(a)], result.Nodes[size_t(b)])
                          .Crossed(gp_Vec(result.Nodes[size_t(a)], result.Nodes[size_t(c)]));
    if (orientation * area.Dot(field.Gradient(center)) < 0)
      std::swap(b, c);
    require(result.ImplicitTriangles.size() + result.CapTriangles.size() < size_t(maximumTriangles),
            "Native mesh triangle budget exhausted");
    result.ImplicitTriangles.push_back({a, b, c});
  };
  // Consistent Freudenthal tetrahedra: adjacent cells share the same diagonals.
  constexpr int offsets[8][3] =
    {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
  constexpr int tets[6][4] =
    {{0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7}, {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}};
  for (int z = 0; z < nz; ++z)
    for (int y = 0; y < ny; ++y)
      for (int x = 0; x < nx; ++x)
      {
        int nodes[8];
        for (int i = 0; i < 8; ++i)
          nodes[i] = id(x + offsets[i][0], y + offsets[i][1], z + offsets[i][2]);
        for (const auto& tet : tets)
        {
          int negative[4], positive[4], n = 0, p = 0;
          for (int i : tet)
            if (values[size_t(nodes[i])] < 0)
              negative[n++] = nodes[i];
            else
              positive[p++] = nodes[i];
          if (n == 0 || p == 0)
            continue;
          if (n == 1 || p == 1)
          {
            const int  a     = n == 1 ? negative[0] : positive[0];
            const int* other = n == 1 ? positive : negative;
            triangle(root(a, other[0]), root(a, other[1]), root(a, other[2]));
          }
          else
          {
            const int a = root(negative[0], positive[0]), b = root(negative[0], positive[1]);
            const int c = root(negative[1], positive[0]), d = root(negative[1], positive[1]);
            triangle(a, b, c);
            triangle(b, d, c);
          }
        }
      }
  require(!result.ImplicitTriangles.empty(), "Empty native surface mesh candidate");
  auto loops = boundaryLoops(result.ImplicitTriangles);
  require(loops.size() == result.Wires.size(), "Native/interior mesh boundary count differs");
  size_t correspondenceCost = 0;
  for (const auto& wire : result.Wires)
    for (const auto& loop : loops)
    {
      correspondenceCost += wire.size() * loop.size();
      require(correspondenceCost <= 16000000, "Native mesh loop matching budget exhausted");
    }
  const double     gate = std::sqrt(std::pow(extent.X() / nx, 2) + std::pow(extent.Y() / ny, 2)
                                    + std::pow(extent.Z() / nz, 2))
                          * 4;
  std::set<size_t> paired;
  for (const auto& wire : result.Wires)
  {
    std::vector<size_t> candidates;
    for (size_t i = 0; i < loops.size(); ++i)
      if (pathDistance(wire, loops[i], result) <= gate)
        candidates.push_back(i);
    require(candidates.size() == 1 && paired.insert(candidates[0]).second,
            "Ambiguous native/interior mesh loop correspondence");
    collar(result,
           wire,
           loops[candidates[0]],
           field,
           hessian,
           tolerance,
           reserve,
           orientation,
           maximumTriangles);
  }
  refine(result,
         field,
         minimum,
         maximum,
         hessian,
         tolerance,
         reserve,
         orientation,
         maximumNodes,
         maximumTriangles);
  std::map<Key, Use> edges;
  for (const auto& t : result.ImplicitTriangles)
  {
    const double distance = bound(result, t, field, hessian, tolerance, reserve, orientation);
    require(distance <= tolerance, "Native surface mesh field-distance bound unresolved");
    result.FieldDistanceBound = std::max(result.FieldDistanceBound, distance);
    for (int i = 0; i < 3; ++i)
      add(edges, t[i], t[(i + 1) % 3]);
  }
  for (const auto& t : result.CapTriangles)
    for (int i = 0; i < 3; ++i)
      add(edges, t.Nodes[i], t.Nodes[(i + 1) % 3]);
  for (const auto& entry : edges)
    require(entry.second.Count == 2 && entry.second.Balance == 0,
            "Native candidate mesh is not closed and oriented");
  // Paired edges alone permit pinched/bow-tie vertices. Each used vertex
  // must have one connected cyclic link, not multiple fans sharing a node.
  std::vector<std::vector<Key>> links(result.Nodes.size());
  const auto                    linkTriangle = [&](const Triangle& t) {
    for (int i = 0; i < 3; ++i)
      links[size_t(t[i])].emplace_back(t[(i + 1) % 3], t[(i + 2) % 3]);
  };
  for (const auto& t : result.ImplicitTriangles)
    linkTriangle(t);
  for (const auto& t : result.CapTriangles)
    linkTriangle(t.Nodes);
  for (const auto& link : links)
  {
    if (link.empty())
      continue;
    std::map<int, std::vector<int>> adjacency;
    for (const auto& e : link)
    {
      adjacency[e.first].push_back(e.second);
      adjacency[e.second].push_back(e.first);
    }
    require(adjacency.size() >= 3, "Degenerate native mesh vertex link");
    for (const auto& entry : adjacency)
      require(entry.second.size() == 2, "Nonmanifold native mesh vertex link");
    std::set<int>    visited;
    std::vector<int> pending{adjacency.begin()->first};
    while (!pending.empty())
    {
      const int node = pending.back();
      pending.pop_back();
      if (!visited.insert(node).second)
        continue;
      for (int other : adjacency.at(node))
        pending.push_back(other);
    }
    require(visited.size() == adjacency.size(), "Pinched native mesh vertex");
  }
  // Drop unused grid roots without disturbing canonical cap/seam indices.
  std::vector<int> mapping(result.Nodes.size(), -1);
  for (size_t i = 0; i < prepared.Nodes.size(); ++i)
    mapping[i] = int(i);
  auto oldNodes = std::move(result.Nodes);
  result.Nodes  = prepared.Nodes;
  for (auto& t : result.ImplicitTriangles)
    for (int& node : t)
    {
      auto& mapped = mapping[size_t(node)];
      if (mapped < 0)
      {
        mapped = int(result.Nodes.size());
        result.Nodes.push_back(oldNodes[size_t(node)]);
      }
      node = mapped;
    }
  result.NativeSurfaceMeshed = true;
  result.MeshClosed          = true;
  result.Cells               = cells;
  return result;
}
} // namespace

BRepMesh_ImplicitBoundary::Result BRepMesh_ImplicitBoundary::Tessellate(
  const BRepClass3d_ImplicitBox::Boundary& source,
  double                                   tolerance,
  int                                      maximumCells,
  int                                      maximumNodes,
  int                                      maximumTriangles)
{
  if (!std::isfinite(tolerance) || tolerance < 1.e-7 || maximumCells < 8 || maximumCells > 256)
    throw Standard_DomainError("Invalid native surface mesh settings");
  const auto prepared = Prepare(source, tolerance, maximumNodes, maximumTriangles);
  const auto face     = occ::down_cast<BRep_TImplicitFace>(source.Faces.back().TShape());
  const auto field    = occ::down_cast<Geom_ImplicitGyroid>(face->Field());
  require(!field.IsNull(), "Native surface mesher requires the supported gyroid");
  const double orientation = source.Faces.back().Orientation() == TopAbs_REVERSED ? -1. : 1.;
  std::string  diagnostics;
  for (int cells = 8; cells <= maximumCells; cells *= 2)
  {
    try
    {
      return extract(prepared,
                     *face,
                     *field,
                     orientation,
                     tolerance,
                     cells,
                     maximumNodes,
                     maximumTriangles);
    }
    catch (const Standard_NotImplemented& error)
    {
      diagnostics += std::to_string(cells) + ": " + error.what() + "; ";
      if (cells * 2 > maximumCells)
        throw Standard_NotImplemented(diagnostics.c_str());
    }
  }
  throw Standard_NotImplemented("Native surface meshing limits exhausted");
}
