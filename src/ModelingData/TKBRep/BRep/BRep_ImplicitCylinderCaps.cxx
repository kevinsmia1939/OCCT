// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRep_ImplicitCylinderSection.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_Tool.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace
{
using Section = BRep_ImplicitCylinderSection;
void checkpoint(const Section::Settings& s) { if (s.WorkCheckpoint) s.WorkCheckpoint(); }

bool traceDisk(const TopoDS_Face& face, const Section::Settings& s, int cap,
               const Section::BoundaryResult& boundary, Section::Result& result,
               std::vector<std::array<size_t, 2>>& curveRoots)
{
  BRep_ImplicitPlaneSection::Settings plane;
  const double height = cap == 0 ? s.HeightMinimum : s.HeightMaximum;
  plane.Frame = gp_Ax2(s.Frame.Location().Translated(gp_Vec(s.Frame.Direction()) * height),
                       s.Frame.Direction(), s.Frame.XDirection());
  plane.UMinimum = plane.VMinimum = -s.Radius;
  plane.UMaximum = plane.VMaximum = s.Radius;
  plane.DistanceTolerance = s.DistanceTolerance * .5;
  plane.MinimumCells = s.MinimumCells; plane.MaximumCells = s.MaximumCells;
  plane.WorkCheckpoint = s.WorkCheckpoint;
  const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  const auto transform = (face.Location() * data->Location()).Transformation();
  bool exact = true;
  for (int row = 1; row <= 3; ++row) for (int col = 1; col <= 4; ++col)
    exact &= transform.Value(row, col) == (row == col ? 1. : 0.);
  for (const auto& dir : {plane.Frame.XDirection(), plane.Frame.YDirection()})
    for (int axis = 1; axis <= 3; ++axis)
      exact &= dir.Coord(axis) == 0 || std::abs(dir.Coord(axis)) == 1;
  plane.AllowExactLatticeNodes = exact;
  const auto graph = BRep_ImplicitPlaneSection::TraceGraph(face, plane);
  result = graph.Section; result.Curves.clear(); result.PointContacts.clear();
  if (!graph.Section.Resolved) return false;
  const gp_Vec x(plane.Frame.XDirection()), y(plane.Frame.YDirection());
  const auto xy = [&](const gp_Pnt& p) {
    const gp_Vec delta(plane.Frame.Location(), p);
    return std::array<double, 2> {{delta.Dot(x), delta.Dot(y)}};
  };
  const double magnitude = 1 + std::max({std::abs(plane.Frame.Location().X()),
    std::abs(plane.Frame.Location().Y()), std::abs(plane.Frame.Location().Z()), s.Radius});
  const double reserve = 1024 * std::numeric_limits<double>::epsilon() * magnitude;
  result.GeometricBound += 2 * reserve;
  if (!std::isfinite(result.GeometricBound) || result.GeometricBound > s.DistanceTolerance)
  { result.Resolved = false; ++result.AmbiguousCells; return false; }
  for (const auto& point : graph.Section.PointContacts)
    if (std::hypot(xy(point)[0], xy(point)[1]) <= s.Radius + reserve)
    { ++result.AmbiguousCells; result.Resolved = false; return false; }
  const int n = graph.Section.Resolution;
  const double step = 2 * s.Radius / n;
  using Key = std::pair<int, int>;
  std::map<Key, std::vector<size_t>> cellRoots;
  std::vector<gp_Pnt> points = graph.Vertices;
  std::vector<size_t> ids(points.size(), 0), uses(boundary.Roots.size(), 0);
  std::vector<size_t> rootVertex(boundary.Roots.size(), 0);
  for (size_t index = 0; index < boundary.Roots.size(); ++index)
  {
    const auto& root = boundary.Roots[index];
    if (root.HeightBoundary != cap) continue;
    const auto point = BRep_Tool::Pnt(root.Vertex);
    const auto coord = xy(point);
    const double error = BRep_Tool::Tolerance(root.Vertex) + reserve;
    std::array<int, 2> cell;
    for (int axis = 0; axis < 2; ++axis)
    {
      const double low = std::floor((coord[axis] - error + s.Radius) / step);
      const double high = std::floor((coord[axis] + error + s.Radius) / step);
      if (low != high || low < 0 || high >= n)
      { ++result.AmbiguousCells; result.Resolved = false; return false; }
      cell[axis] = int(low);
    }
    rootVertex[index] = points.size(); points.push_back(point); ids.push_back(index + 1);
    cellRoots[{cell[0], cell[1]}].push_back(index);
  }
  std::vector<std::pair<size_t, size_t>> segments;
  for (const auto& arc : graph.Arcs)
  {
    checkpoint(s);
    const double lowX = -s.Radius + step * arc.U, highX = lowX + step;
    const double lowY = -s.Radius + step * arc.V, highY = lowY + step;
    const double far = std::hypot(std::max(std::abs(lowX), std::abs(highX)),
                                  std::max(std::abs(lowY), std::abs(highY)));
    const double near = std::hypot(std::max({lowX, -highX, 0.}), std::max({lowY, -highY, 0.}));
    const auto found = cellRoots.find({arc.U, arc.V});
    if (far < s.Radius - reserve)
    {
      if (found != cellRoots.end()) { ++result.GraphDefects; result.Resolved = false; return false; }
      segments.emplace_back(arc.First, arc.Last); continue;
    }
    if (near > s.Radius + reserve)
    {
      if (found != cellRoots.end()) { ++result.GraphDefects; result.Resolved = false; return false; }
      continue;
    }
    std::vector<size_t> retained;
    for (size_t vertex : {arc.First, arc.Last})
    {
      const auto coord = xy(points[vertex]);
      const double radial = std::hypot(coord[0], coord[1]);
      const double error = graph.VertexRadius + reserve;
      if (radial + error < s.Radius) retained.push_back(vertex);
      else if (!(radial - error > s.Radius))
      { ++result.AmbiguousCells; result.Resolved = false; return false; }
    }
    if (found != cellRoots.end()) for (size_t index : found->second)
    { retained.push_back(rootVertex[index]); ++uses[index]; }
    if (retained.size() == 2) segments.emplace_back(retained[0], retained[1]);
    else if (!retained.empty())
    { ++result.GraphDefects; result.Resolved = false; return false; }
  }
  for (size_t index = 0; index < boundary.Roots.size(); ++index)
    if (boundary.Roots[index].HeightBoundary == cap && uses[index] != 1)
    { ++result.GraphDefects; result.Resolved = false; return false; }
  std::vector<std::vector<size_t>> adjacent(points.size());
  for (size_t index = 0; index < segments.size(); ++index)
  { adjacent[segments[index].first].push_back(index); adjacent[segments[index].second].push_back(index); }
  for (size_t index = 0; index < adjacent.size(); ++index)
    if ((!adjacent[index].empty() && adjacent[index].size() != (ids[index] ? 1 : 2))
        || (ids[index] && adjacent[index].empty())) ++result.GraphDefects;
  if (result.GraphDefects) { result.Resolved = false; return false; }
  std::vector<bool> used(segments.size(), false);
  const auto walk = [&](size_t start, size_t next) {
    BRep_ImplicitPlaneSection::Curve curve;
    size_t current = start; curve.Points.push_back(points[start]);
    while (!used[next])
    {
      checkpoint(s); used[next] = true;
      current = segments[next].first == current ? segments[next].second : segments[next].first;
      curve.Points.push_back(points[current]);
      if (current == start) { curve.Closed = true; break; }
      bool found = false;
      for (auto candidate : adjacent[current]) if (!used[candidate])
      { next = candidate; found = true; break; }
      if (!found) break;
    }
    curveRoots.push_back({ids[start], ids[current]});
    result.Curves.push_back(std::move(curve));
  };
  for (size_t i = 0; i < adjacent.size(); ++i)
    if (adjacent[i].size() == 1 && !used[adjacent[i][0]]) walk(i, adjacent[i][0]);
  for (size_t i = 0; i < segments.size(); ++i) if (!used[i]) walk(segments[i].first, i);
  return true;
}
} // namespace

Section::DiskResult Section::TraceEndCapDisks(const TopoDS_Face& face, const Settings& s)
{
  DiskResult result;
  result.Boundary = TracePeriodicBoundary(face, s);
  if (!result.Boundary.BoundaryIncidenceVerified) return result;
  for (int cap = 0; cap < 2; ++cap)
    if (!traceDisk(face, s, cap, result.Boundary, result.Disks[cap], result.CurveRoots[cap]))
    {
      for (int i = 0; i < 2; ++i)
      { result.Disks[i].Resolved = false; result.Disks[i].Curves.clear(); result.CurveRoots[i].clear(); }
      return result;
    }
  result.Verified = true;
  return result;
}
