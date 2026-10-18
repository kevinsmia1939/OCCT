// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitBox.hxx>
#include <BRepClass3d_ImplicitBoxIO.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_Tool.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Standard_DomainError.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <sstream>

namespace
{
using Edge = std::pair<int, int>;

struct Group
{
  int              Side = -1, Shell = -1, Region = -1;
  std::vector<int> Triangles;
};

int capSide(const TopoDS_Face& face)
{
  const auto normal = BRepAdaptor_Surface(face).Plane().Axis().Direction();
  int        axis   = 0;
  for (int j = 1; j < 3; ++j)
    if (std::abs(normal.Coord(j + 1)) > std::abs(normal.Coord(axis + 1)))
      axis = j;
  return 2 * axis + int(normal.Coord(axis + 1) > 0);
}

double segmentDistance(const gp_Pnt& p, const gp_Pnt& a, const gp_Pnt& b)
{
  const gp_Vec ab(a, b);
  const double length = ab.SquareMagnitude();
  const double t      = length > 0 ? std::clamp(gp_Vec(a, p).Dot(ab) / length, 0., 1.) : 0;
  return p.Distance(a.Translated(ab * t));
}
} // namespace

BRepClass3d_ImplicitBox::CapAttachment BRepClass3d_ImplicitBox::CheckCapAttachment(
  const Boundary& source,
  int             resolution,
  int             maximumDepth,
  size_t          maximumCells,
  size_t          maximumTests)
{
  if (resolution < 4 || resolution > 64 || (resolution & (resolution - 1)) || maximumDepth < 0
      || maximumDepth > 16 || !maximumCells || maximumCells > 10000000 || !maximumTests
      || maximumTests > 10000000)
    throw Standard_DomainError("Invalid cap attachment witness budget");
  CapAttachment report;
  const auto    graph        = CheckCorrespondence(source);
  report.StoredGraphResolved = graph.Resolved;
  if (!graph.Resolved)
  {
    report.Reason = "stored_graph_correspondence_unresolved: " + graph.Reason;
    return report;
  }
  std::ostringstream record;
  BRepClass3d_ImplicitBoxIO::Write(source, record);
  std::istringstream input(record.str());
  Boundary           copy;
  BRepClass3d_ImplicitBoxIO::Read(input, copy);
  auto reference = CheckRegionTopology(copy.Faces.back(), resolution, maximumDepth, maximumCells, copy.Policy);
  report.RegionReferenceResolved = reference.ClosedReferenceManifold;
  if (!report.RegionReferenceResolved)
  {
    report.Reason = "region_reference_unresolved: " + reference.Reason;
    return report;
  }
  resolution = reference.Surface.Resolution;
  report.ReferenceShellCount             = int(reference.ComponentEuler.size());
  report.ReferenceRegionCount            = reference.RegionCount;
  const int                        count = int(reference.Triangles.size());
  std::vector<int>                 parent(size_t(count), -1), depth(size_t(count), -1);
  std::vector<std::vector<int>>    neighbors(size_t(count), std::vector<int>{});
  std::map<Edge, std::vector<int>> uses;
  for (int i = 0; i < count; ++i)
    if (reference.TriangleSides[size_t(i)] >= 0)
    {
      parent[size_t(i)] = i;
      const auto& t     = reference.Triangles[size_t(i)];
      for (int j = 0; j < 3; ++j)
        uses[std::minmax(t[size_t(j)], t[size_t((j + 1) % 3)])].push_back(i);
    }
  auto root = [&](int i) {
    while (parent[size_t(i)] != i)
    {
      parent[size_t(i)] = parent[size_t(parent[size_t(i)])];
      i                 = parent[size_t(i)];
    }
    return i;
  };
  std::queue<int> queue;
  for (const auto& entry : uses)
  {
    const auto& adjacent = entry.second;
    if (adjacent.size() == 2
        && reference.TriangleSides[size_t(adjacent[0])]
             == reference.TriangleSides[size_t(adjacent[1])])
    {
      parent[size_t(root(adjacent[0]))] = root(adjacent[1]);
      neighbors[size_t(adjacent[0])].push_back(adjacent[1]);
      neighbors[size_t(adjacent[1])].push_back(adjacent[0]);
    }
    else
      for (int i : adjacent)
        if (depth[size_t(i)] < 0)
        {
          depth[size_t(i)] = 0;
          queue.push(i);
        }
  }
  while (!queue.empty())
  {
    const int i = queue.front();
    queue.pop();
    for (int j : neighbors[size_t(i)])
      if (depth[size_t(j)] < 0)
      {
        depth[size_t(j)] = depth[size_t(i)] + 1;
        queue.push(j);
      }
  }
  std::map<int, Group> groups;
  for (int i = 0; i < count; ++i)
    if (parent[size_t(i)] >= 0)
    {
      auto&     group = groups[root(i)];
      const int shell = reference.TriangleComponents[size_t(i)];
      if (group.Side >= 0
          && (group.Side != reference.TriangleSides[size_t(i)] || group.Shell != shell))
      {
        report.Reason = "reference_cap_component_inconsistent";
        return report;
      }
      group.Side   = reference.TriangleSides[size_t(i)];
      group.Shell  = shell;
      group.Region = reference.ComponentRegions[size_t(shell)];
      group.Triangles.push_back(i);
    }
  report.ReferenceCapCount = int(groups.size());
  report.ReferenceTriangleCaps.assign(size_t(count), -1);
  const size_t caps = copy.Faces.size() - 1;
  if (caps != groups.size())
  {
    report.Reason = "stored_reference_cap_region_counts_differ";
    return report;
  }
  report.ReferenceCaps.assign(caps, -1);
  report.ReferenceShells.assign(caps, -1);
  report.ReferenceRegions.assign(caps, -1);
  report.Witnesses.resize(caps);
  report.Clearances.resize(caps, 0);
  const auto  data  = occ::down_cast<BRep_TImplicitFace>(copy.Faces.back().TShape());
  const auto& field = *static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  const auto  low = data->Minimum(), high = data->Maximum();
  double      magnitude = 1;
  for (int a = 1; a <= 3; ++a)
    magnitude = std::max({magnitude, std::abs(low.Coord(a)), std::abs(high.Coord(a))});
  const double reserve = 8192 * std::numeric_limits<double>::epsilon() * magnitude;
  const double fieldReserve =
    4096 * std::numeric_limits<double>::epsilon() * (1 + std::abs(field.IsoValue()))
    + field.LipschitzBound() * reserve;
  report.RequiredClearance = graph.ReferenceSeamBound + graph.CurveDeviationBound
                             + graph.CapParameterDeviationBound + graph.MatchingTolerance;
  report.MatchingTolerance = graph.MatchingTolerance;
  std::vector<std::vector<std::pair<gp_Pnt, gp_Pnt>>> segments(caps);
  std::vector<int>                                    sides;
  for (size_t cap = 0; cap < caps; ++cap)
  {
    sides.push_back(capSide(copy.Faces[cap]));
    for (TopExp_Explorer it(copy.Faces[cap], TopAbs_EDGE); it.More(); it.Next())
    {
      TopLoc_Location location;
      double          first, last;
      const auto      curve = occ::down_cast<Geom_BSplineCurve>(
        BRep_Tool::Curve(TopoDS::Edge(it.Current()), location, first, last));
      for (int pole = 1; pole < curve->NbPoles(); ++pole)
        segments[cap].emplace_back(curve->Pole(pole), curve->Pole(pole + 1));
    }
  }
  const bool negative    = copy.Faces.back().Orientation() == TopAbs_FORWARD;
  auto       linearValue = [&](const gp_Pnt& p, int side) {
    const int    axis = side / 2, u = (axis + 1) % 3, v = (axis + 2) % 3;
    const double su =
      (p.Coord(u + 1) - low.Coord(u + 1)) * resolution / (high.Coord(u + 1) - low.Coord(u + 1));
    const double sv =
      (p.Coord(v + 1) - low.Coord(v + 1)) * resolution / (high.Coord(v + 1) - low.Coord(v + 1));
    const int    iu = std::clamp(int(std::floor(su)), 0, resolution - 1),
                 iv = std::clamp(int(std::floor(sv)), 0, resolution - 1);
    const double fu = std::clamp(su - iu, 0., 1.), fv = std::clamp(sv - iv, 0., 1.);
    auto         sample = [&](int du, int dv) {
      auto q = low;
      q.SetCoord(axis + 1, side % 2 ? high.Coord(axis + 1) : low.Coord(axis + 1));
      q.SetCoord(u + 1,
                 iu + du == resolution
                   ? high.Coord(u + 1)
                   : low.Coord(u + 1)
                       + (high.Coord(u + 1) - low.Coord(u + 1)) * (iu + du) / resolution);
      q.SetCoord(v + 1,
                 iv + dv == resolution
                   ? high.Coord(v + 1)
                   : low.Coord(v + 1)
                       + (high.Coord(v + 1) - low.Coord(v + 1)) * (iv + dv) / resolution);
      double f = field.Value(q);
      if (copy.Policy == CornerPolicy::Regularized)
      {
        bool corner = true;
        for (int axis = 1; axis <= 3; ++axis)
          corner &= q.Coord(axis) == low.Coord(axis) || q.Coord(axis) == high.Coord(axis);
        double exact = 0;
        gp_Vec gradient;
        if (corner && field.LatticeJet(q, exact, gradient) && exact == 0)
          f = 0; // Same fixed node as L, never a rounded near-zero substitution.
      }
      return f;
    };
    return fu >= fv ? (1 - fu) * sample(0, 0) + (fu - fv) * sample(1, 0) + fv * sample(1, 1)
                    : (1 - fv) * sample(0, 0) + (fv - fu) * sample(0, 1) + fu * sample(1, 1);
  };
  int groupIndex = 0;
  for (auto& entry : groups)
  {
    auto& group = entry.second;
    report.ReferenceCapSides.push_back(group.Side);
    for (int index : group.Triangles)
      report.ReferenceTriangleCaps[size_t(index)] = groupIndex;
    std::stable_sort(group.Triangles.begin(), group.Triangles.end(), [&](int a, int b) {
      return depth[size_t(a)] > depth[size_t(b)];
    });
    int matched = -1;
    for (int index : group.Triangles)
    {
      if (report.WitnessTests >= maximumTests)
      {
        report.Reason = "cap_witness_budget_exhausted";
        return report;
      }
      ++report.WitnessTests;
      const auto& t = reference.Triangles[size_t(index)];
      gp_Pnt      p((reference.Nodes[size_t(t[0])].XYZ() + reference.Nodes[size_t(t[1])].XYZ()
                     + reference.Nodes[size_t(t[2])].XYZ())
                    / 3);
      p.SetCoord(group.Side / 2 + 1,
                 group.Side % 2 ? high.Coord(group.Side / 2 + 1) : low.Coord(group.Side / 2 + 1));
      const double sign = negative ? -1 : 1;
      if (!(sign * field.Value(p) > fieldReserve
            && sign * linearValue(p, group.Side) > fieldReserve))
        continue;
      for (size_t cap = 0; cap < caps; ++cap)
      {
        if (sides[cap] != group.Side
            || BRepClass_FaceClassifier(copy.Faces[cap], p, reserve).State() != TopAbs_IN)
          continue;
        double clearance = std::numeric_limits<double>::infinity();
        for (const auto& edge : segments[cap])
        {
          if (report.SegmentTests >= maximumTests)
          {
            report.Reason = "cap_clearance_budget_exhausted";
            return report;
          }
          ++report.SegmentTests;
          clearance = std::min(clearance, segmentDistance(p, edge.first, edge.second) - reserve);
          if (clearance <= report.RequiredClearance)
            break;
        }
        if (!(clearance > report.RequiredClearance))
          continue;
        if (matched >= 0 || report.ReferenceCaps[cap] >= 0)
        {
          report.Reason = "cap_attachment_ambiguous_or_nonbijective";
          return report;
        }
        matched                      = int(cap);
        report.ReferenceCaps[cap]    = groupIndex;
        report.ReferenceShells[cap]  = group.Shell;
        report.ReferenceRegions[cap] = group.Region;
        report.Witnesses[cap]        = p;
        report.Clearances[cap]       = clearance;
      }
      if (matched >= 0)
        break;
    }
    if (matched < 0)
    {
      report.Reason = "cap_common_retained_witness_unresolved";
      return report;
    }
    ++groupIndex;
  }
  report.Resolved  = true;
  report.Reason    = "cap_regions_attached";
  report.Reference = std::move(reference);
  return report;
}
