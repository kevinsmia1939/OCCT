// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitCylinder.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_Plane.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Circle.hxx>
#include <NCollection_Array1.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS_Wire.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_AbortiveTransaction.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace
{
using Settings = BRep_ImplicitCylinderSection::Settings;
using Caps = BRepClass3d_ImplicitCylinder::Caps;
void require(bool value, const char* reason) { if (!value) throw Standard_DomainError(reason); }
void checkpoint(const Settings& s) { if (s.WorkCheckpoint) s.WorkCheckpoint(); }
struct Link { TopoDS_Edge Edge; double Area = 0, Length = 0, Tolerance = 0; };
struct Cycle
{
  TopoDS_Wire Wire;
  TopoDS_Face Interior;
  gp_Pnt Witness;
  double Area = 0, Reserve = 0;
  int Parent = -1;
};

void buildCap(const TopoDS_Face& source, const Settings& s, int cap, bool negative, Caps& result,
              size_t& nestingTests)
{
  const auto& disk = result.Source.Disks.Disks[cap];
  const auto& roots = result.Source.Disks.Boundary.Roots;
  const auto& curveRoots = result.Source.Disks.CurveRoots[cap];
  const auto origin = s.Frame.Location().Translated(gp_Vec(s.Frame.Direction())
                       * (cap == 0 ? s.HeightMinimum : s.HeightMaximum));
  const gp_Dir normal = cap == 0 ? s.Frame.Direction().Reversed() : s.Frame.Direction();
  const gp_Ax2 frame(origin, normal, s.Frame.XDirection());
  const occ::handle<Geom_Plane> plane = new Geom_Plane(frame);
  const gp_Vec x(frame.XDirection()), y(frame.YDirection()), n(normal);
  const auto uv = [&](const gp_Pnt& p) {
    const gp_Vec d(origin, p); return gp_Pnt2d(d.Dot(x), d.Dot(y));
  };
  const double magnitude = 1 + std::max({std::abs(origin.X()), std::abs(origin.Y()),
                                         std::abs(origin.Z()), s.Radius});
  const double precision = 2048 * std::numeric_limits<double>::epsilon() * magnitude;
  require(precision < s.DistanceTolerance, "cap_coordinate_precision");
  const auto data = occ::down_cast<BRep_TImplicitFace>(source.TShape());
  const auto* field = static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  const double scale = std::abs((source.Location() * data->Location()).Transformation().ScaleFactor());
  const double k = 2 * std::acos(-1.) / field->Period() / scale;
  const double hessian = 4 * k * k;
  BRep_Builder builder;
  std::vector<Link> links;
  for (const auto& arc : result.Source.Rings.Rings[cap])
  {
    checkpoint(s);
    if (arc.Negative != negative) continue;
    double first, last; BRep_Tool::Range(arc.Edge, first, last);
    const occ::handle<Geom2d_Circle> pc = new Geom2d_Circle(
      gp_Ax2d(gp_Pnt2d(0, 0), gp_Dir2d(1, 0)), s.Radius, cap == 1);
    builder.UpdateEdge(arc.Edge, pc, plane, TopLoc_Location(), BRep_Tool::Tolerance(arc.Edge));
    builder.Range(arc.Edge, plane, TopLoc_Location(), first, last);
    builder.SameRange(arc.Edge, true); builder.SameParameter(arc.Edge, true);
    const double direction = arc.Edge.Orientation() == TopAbs_REVERSED ? -1 : 1;
    const double sense = cap == 0 ? -1 : 1;
    links.push_back({arc.Edge, .5 * s.Radius * s.Radius * (last - first) * direction * sense,
                     s.Radius * (last - first), BRep_Tool::Tolerance(arc.Edge)});
  }
  for (size_t index = 0; index < disk.Curves.size(); ++index)
  {
    checkpoint(s);
    const auto& curve = disk.Curves[index];
    const int count = int(curve.Points.size());
    require(count >= 2, "cap_empty_contour");
    NCollection_Array1<gp_Pnt> poles(1, count);
    NCollection_Array1<gp_Pnt2d> poles2d(1, count);
    NCollection_Array1<double> knots(1, count);
    NCollection_Array1<int> mults(1, count);
    double area = 0, length = 0;
    int direction = 0;
    for (int i = 1; i <= count; ++i)
    {
      const auto& p = curve.Points[size_t(i - 1)];
      poles.SetValue(i, p); poles2d.SetValue(i, uv(p));
      knots.SetValue(i, i - 1); mults.SetValue(i, i == 1 || i == count ? 2 : 1);
      if (i == 1) continue;
      const auto& a = curve.Points[size_t(i - 2)];
      const gp_Vec tangent(a, p);
      const double distance = tangent.Magnitude();
      require(distance > precision, "cap_collapsed_contour_segment");
      const gp_Pnt middle((a.XYZ() + p.XYZ()) * .5);
      const auto gradient = BRep_ImplicitTool::Gradient(source, middle);
      const double orientation = tangent.Dot(n.Crossed(gradient)) * (negative ? 1 : -1);
      const double uncertainty = hessian * disk.GeometricBound * distance
        + 2048 * std::numeric_limits<double>::epsilon() * (1 + gradient.Magnitude() * distance);
      require(std::abs(orientation) > uncertainty, "cap_contour_orientation_unresolved");
      const int sign = orientation > 0 ? 1 : -1;
      require(direction == 0 || direction == sign, "cap_inconsistent_contour_orientation");
      direction = sign;
      const auto pa = uv(a), pb = uv(p);
      area += .5 * (pa.X() * pb.Y() - pb.X() * pa.Y()); length += distance;
    }
    TopoDS_Vertex first, last;
    const auto& ends = curveRoots[index];
    if (curve.Closed)
    {
      require(ends[0] == 0 && ends[1] == 0, "cap_closed_contour_root_ids");
      builder.MakeVertex(first, curve.Points.front(), precision); last = first;
    }
    else
    {
      require(ends[0] && ends[1] && ends[0] <= roots.size() && ends[1] <= roots.size(), "cap_missing_root");
      first = roots[ends[0] - 1].Vertex; last = roots[ends[1] - 1].Vertex;
    }
    const double tolerance = std::max({precision, BRep_Tool::Tolerance(first), BRep_Tool::Tolerance(last)});
    require(tolerance <= s.DistanceTolerance, "cap_edge_tolerance");
    const occ::handle<Geom_BSplineCurve> geometry = new Geom_BSplineCurve(poles, knots, mults, 1);
    const occ::handle<Geom2d_BSplineCurve> pc = new Geom2d_BSplineCurve(poles2d, knots, mults, 1);
    TopoDS_Edge edge;
    builder.MakeEdge(edge, geometry, tolerance);
    first.Orientation(TopAbs_FORWARD); last.Orientation(TopAbs_REVERSED);
    builder.Add(edge, first); builder.Add(edge, last);
    builder.UpdateVertex(first, 0, edge, BRep_Tool::Tolerance(first));
    builder.UpdateVertex(last, count - 1, edge, BRep_Tool::Tolerance(last));
    builder.Range(edge, 0, count - 1);
    builder.UpdateEdge(edge, pc, plane, TopLoc_Location(), tolerance);
    builder.Range(edge, plane, TopLoc_Location(), 0, count - 1);
    builder.SameRange(edge, true); builder.SameParameter(edge, true);
    if (direction < 0) edge.Reverse();
    result.FieldEdges[cap].push_back(edge);
    links.push_back({edge, area * direction, length, tolerance});
  }

  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> vertices;
  std::vector<int> begin, end;
  for (const auto& link : links)
  {
    TopoDS_Vertex a, b; TopExp::Vertices(link.Edge, a, b, true);
    require(!a.IsNull() && !b.IsNull(), "cap_edge_without_endpoints");
    begin.push_back(vertices.Add(a) - 1); end.push_back(vertices.Add(b) - 1);
  }
  std::vector<int> incoming(size_t(vertices.Size()), 0), outgoing(size_t(vertices.Size()), -1);
  for (size_t i = 0; i < links.size(); ++i)
  {
    require(outgoing[size_t(begin[i])] == -1, "cap_branched_directed_cycle");
    outgoing[size_t(begin[i])] = int(i); ++incoming[size_t(end[i])];
  }
  for (size_t i = 0; i < incoming.size(); ++i)
    require(incoming[i] == 1 && outgoing[i] >= 0, "cap_open_or_branched_cycle");
  std::vector<bool> used(links.size(), false);
  std::vector<Cycle> cycles;
  for (size_t first = 0; first < links.size(); ++first)
  {
    if (used[first]) continue;
    Cycle cycle; builder.MakeWire(cycle.Wire);
    double length = 0, tolerance = precision;
    size_t current = first;
    do
    {
      checkpoint(s);
      require(!used[current], "cap_cycle_joined_foreign_cycle"); used[current] = true;
      builder.Add(cycle.Wire, links[current].Edge);
      cycle.Area += links[current].Area; length += links[current].Length;
      tolerance = std::max(tolerance, links[current].Tolerance);
      current = size_t(outgoing[size_t(end[current])]);
    } while (current != first);
    cycle.Wire.Closed(true);
    cycle.Reserve = (precision + tolerance) * length;
    require(std::isfinite(cycle.Area) && std::isfinite(cycle.Reserve)
              && std::isfinite(length) && std::abs(cycle.Area) > cycle.Reserve,
            "cap_cycle_area_unresolved");
    double a, b; const auto geometry = BRep_Tool::Curve(links[first].Edge, a, b);
    cycle.Witness = geometry->Value(a + (b - a) * .5);
    auto positiveWire = cycle.Wire;
    if (cycle.Area < 0) positiveWire.Reverse();
    BRepBuilderAPI_MakeFace interior(plane, positiveWire, true);
    require(interior.IsDone(), "cap_cycle_face_construction");
    cycle.Interior = interior.Face();
    require(BRepCheck_Analyzer(cycle.Interior).IsValid(), "cap_cycle_invalid_wire_geometry");
    cycles.push_back(std::move(cycle));
  }
  for (size_t i = 0; i < cycles.size(); ++i)
  {
    double parentArea = std::numeric_limits<double>::infinity();
    for (size_t j = 0; j < cycles.size(); ++j)
    {
      checkpoint(s);
      require(++nestingTests <= result.MaximumNestingTests, "cap_nesting_budget");
      if (i == j || std::abs(cycles[j].Area) <= std::abs(cycles[i].Area)) continue;
      BRepClass_FaceClassifier classifier(cycles[j].Interior, cycles[i].Witness, precision);
      require(classifier.State() != TopAbs_ON, "cap_touching_nested_cycles");
      if (classifier.State() == TopAbs_IN && std::abs(cycles[j].Area) < parentArea)
      { parentArea = std::abs(cycles[j].Area); cycles[i].Parent = int(j); }
    }
    if (cycles[i].Parent < 0) require(cycles[i].Area > 0, "cap_hole_without_outer_cycle");
    else require((cycles[i].Area > 0) != (cycles[size_t(cycles[i].Parent)].Area > 0),
                 "cap_inconsistent_material_nesting");
  }
  for (size_t i = 0; i < cycles.size(); ++i)
  {
    if (cycles[i].Area < 0) continue;
    BRepBuilderAPI_MakeFace maker(plane, cycles[i].Wire, true);
    for (const auto& hole : cycles)
      if (hole.Parent == int(i) && hole.Area < 0) maker.Add(hole.Wire);
    require(maker.IsDone(), "cap_region_construction");
    const auto face = maker.Face();
    require(BRepCheck_Analyzer(face).IsValid(), "cap_region_invalid_topology");
    result.Faces[cap].push_back(face);
  }
}
} // namespace

static Caps buildCapsFromSkeleton(const TopoDS_Face& face, const Settings& s, bool negative,
                                  size_t maximumNestingTests,
                                  BRep_ImplicitCylinderSection::CapSkeleton skeleton)
{
  Caps result;
  require(maximumNestingTests > 0 && maximumNestingTests <= 10000000, "invalid_cap_nesting_budget");
  result.MaximumNestingTests = maximumNestingTests;
  result.Field = face; result.Parameters = s; result.Parameters.WorkCheckpoint = nullptr;
  result.Negative = negative;
  result.Source = std::move(skeleton);
  if (!result.Source.Verified) { result.Reason = "unresolved_cap_skeleton"; return result; }
  result.ContourBound = std::max(result.Source.Disks.Disks[0].GeometricBound,
                                 result.Source.Disks.Disks[1].GeometricBound);
  try
  {
    size_t nestingTests = 0;
    for (int cap = 0; cap < 2; ++cap) buildCap(face, s, cap, negative, result, nestingTests);
    result.TopologyValidated = true; result.Reason = "cap_topology_validated";
  }
  catch (const Standard_AbortiveTransaction&) { throw; }
  catch (const Standard_Failure& error)
  {
    result.Reason = error.what() ? error.what() : "cap_construction_failed";
    for (int cap = 0; cap < 2; ++cap) { result.Faces[cap].clear(); result.FieldEdges[cap].clear(); }
  }
  return result;
}

Caps BRepClass3d_ImplicitCylinder::BuildCaps(const TopoDS_Face& face, const Settings& s, bool negative,
                                           size_t maximumNestingTests)
{
  require(maximumNestingTests > 0 && maximumNestingTests <= 10000000, "invalid_cap_nesting_budget");
  return buildCapsFromSkeleton(face, s, negative, maximumNestingTests,
    BRep_ImplicitCylinderSection::BuildCapSkeleton(face, s));
}

BRepClass3d_ImplicitCylinder::ChartCaps BRepClass3d_ImplicitCylinder::BuildChartCaps(
  const TopoDS_Face& face, const Settings& s, bool negative, size_t maximumNestingTests)
{
  require(maximumNestingTests > 0 && maximumNestingTests <= 10000000, "invalid_cap_nesting_budget");
  auto source = BRep_ImplicitCylinderSection::BuildChartCapSkeleton(face, s);
  ChartCaps result;
  result.Chart = std::move(source.Chart);
  result.Regions = buildCapsFromSkeleton(face, s, negative, maximumNestingTests, std::move(source.Caps));
  return result;
}

BRepClass3d_ImplicitCylinder::SideBoundary BRepClass3d_ImplicitCylinder::BuildSideBoundary(
  const TopoDS_Face& face, const Settings& s, bool negative, size_t maximumNestingTests)
{
  SideBoundary result;
  result.Source = BuildChartCaps(face, s, negative, maximumNestingTests);
  if (!result.Source.Regions.TopologyValidated)
  { result.Reason = result.Source.Regions.Reason; return result; }
  const auto& graph = result.Source.Chart.Graph;
  const auto surface = result.Source.Chart.Surface;
  const double turn = 2 * std::acos(-1.);
  const double magnitude = 1 + s.Radius + std::abs(s.HeightMinimum) + std::abs(s.HeightMaximum)
    + s.Frame.Location().Distance(gp_Pnt());
  const double precision = 4096 * std::numeric_limits<double>::epsilon() * magnitude
    * (1 + std::abs(s.AngleMinimum) + std::abs(s.AngleMaximum));
  BRep_Builder builder;
  try
  {
    require(std::isfinite(precision) && precision < s.DistanceTolerance, "side_coordinate_precision");
    const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
    const auto* field = static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
    const double scale = std::abs((face.Location() * data->Location()).Transformation().ScaleFactor());
    const double k = 2 * std::acos(-1.) / field->Period() / scale;
    const double hessian = 4 * k * k;
    const double lipschitz = field->LipschitzBound() / scale;
    for (size_t c = 0; c < graph.CurveVertices.size(); ++c)
    {
      checkpoint(s);
      auto ids = graph.CurveVertices[c];
      require(ids.size() >= 2, "side_empty_contour");
      const bool closed = graph.Boundary.Section.Curves[c].Closed;
      if (closed)
      {
        ids.pop_back();
        const auto seam = std::find_if(ids.begin(), ids.end(), [&](size_t id) {
          return graph.Vertices[id].Seam;
        });
        if (seam != ids.end()) std::rotate(ids.begin(), seam, ids.end());
        ids.push_back(ids.front());
      }
      size_t begin = 0;
      while (begin + 1 < ids.size())
      {
        size_t end = begin + 1;
        while (end + 1 < ids.size() && !graph.Vertices[ids[end]].Seam) ++end;
        const int count = int(end - begin + 1);
        std::vector<gp_Pnt2d> uv;
        for (size_t i = begin; i <= end; ++i)
        {
          auto p = graph.Vertices[ids[i]].UV;
          if (!uv.empty()) p.SetX(p.X() + turn * std::round((uv.back().X() - p.X()) / turn));
          uv.push_back(p);
        }
        double mean = 0;
        for (const auto& p : uv) mean += (p.X() - s.AngleMinimum) / count;
        const double shift = turn * std::floor(mean / turn);
        for (auto& p : uv)
        {
          p.SetX(p.X() - shift);
          require(p.X() >= s.AngleMinimum - precision / s.Radius
                    && p.X() <= s.AngleMaximum + precision / s.Radius,
                  "side_unsplit_chart_seam");
        }
        NCollection_Array1<gp_Pnt> poles(1, count);
        NCollection_Array1<gp_Pnt2d> poles2d(1, count);
        NCollection_Array1<double> knots(1, count);
        NCollection_Array1<int> mults(1, count);
        int direction = 0;
        double deviation = 0;
        for (int i = 1; i <= count; ++i)
        {
          checkpoint(s);
          const auto& p = graph.Vertices[ids[begin + size_t(i - 1)]].Point;
          poles.SetValue(i, p); poles2d.SetValue(i, uv[size_t(i - 1)]);
          knots.SetValue(i, i - 1); mults.SetValue(i, i == 1 || i == count ? 2 : 1);
          if (i == 1) continue;
          const auto& a = poles.Value(i - 1);
          const gp_Vec tangent(a, p);
          const double length = tangent.Magnitude();
          require(length > precision, "side_collapsed_segment");
          const double angle = (uv[size_t(i - 2)].X() + uv[size_t(i - 1)].X()) * .5;
          const double delta = std::abs(uv[size_t(i - 1)].X() - uv[size_t(i - 2)].X());
          deviation = std::max(deviation, s.Radius * delta * delta / 8 + precision);
          const gp_Vec normal = gp_Vec(s.Frame.XDirection()) * std::cos(angle)
                              + gp_Vec(s.Frame.YDirection()) * std::sin(angle);
          const auto gradient = BRep_ImplicitTool::Gradient(face, gp_Pnt((a.XYZ() + p.XYZ()) * .5));
          const double orientation = tangent.Dot(normal.Crossed(gradient)) * (negative ? 1 : -1);
          const double uncertainty = length * (hessian * graph.Boundary.Section.GeometricBound
            + lipschitz * delta * .5) + precision * (1 + gradient.Magnitude());
          require(std::abs(orientation) > uncertainty, "side_orientation_unresolved");
          const int sign = orientation > 0 ? 1 : -1;
          require(direction == 0 || direction == sign, "side_inconsistent_orientation");
          direction = sign;
        }
        auto first = graph.Vertices[ids[begin]].Vertex, last = graph.Vertices[ids[end]].Vertex;
        if (first.IsNull())
        {
          require(closed && begin == 0 && end + 1 == ids.size(), "side_missing_endpoint");
          builder.MakeVertex(first, poles.Value(1), precision); last = first;
        }
        require(!last.IsNull(), "side_missing_endpoint");
        const double tolerance = std::max({deviation, BRep_Tool::Tolerance(first), BRep_Tool::Tolerance(last)});
        require(tolerance <= s.DistanceTolerance, "side_edge_tolerance");
        const occ::handle<Geom_BSplineCurve> geometry = new Geom_BSplineCurve(poles, knots, mults, 1);
        const occ::handle<Geom2d_BSplineCurve> pc = new Geom2d_BSplineCurve(poles2d, knots, mults, 1);
        TopoDS_Edge edge; builder.MakeEdge(edge, geometry, tolerance);
        first.Orientation(TopAbs_FORWARD); last.Orientation(TopAbs_REVERSED);
        builder.Add(edge, first); builder.Add(edge, last);
        builder.UpdateVertex(first, 0, edge, BRep_Tool::Tolerance(first));
        builder.UpdateVertex(last, count - 1, edge, BRep_Tool::Tolerance(last));
        builder.Range(edge, 0, count - 1);
        builder.UpdateEdge(edge, pc, surface, TopLoc_Location(), tolerance);
        builder.Range(edge, surface, TopLoc_Location(), 0, count - 1);
        builder.SameRange(edge, true); builder.SameParameter(edge, true);
        if (direction < 0) edge.Reverse();
        result.FieldEdges.push_back(edge);
        result.SurfaceDeviationBound = std::max(result.SurfaceDeviationBound, deviation);
        begin = end;
      }
    }
    result.Verified = true; result.Reason = "side_boundary_verified";
  }
  catch (const Standard_AbortiveTransaction&) { throw; }
  catch (const Standard_Failure& error)
  {
    result.FieldEdges.clear(); result.SurfaceDeviationBound = 0;
    result.Reason = error.what() ? error.what() : "side_boundary_failed";
  }
  return result;
}

BRepClass3d_ImplicitCylinder::Sides BRepClass3d_ImplicitCylinder::BuildSides(
  const TopoDS_Face& source, const Settings& s, bool negative, size_t maximumNestingTests)
{
  Sides result;
  result.Boundary = BuildSideBoundary(source, s, negative, maximumNestingTests);
  if (!result.Boundary.Verified) { result.Reason = result.Boundary.Reason; return result; }
  const auto& chart = result.Boundary.Source.Chart;
  const auto surface = chart.Surface;
  BRep_Builder builder;
  struct SideLink { TopoDS_Edge Edge; double Area, Length, Tolerance; gp_Pnt2d Witness; int Begin, End; };
  struct SideCycle { TopoDS_Wire Wire; TopoDS_Face Interior; gp_Pnt2d Witness; double Area; int Parent = -1; };
  try
  {
    const double precision = 8192 * std::numeric_limits<double>::epsilon()
      * (1 + s.Radius + std::abs(s.HeightMinimum) + std::abs(s.HeightMaximum)
         + s.Frame.Location().Distance(gp_Pnt()))
      * (1 + std::abs(s.AngleMinimum) + std::abs(s.AngleMaximum));
    require(precision < s.DistanceTolerance, "side_region_precision");
    std::vector<TopoDS_Edge> edges = result.Boundary.FieldEdges;
    for (const auto& ring : chart.Rings.Rings)
      for (const auto& arc : ring) if (arc.Negative == negative)
      { auto edge = arc.Edge; edge.Reverse(); edges.push_back(edge); }
    for (const auto& edge : chart.SeamEdges)
    {
      checkpoint(s);
      double first, last; BRep_Tool::Range(edge, first, last);
      const auto p = surface->Value(s.AngleMinimum, (first + last) * .5);
      const double value = BRep_ImplicitTool::Value(source, p);
      const auto gradient = BRep_ImplicitTool::Gradient(source, p);
      require(std::isfinite(value) && std::abs(value) > precision * (1 + gradient.Magnitude()),
              "side_seam_material_unresolved");
      if ((value < 0) == negative)
      { edges.push_back(edge); auto reverse = edge; reverse.Reverse(); edges.push_back(reverse); }
    }
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> vertices;
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> seamVertices;
    seamVertices.Add(chart.SeamEnds[0]); seamVertices.Add(chart.SeamEnds[1]);
    for (const auto& node : chart.Graph.Vertices)
    { checkpoint(s); if (node.Seam) seamVertices.Add(node.Vertex); }
    std::map<std::pair<int, int>, int> occurrences;
    const auto occurrence = [&](const TopoDS_Vertex& v, const gp_Pnt2d& uv) {
      const bool seam = seamVertices.Contains(v);
      const int branch = seam ? (std::abs(uv.X() - s.AngleMinimum)
                              < std::abs(uv.X() - s.AngleMaximum) ? 1 : 2) : 0;
      const auto key = std::make_pair(vertices.Add(v), branch);
      const auto found = occurrences.find(key);
      if (found != occurrences.end()) return found->second;
      const int id = int(occurrences.size()); occurrences.emplace(key, id); return id;
    };
    std::vector<SideLink> links;
    for (const auto& edge : edges)
    {
      checkpoint(s);
      double first, last;
      const auto pc = BRep_Tool::CurveOnSurface(edge, surface, TopLoc_Location(), first, last);
      require(!pc.IsNull(), "side_missing_pcurve");
      std::vector<gp_Pnt2d> points;
      const auto bs = occ::down_cast<Geom2d_BSplineCurve>(pc);
      if (!bs.IsNull())
        for (int i = 1; i <= bs->NbPoles(); ++i) points.push_back(bs->Pole(i));
      else { points.push_back(pc->Value(first)); points.push_back(pc->Value(last)); }
      if (edge.Orientation() == TopAbs_REVERSED) std::reverse(points.begin(), points.end());
      double area = 0, length = 0;
      for (size_t i = 1; i < points.size(); ++i)
      {
        const auto& a = points[i - 1]; const auto& b = points[i];
        area += .5 * s.Radius * (a.X() * b.Y() - b.X() * a.Y());
        length += std::hypot(s.Radius * (b.X() - a.X()), b.Y() - a.Y());
      }
      TopoDS_Vertex a, b; TopExp::Vertices(edge, a, b, true);
      require(!a.IsNull() && !b.IsNull(), "side_missing_vertices");
      links.push_back({edge, area, length, BRep_Tool::Tolerance(edge),
        pc->Value((first + last) * .5), occurrence(a, points.front()), occurrence(b, points.back())});
    }
    std::vector<int> outgoing(occurrences.size(), -1), incoming(occurrences.size(), 0);
    for (size_t i = 0; i < links.size(); ++i)
    {
      require(outgoing[size_t(links[i].Begin)] == -1, "side_branched_material_cycle");
      outgoing[size_t(links[i].Begin)] = int(i); ++incoming[size_t(links[i].End)];
    }
    for (size_t i = 0; i < incoming.size(); ++i)
      require(incoming[i] == 1 && outgoing[i] >= 0, "side_open_material_cycle");
    std::vector<bool> used(links.size(), false);
    std::vector<SideCycle> cycles;
    for (size_t first = 0; first < links.size(); ++first)
    {
      if (used[first]) continue;
      SideCycle cycle; builder.MakeWire(cycle.Wire); cycle.Area = 0; cycle.Witness = links[first].Witness;
      size_t current = first;
      double length = 0, tolerance = precision;
      do
      {
        checkpoint(s);
        require(!used[current], "side_joined_foreign_cycle"); used[current] = true;
        const auto& link = links[current]; builder.Add(cycle.Wire, link.Edge);
        cycle.Area += link.Area; length += link.Length; tolerance = std::max(tolerance, link.Tolerance);
        current = size_t(outgoing[size_t(link.End)]);
      } while (current != first);
      cycle.Wire.Closed(true);
      require(std::isfinite(cycle.Area) && std::isfinite(length)
                && std::abs(cycle.Area) > (precision + tolerance) * length,
              "side_cycle_area_unresolved");
      auto positive = cycle.Wire; if (cycle.Area < 0) positive.Reverse();
      builder.MakeFace(cycle.Interior, surface, precision); builder.Add(cycle.Interior, positive);
      require(BRepCheck_Analyzer(cycle.Interior).IsValid(), "side_cycle_invalid_topology");
      cycles.push_back(cycle);
    }
    size_t tests = 0;
    for (size_t i = 0; i < cycles.size(); ++i)
    {
      double parentArea = std::numeric_limits<double>::infinity();
      for (size_t j = 0; j < cycles.size(); ++j)
      {
        checkpoint(s); require(++tests <= maximumNestingTests, "side_nesting_budget");
        if (i == j || std::abs(cycles[j].Area) <= std::abs(cycles[i].Area)) continue;
        const auto state = BRepClass_FaceClassifier(cycles[j].Interior, cycles[i].Witness, precision).State();
        require(state != TopAbs_ON, "side_touching_cycles");
        if (state == TopAbs_IN && std::abs(cycles[j].Area) < parentArea)
        { parentArea = std::abs(cycles[j].Area); cycles[i].Parent = int(j); }
      }
      if (cycles[i].Parent < 0) require(cycles[i].Area > 0, "side_hole_without_outer_cycle");
      else require((cycles[i].Area > 0) != (cycles[size_t(cycles[i].Parent)].Area > 0),
                    "side_inconsistent_material_nesting");
    }
    for (size_t i = 0; i < cycles.size(); ++i)
    {
      if (cycles[i].Area < 0) continue;
      TopoDS_Face face; builder.MakeFace(face, surface, precision); builder.Add(face, cycles[i].Wire);
      for (const auto& hole : cycles)
        if (hole.Parent == int(i) && hole.Area < 0) builder.Add(face, hole.Wire);
      require(BRepCheck_Analyzer(face).IsValid(), "side_region_invalid_topology");
      result.Faces.push_back(face);
    }
    result.TopologyValidated = true; result.Reason = "side_topology_validated";
  }
  catch (const Standard_AbortiveTransaction&) { throw; }
  catch (const Standard_Failure& error)
  { result.Faces.clear(); result.Reason = error.what() ? error.what() : "side_region_failed"; }
  return result;
}

BRepClass3d_ImplicitCylinder::TrimBoundary BRepClass3d_ImplicitCylinder::BuildTrimBoundary(
  const TopoDS_Face& face, const Settings& s, bool negative, size_t maximumNestingTests)
{
  TrimBoundary result;
  result.Source = BuildSides(face, s, negative, maximumNestingTests);
  if (!result.Source.TopologyValidated) { result.Reason = result.Source.Reason; return result; }
  try
  {
    std::vector<TopoDS_Edge> field = result.Source.Boundary.FieldEdges;
    for (const auto& cap : result.Source.Boundary.Source.Regions.FieldEdges)
      field.insert(field.end(), cap.begin(), cap.end());
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> expected, allEdges, vertices;
    for (const auto& edge : field)
    {
      checkpoint(s);
      require(!expected.Contains(edge), "cylinder_duplicate_field_edge"); expected.Add(edge);
    }
    struct Uses { int Forward = 0, Reversed = 0; };
    std::vector<Uses> uses;
    std::vector<TopoDS_Face> conventional = result.Source.Faces;
    for (const auto& cap : result.Source.Boundary.Source.Regions.Faces)
      conventional.insert(conventional.end(), cap.begin(), cap.end());
    for (const auto& conventionalFace : conventional)
      for (TopExp_Explorer it(conventionalFace, TopAbs_EDGE); it.More(); it.Next())
      {
        checkpoint(s);
        const auto& edge = it.Current();
        const int id = allEdges.Add(edge) - 1;
        if (size_t(id) >= uses.size()) uses.resize(size_t(id) + 1);
        require(edge.Orientation() == TopAbs_FORWARD || edge.Orientation() == TopAbs_REVERSED,
                "cylinder_nonmanifold_edge_orientation");
        if (edge.Orientation() == TopAbs_FORWARD) ++uses[size_t(id)].Forward;
        else ++uses[size_t(id)].Reversed;
      }
    for (size_t id = 1; id <= allEdges.Size(); ++id)
    {
      checkpoint(s);
      const auto& use = uses[size_t(id - 1)];
      const int fieldId = expected.FindIndex(allEdges.FindKey(id));
      if (fieldId)
      {
        const bool forward = field[size_t(fieldId - 1)].Orientation() == TopAbs_FORWARD;
        require(use.Forward == (forward ? 1 : 0) && use.Reversed == (forward ? 0 : 1),
                "cylinder_field_edge_incidence");
      }
      else require(use.Forward == 1 && use.Reversed == 1, "cylinder_unpaired_conventional_edge");
    }
    std::vector<int> begin, end;
    for (auto& edge : field)
    {
      require(allEdges.Contains(edge), "cylinder_unused_field_edge");
      edge.Reverse();
      TopoDS_Vertex a, b; TopExp::Vertices(edge, a, b, true);
      require(!a.IsNull() && !b.IsNull(), "cylinder_trim_missing_vertex");
      begin.push_back(vertices.Add(a) - 1); end.push_back(vertices.Add(b) - 1);
    }
    std::vector<int> outgoing(size_t(vertices.Size()), -1), incoming(size_t(vertices.Size()), 0);
    for (size_t i = 0; i < field.size(); ++i)
    {
      require(outgoing[size_t(begin[i])] < 0, "cylinder_trim_branched_vertex");
      outgoing[size_t(begin[i])] = int(i); ++incoming[size_t(end[i])];
    }
    for (size_t i = 0; i < incoming.size(); ++i)
      require(incoming[i] == 1 && outgoing[i] >= 0, "cylinder_trim_open_vertex");
    result.EdgeWires.assign(field.size(), -1);
    BRep_Builder builder;
    for (size_t first = 0; first < field.size(); ++first)
    {
      if (result.EdgeWires[first] >= 0) continue;
      TopoDS_Wire wire; builder.MakeWire(wire);
      size_t current = first;
      do
      {
        checkpoint(s);
        require(result.EdgeWires[current] < 0, "cylinder_trim_joined_foreign_cycle");
        result.EdgeWires[current] = int(result.Wires.size());
        builder.Add(wire, field[current]);
        current = size_t(outgoing[size_t(end[current])]);
      } while (current != first);
      wire.Closed(true);
      require(BRepCheck_Analyzer(wire).IsValid(), "cylinder_trim_invalid_wire");
      result.Wires.push_back(wire);
    }
    result.IncidenceVerified = true; result.Reason = "cylinder_trim_incidence_verified";
  }
  catch (const Standard_AbortiveTransaction&) { throw; }
  catch (const Standard_Failure& error)
  {
    result.Wires.clear(); result.EdgeWires.clear();
    result.Reason = error.what() ? error.what() : "cylinder_trim_incidence_failed";
  }
  return result;
}
