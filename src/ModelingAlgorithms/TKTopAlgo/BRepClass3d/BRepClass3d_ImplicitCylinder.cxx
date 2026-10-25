// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitCylinder.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_TImplicitCylinderPatch.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Circle.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Circle.hxx>
#include <Geom2d_Line.hxx>
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
#include <set>
#include <TopoDS_Iterator.hxx>

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

BRepClass3d_ImplicitCylinder::AxialSheet BRepClass3d_ImplicitCylinder::CheckAxialSheet(
  const TopoDS_Face& face, const Settings& s, bool negative, size_t maximumNestingTests,
  int maximumDepth, size_t maximumCells)
{
  require(maximumDepth >= 0 && maximumDepth <= 32 && maximumCells > 0 && maximumCells <= 2000000,
          "invalid_axial_sheet_budget");
  AxialSheet result;
  result.Boundary = BuildTrimBoundary(face, s, negative, maximumNestingTests);
  if (!result.Boundary.IncidenceVerified) { result.Reason = result.Boundary.Reason; return result; }
  try
  {
    const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
    const auto* field = static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
    const double scale = std::abs((face.Location() * data->Location()).Transformation().ScaleFactor());
    const double k = 2 * std::acos(-1.) / field->Period() / scale;
    const double hessian = 4 * k * k, lipschitz = field->LipschitzBound() / scale;
    const double magnitude = 1 + s.Radius + std::abs(s.HeightMinimum) + std::abs(s.HeightMaximum)
      + s.Frame.Location().Distance(gp_Pnt());
    const double precision = 4096 * std::numeric_limits<double>::epsilon() * magnitude;
    const double reserve = 4096 * std::numeric_limits<double>::epsilon() * (1 + hessian + lipschitz)
      + (hessian + lipschitz) * precision;
    require(std::isfinite(reserve) && std::isfinite(hessian), "axial_sheet_precision");
    const gp_Vec axis(s.Frame.Direction()), x(s.Frame.XDirection()), y(s.Frame.YDirection());
    const auto point = [&](const gp_Pnt& local) {
      return s.Frame.Location().Translated(x * local.X() + y * local.Y() + axis * local.Z());
    };
    const double middle = s.HeightMinimum + (s.HeightMaximum - s.HeightMinimum) * .5;
    const double derivative = BRep_ImplicitTool::Gradient(face, point(gp_Pnt(0, 0, middle))).Dot(axis);
    require(std::isfinite(derivative) && std::abs(derivative) > reserve, "axial_sheet_center_tangent");
    const int direction = derivative > 0 ? 1 : -1;
    struct Cell { gp_Pnt Low, High; int Depth; };
    std::vector<Cell> pending {{gp_Pnt(-s.Radius, -s.Radius, s.HeightMinimum),
                               gp_Pnt(s.Radius, s.Radius, s.HeightMaximum), 0}};
    double lowerBound = std::numeric_limits<double>::infinity();
    while (!pending.empty())
    {
      checkpoint(s);
      require(result.Cells < maximumCells, "axial_sheet_cell_budget"); ++result.Cells;
      const auto cell = pending.back(); pending.pop_back();
      // Only cells strictly separated from the disk are skipped. Near-boundary
      // boxes remain, preserving a conservative cover with coordinate reserves.
      const auto nearest = [](double low, double high) {
        return low > 0 ? low : (high < 0 ? -high : 0.);
      };
      const double radial = std::hypot(nearest(cell.Low.X(), cell.High.X()),
                                       nearest(cell.Low.Y(), cell.High.Y()));
      if (radial > s.Radius + precision) continue;
      const gp_Pnt center((cell.Low.XYZ() + cell.High.XYZ()) * .5);
      const gp_Vec half(cell.Low, cell.High);
      const double radius = half.Magnitude() * .5 + precision;
      const double value = direction * BRep_ImplicitTool::Gradient(face, point(center)).Dot(axis);
      const double uncertainty = hessian * radius + reserve;
      require(std::isfinite(value) && std::isfinite(uncertainty), "axial_sheet_nonfinite_cone");
      const double bound = value - uncertainty;
      if (bound > 0) { lowerBound = std::min(lowerBound, bound); continue; }
      require(value + uncertainty > 0, "axial_sheet_direction_changes");
      require(cell.Depth < maximumDepth, "axial_sheet_depth_budget");
      int split = 1;
      for (int a = 2; a <= 3; ++a) if (half.Coord(a) > half.Coord(split)) split = a;
      const double cut = cell.Low.Coord(split) + half.Coord(split) * .5;
      require(cut > cell.Low.Coord(split) && cut < cell.High.Coord(split), "axial_sheet_split_precision");
      auto first = cell, second = cell;
      first.High.SetCoord(split, cut); second.Low.SetCoord(split, cut);
      first.Depth = second.Depth = cell.Depth + 1;
      pending.push_back(second); pending.push_back(first);
    }
    require(std::isfinite(lowerBound) && lowerBound > 0, "axial_sheet_missing_cover");
    result.Direction = direction; result.DerivativeLowerBound = lowerBound;
    result.MonotonicityVerified = true; result.Reason = "axial_sheet_monotonicity_verified";
  }
  catch (const Standard_AbortiveTransaction&) { throw; }
  catch (const Standard_Failure& error)
  { result.Reason = error.what() ? error.what() : "axial_sheet_unresolved"; }
  return result;
}

BRepClass3d_ImplicitCylinder::ProjectedRegions BRepClass3d_ImplicitCylinder::BuildProjectedRegions(
  const TopoDS_Face& face, const Settings& s, bool negative, size_t maximumNestingTests,
  int maximumDepth, size_t maximumCells)
{
  ProjectedRegions result;
  result.Source = CheckAxialSheet(face, s, negative, maximumNestingTests, maximumDepth, maximumCells);
  if (!result.Source.MonotonicityVerified) { result.Reason = result.Source.Reason; return result; }
  try
  {
    const auto& boundary = result.Source.Boundary;
    const auto& chart = boundary.Source.Boundary.Source.Chart;
    const occ::handle<Geom_Plane> plane = new Geom_Plane(s.Frame);
    const occ::handle<Geom_Circle> circle = new Geom_Circle(s.Frame, s.Radius);
    const gp_Vec x(s.Frame.XDirection()), y(s.Frame.YDirection());
    const auto uv = [&](const gp_Pnt& p) {
      const gp_Vec d(s.Frame.Location(), p); return gp_Pnt2d(d.Dot(x), d.Dot(y));
    };
    const auto project = [&](const gp_Pnt& p) { const auto q = uv(p); return plane->Value(q.X(), q.Y()); };
    const double precision = 8192 * std::numeric_limits<double>::epsilon()
      * (1 + s.Radius + s.Frame.Location().Distance(gp_Pnt()));
    require(std::isfinite(precision) && precision < s.DistanceTolerance, "projected_coordinate_precision");
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> sideEdges, nativeVertices;
    for (const auto& edge : boundary.Source.Boundary.FieldEdges) sideEdges.Add(edge);
    std::vector<TopoDS_Vertex> projectedVertices;
    BRep_Builder builder;
    const auto vertex = [&](const TopoDS_Vertex& v) {
      const int id = nativeVertices.Add(v) - 1;
      if (size_t(id) == projectedVertices.size())
      {
        TopoDS_Vertex q; builder.MakeVertex(q, project(BRep_Tool::Pnt(v)),
                                           std::max(precision, BRep_Tool::Tolerance(v)));
        projectedVertices.push_back(q);
      }
      return projectedVertices[size_t(id)];
    };
    std::vector<Cycle> cycles;
    double maximumTolerance = precision;
    const int sense = result.Source.Direction * (negative ? 1 : -1);
    for (const auto& wire : boundary.Wires)
    {
      Cycle cycle; builder.MakeWire(cycle.Wire);
      double length = 0, tolerance = precision;
      bool witness = false;
      for (TopExp_Explorer it(wire, TopAbs_EDGE); it.More(); it.Next())
      {
        checkpoint(s);
        const auto original = TopoDS::Edge(it.Current());
        TopoDS_Vertex a, b; TopExp::Vertices(original, a, b);
        auto firstVertex = vertex(a), lastVertex = vertex(b);
        double first, last, area = 0;
        TopoDS_Edge edge;
        int direction = original.Orientation() == TopAbs_REVERSED ? -1 : 1;
        const double tol = std::max({precision, BRep_Tool::Tolerance(original),
                                     BRep_Tool::Tolerance(firstVertex), BRep_Tool::Tolerance(lastVertex)});
        if (sideEdges.Contains(original))
        {
          const auto pc = occ::down_cast<Geom2d_BSplineCurve>(BRep_Tool::CurveOnSurface(
            original, chart.Surface, TopLoc_Location(), first, last));
          require(!pc.IsNull() && pc->Degree() == 1, "projected_missing_side_correspondence");
          first = pc->Pole(1).X(); last = pc->Pole(pc->NbPoles()).X();
          const int angular = last > first ? 1 : -1;
          require(std::abs(last - first) * s.Radius > precision, "projected_collapsed_side_arc");
          for (int i = 2; i <= pc->NbPoles(); ++i)
            require(angular * (pc->Pole(i).X() - pc->Pole(i - 1).X()) * s.Radius > precision,
                    "projected_nonmonotone_side_arc");
          if (angular < 0) { std::swap(first, last); std::swap(firstVertex, lastVertex); }
          direction *= angular;
          builder.MakeEdge(edge, circle, tol);
          const occ::handle<Geom2d_Circle> curve2d = new Geom2d_Circle(
            gp_Ax2d(gp_Pnt2d(0, 0), gp_Dir2d(1, 0)), s.Radius);
          builder.UpdateEdge(edge, curve2d, plane, TopLoc_Location(), tol);
          area = .5 * s.Radius * s.Radius * (last - first); length += s.Radius * (last - first);
        }
        else
        {
          const auto originalCurve = occ::down_cast<Geom_BSplineCurve>(BRep_Tool::Curve(original, first, last));
          require(!originalCurve.IsNull() && originalCurve->Degree() == 1
                    && originalCurve->NbKnots() == originalCurve->NbPoles(), "projected_cap_curve_type");
          const int count = originalCurve->NbPoles();
          NCollection_Array1<gp_Pnt> poles(1, count);
          NCollection_Array1<gp_Pnt2d> poles2d(1, count);
          NCollection_Array1<double> knots(1, count);
          NCollection_Array1<int> mults(1, count);
          for (int i = 1; i <= count; ++i)
          {
            poles.SetValue(i, project(originalCurve->Pole(i))); poles2d.SetValue(i, uv(originalCurve->Pole(i)));
            knots.SetValue(i, originalCurve->Knot(i)); mults.SetValue(i, originalCurve->Multiplicity(i));
            if (i > 1)
            {
              const auto p = poles2d.Value(i - 1), q = poles2d.Value(i);
              area += .5 * (p.X() * q.Y() - q.X() * p.Y()); length += p.Distance(q);
            }
          }
          const occ::handle<Geom_BSplineCurve> geometry = new Geom_BSplineCurve(poles, knots, mults, 1);
          const occ::handle<Geom2d_BSplineCurve> curve2d = new Geom2d_BSplineCurve(poles2d, knots, mults, 1);
          builder.MakeEdge(edge, geometry, tol);
          builder.UpdateEdge(edge, curve2d, plane, TopLoc_Location(), tol);
        }
        firstVertex.Orientation(TopAbs_FORWARD); lastVertex.Orientation(TopAbs_REVERSED);
        builder.Add(edge, firstVertex); builder.Add(edge, lastVertex);
        builder.UpdateVertex(firstVertex, first, edge, BRep_Tool::Tolerance(firstVertex));
        builder.UpdateVertex(lastVertex, last, edge, BRep_Tool::Tolerance(lastVertex));
        builder.Range(edge, first, last); builder.Range(edge, plane, TopLoc_Location(), first, last);
        builder.SameRange(edge, true); builder.SameParameter(edge, true);
        if (direction * sense < 0) edge.Reverse();
        builder.Add(cycle.Wire, edge); cycle.Area += area * direction * sense;
        tolerance = std::max(tolerance, tol);
        if (!witness)
        { cycle.Witness = BRep_Tool::Curve(edge, first, last)->Value(first + (last - first) * .5); witness = true; }
      }
      cycle.Wire.Closed(true); cycle.Reserve = (precision + tolerance) * length;
      maximumTolerance = std::max(maximumTolerance, tolerance);
      require(witness && std::isfinite(cycle.Area) && std::isfinite(cycle.Reserve)
                && std::abs(cycle.Area) > cycle.Reserve, "projected_area_unresolved");
      auto positive = cycle.Wire; if (cycle.Area < 0) positive.Reverse();
      BRepBuilderAPI_MakeFace maker(plane, positive, true);
      require(maker.IsDone(), "projected_cycle_construction"); cycle.Interior = maker.Face();
      require(BRepCheck_Analyzer(cycle.Interior).IsValid(), "projected_invalid_cycle");
      cycles.push_back(cycle);
    }
    size_t tests = 0;
    for (size_t i = 0; i < cycles.size(); ++i)
    {
      double parentArea = std::numeric_limits<double>::infinity();
      for (size_t j = 0; j < cycles.size(); ++j)
      {
        checkpoint(s); require(++tests <= maximumNestingTests, "projected_nesting_budget");
        if (j > i)
        {
          BRepExtrema_DistShapeShape distance(cycles[i].Wire, cycles[j].Wire);
          require(distance.IsDone() && distance.Value() > 2 * maximumTolerance,
                  "projected_intersecting_or_uncertain_cycles");
        }
        if (i == j || std::abs(cycles[j].Area) <= std::abs(cycles[i].Area)) continue;
        const auto state = BRepClass_FaceClassifier(cycles[j].Interior, cycles[i].Witness, precision).State();
        require(state != TopAbs_ON, "projected_touching_cycles");
        if (state == TopAbs_IN && std::abs(cycles[j].Area) < parentArea)
        { parentArea = std::abs(cycles[j].Area); cycles[i].Parent = int(j); }
      }
      if (cycles[i].Parent < 0) require(cycles[i].Area > 0, "projected_hole_without_outer");
      else require((cycles[i].Area > 0) != (cycles[size_t(cycles[i].Parent)].Area > 0),
                    "projected_inconsistent_nesting");
    }
    result.WireRegions.assign(cycles.size(), -1);
    for (size_t i = 0; i < cycles.size(); ++i)
    {
      if (cycles[i].Area < 0) continue;
      const int region = int(result.Faces.size()); result.WireRegions[i] = region;
      BRepBuilderAPI_MakeFace maker(plane, cycles[i].Wire, true);
      for (size_t j = 0; j < cycles.size(); ++j)
        if (cycles[j].Parent == int(i) && cycles[j].Area < 0)
        { maker.Add(cycles[j].Wire); result.WireRegions[j] = region; }
      require(maker.IsDone(), "projected_region_construction");
      require(BRepCheck_Analyzer(maker.Face()).IsValid(), "projected_region_invalid");
      result.Faces.push_back(maker.Face());
    }
    for (int id : result.WireRegions) require(id >= 0, "projected_unassigned_wire");
    result.ProxyTopologyValidated = true; result.Reason = "projected_proxy_topology_validated";
  }
  catch (const Standard_AbortiveTransaction&) { throw; }
  catch (const Standard_Failure& error)
  { result.Faces.clear(); result.WireRegions.clear(); result.Reason = error.what() ? error.what() : "projected_proxy_failed"; }
  return result;
}

BRepClass3d_ImplicitCylinder::ProjectedRegions BRepClass3d_ImplicitCylinder::CheckProjectedCorrespondence(
  const TopoDS_Face& face, const Settings& s, bool negative, size_t maximumNestingTests,
  int maximumDepth, size_t maximumCells)
{
  auto result = BuildProjectedRegions(face, s, negative, maximumNestingTests, maximumDepth, maximumCells);
  if (!result.ProxyTopologyValidated) return result;
  try
  {
    const auto& boundary = result.Source.Boundary;
    const auto& sides = boundary.Source.Boundary;
    require(result.Faces.size() == 1 && boundary.Wires.size() == 1
              && sides.FieldEdges.size() == 1
              && sides.Source.Regions.FieldEdges[0].empty() && sides.Source.Regions.FieldEdges[1].empty(),
            "projected_full_disk_not_qualified");
    const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
    const auto* field = static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
    const double scale = std::abs((face.Location() * data->Location()).Transformation().ScaleFactor());
    const double lipschitz = field->LipschitzBound() / scale;
    const double precision = 8192 * std::numeric_limits<double>::epsilon()
      * (1 + s.Radius + std::abs(s.HeightMinimum) + std::abs(s.HeightMaximum)
         + s.Frame.Location().Distance(gp_Pnt()));
    const double reserve = precision * (1 + lipschitz)
      + 8192 * std::numeric_limits<double>::epsilon() * (1 + std::abs(field->IsoValue()));
    require(std::isfinite(reserve), "projected_endpoint_precision");
    double first, last;
    const auto pc = BRep_Tool::CurveOnSurface(sides.FieldEdges.front(), sides.Source.Chart.Surface,
                                             TopLoc_Location(), first, last);
    require(!pc.IsNull() && std::abs(std::abs(pc->Value(last).X() - pc->Value(first).X())
              - 2 * std::acos(-1.)) * s.Radius < precision, "projected_incomplete_circle");
    struct Cell { double X0, X1, Y0, Y1; int Depth; };
    std::vector<Cell> pending {{-s.Radius, s.Radius, -s.Radius, s.Radius, 0}};
    const gp_Vec x(s.Frame.XDirection()), y(s.Frame.YDirection()), axis(s.Frame.Direction());
    while (!pending.empty())
    {
      checkpoint(s);
      require(result.CorrespondenceCells < maximumCells, "projected_endpoint_cell_budget");
      ++result.CorrespondenceCells;
      const auto cell = pending.back(); pending.pop_back();
      const auto nearest = [](double low, double high) { return low > 0 ? low : (high < 0 ? -high : 0.); };
      if (std::hypot(nearest(cell.X0, cell.X1), nearest(cell.Y0, cell.Y1)) > s.Radius + precision) continue;
      const double dx = cell.X1 - cell.X0, dy = cell.Y1 - cell.Y0;
      const auto center = s.Frame.Location().Translated(x * (cell.X0 + dx * .5) + y * (cell.Y0 + dy * .5));
      const double lo = result.Source.Direction * BRep_ImplicitTool::Value(face, center.Translated(axis * s.HeightMinimum));
      const double hi = result.Source.Direction * BRep_ImplicitTool::Value(face, center.Translated(axis * s.HeightMaximum));
      const double uncertainty = lipschitz * (std::hypot(dx, dy) * .5 + precision) + reserve;
      require(std::isfinite(lo) && std::isfinite(hi) && std::isfinite(uncertainty), "projected_endpoint_nonfinite");
      if (lo < -uncertainty && hi > uncertainty) continue;
      require(lo - uncertainty < 0 && hi + uncertainty > 0, "projected_endpoint_signs_fail");
      require(cell.Depth < maximumDepth, "projected_endpoint_depth_budget");
      auto a = cell, b = cell; a.Depth = b.Depth = cell.Depth + 1;
      if (dx >= dy) { const double cut = cell.X0 + dx * .5; a.X1 = b.X0 = cut; }
      else { const double cut = cell.Y0 + dy * .5; a.Y1 = b.Y0 = cut; }
      pending.push_back(b); pending.push_back(a);
    }
    result.MathematicalCorrespondenceVerified = true;
    result.Reason = "projected_full_disk_correspondence_verified";
  }
  catch (const Standard_AbortiveTransaction&) { throw; }
  catch (const Standard_Failure& error)
  { result.Reason = error.what() ? error.what() : "projected_correspondence_unresolved"; }
  return result;
}

BRepClass3d_ImplicitCylinder::NativeFaces BRepClass3d_ImplicitCylinder::BuildImplicitFaces(
  const TopoDS_Face& source, const Settings& s, bool negative, size_t maximumNestingTests,
  int maximumDepth, size_t maximumCells)
{
  NativeFaces result;
  result.Reference = CheckProjectedCorrespondence(source, s, negative, maximumNestingTests,
                                                  maximumDepth, maximumCells);
  if (!result.Reference.MathematicalCorrespondenceVerified)
  { result.Reason = result.Reference.Reason; return result; }
  try
  {
    checkpoint(s);
    const auto data = occ::down_cast<BRep_TImplicitFace>(source.TShape());
    const auto transform = (source.Location() * data->Location()).Transformation();
    const auto inverse = transform.Inverted();
    const double scale = std::abs(transform.ScaleFactor());
    const gp_Ax2 local(s.Frame.Location().Transformed(inverse),
      gp_Dir(gp_Vec(s.Frame.Direction()).Transformed(inverse)),
      gp_Dir(gp_Vec(s.Frame.XDirection()).Transformed(inverse)));
    auto face = BRep_TImplicitCylinderPatch::Restrict(source, local, s.Radius / scale,
      s.HeightMinimum / scale, s.HeightMaximum / scale);
    face.Orientation(TopAbs_FORWARD);
    BRep_Builder builder;
    const auto& side = result.Reference.Source.Boundary.Source.Boundary;
    const double wireBound = side.Source.Chart.Graph.Boundary.Section.GeometricBound
      + side.SurfaceDeviationBound;
    require(std::isfinite(wireBound) && wireBound > 0 && wireBound <= s.DistanceTolerance,
            "implicit_cylinder_trim_tolerance_budget");
    // A new field face has its own trim policy; never mutate the primitive's
    // tolerance or claim a smaller tolerance than the retained contour bound.
    builder.UpdateFace(face, s.DistanceTolerance);
    auto wire = result.Reference.Source.Boundary.Wires.front();
    if (!negative) wire.Reverse(); // Compensate the face orientation below.
    builder.Add(face, wire);
    face.Orientation(negative ? TopAbs_FORWARD : TopAbs_REVERSED);
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> expected;
    for (TopExp_Explorer it(result.Reference.Source.Boundary.Wires.front(), TopAbs_EDGE); it.More(); it.Next())
      expected.Add(it.Current());
    size_t count = 0;
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
    {
      checkpoint(s);
      const int id = expected.FindIndex(it.Current());
      require(id > 0 && it.Current().IsEqual(expected.FindKey(id)), "implicit_cylinder_trim_incidence"); ++count;
    }
    require(count == expected.Size(), "implicit_cylinder_missing_trim_edges");
    result.Faces.push_back(face); result.FacesAttached = true;
    result.Reason = "implicit_cylinder_full_disk_face_attached";
  }
  catch (const Standard_AbortiveTransaction&) { throw; }
  catch (const Standard_Failure& error)
  { result.Faces.clear(); result.Reason = error.what() ? error.what() : "implicit_cylinder_attachment_failed"; }
  return result;
}

BRepClass3d_ImplicitCylinder::NativeShell BRepClass3d_ImplicitCylinder::BuildNativeShell(
  const TopoDS_Face& source, const Settings& s, bool negative, size_t maximumNestingTests,
  int maximumDepth, size_t maximumCells)
{
  NativeShell result;
  result.Source = BuildImplicitFaces(source, s, negative, maximumNestingTests, maximumDepth, maximumCells);
  if (!result.Source.FacesAttached) { result.Reason = result.Source.Reason; return result; }
  try
  {
    const auto& reference = result.Source.Reference;
    const auto& sides = reference.Source.Boundary.Source;
    const int retained = negative == (reference.Source.Direction > 0) ? 0 : 1;
    require(result.Source.Faces.size() == 1 && sides.Faces.size() == 1
              && sides.Boundary.Source.Regions.Faces[retained].size() == 1
              && sides.Boundary.Source.Regions.Faces[1 - retained].empty(),
            "cylinder_shell_full_disk_face_partition");
    std::vector<TopoDS_Face> faces {sides.Boundary.Source.Regions.Faces[retained][0],
                                  sides.Faces[0], result.Source.Faces[0]};
    TopoDS_Shell shell; BRep_Builder builder; builder.MakeShell(shell);
    for (const auto& face : faces) builder.Add(shell, face);
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> edges, vertices;
    TopExp::MapShapes(shell, TopAbs_EDGE, edges); TopExp::MapShapes(shell, TopAbs_VERTEX, vertices);
    struct Use { int Forward = 0, Reversed = 0; std::vector<int> Faces; };
    std::vector<Use> uses(edges.Size());
    std::vector<std::map<int, std::vector<int>>> links(vertices.Size());
    std::vector<std::vector<int>> adjacent(faces.size());
    for (size_t f = 0; f < faces.size(); ++f)
    {
      int wires = 0;
      for (TopoDS_Iterator wi(faces[f]); wi.More(); wi.Next())
      {
        checkpoint(s); ++wires;
        require(wi.Value().ShapeType() == TopAbs_WIRE, "cylinder_shell_nonwire_face_child");
        std::vector<TopoDS_Edge> ordered;
        for (TopoDS_Iterator ei(wi.Value()); ei.More(); ei.Next())
        { require(ei.Value().ShapeType() == TopAbs_EDGE, "cylinder_shell_nonedge_wire_child");
          ordered.push_back(TopoDS::Edge(ei.Value())); }
        // These private builders insert each directed cycle in traversal order.
        // Verify every join rather than assuming this for arbitrary user wires.
        if (wi.Value().Orientation() == TopAbs_REVERSED) std::reverse(ordered.begin(), ordered.end());
        require(!ordered.empty(), "cylinder_shell_empty_wire");
        std::vector<int> starts, ends, startGerms, endGerms;
        for (const auto& edge : ordered)
        {
          checkpoint(s);
          const int id = edges.FindIndex(edge) - 1;
          require(id >= 0, "cylinder_shell_missing_edge");
          const bool forward = edge.Orientation() == TopAbs_FORWARD;
          require(forward || edge.Orientation() == TopAbs_REVERSED, "cylinder_shell_edge_orientation");
          auto& use = uses[size_t(id)]; use.Faces.push_back(int(f));
          if (forward) ++use.Forward; else ++use.Reversed;
          TopoDS_Vertex a, b;
          TopExp::Vertices(TopoDS::Edge(edge.Oriented(TopAbs_FORWARD)), a, b);
          const int av = vertices.FindIndex(a) - 1, bv = vertices.FindIndex(b) - 1;
          require(av >= 0 && bv >= 0, "cylinder_shell_missing_vertex");
          starts.push_back(forward ? av : bv); ends.push_back(forward ? bv : av);
          startGerms.push_back(2 * id + (forward ? 0 : 1));
          endGerms.push_back(2 * id + (forward ? 1 : 0));
        }
        for (size_t i = 0; i < ordered.size(); ++i)
        {
          const size_t next = (i + 1) % ordered.size();
          require(ends[i] == starts[next], "cylinder_shell_unordered_or_open_wire");
          const int a = endGerms[i], b = startGerms[next];
          require(a != b, "cylinder_shell_collapsed_corner");
          links[size_t(ends[i])][a].push_back(b); links[size_t(ends[i])][b].push_back(a);
        }
      }
      require(wires == 1, "cylinder_shell_expected_disk_cell");
    }
    for (size_t id = 0; id < uses.size(); ++id)
    {
      checkpoint(s);
      const auto& use = uses[id];
      require(use.Forward == 1 && use.Reversed == 1 && use.Faces.size() == 2,
              "cylinder_shell_unpaired_edge");
      const int a = use.Faces[0], b = use.Faces[1];
      if (a == b)
        require(!BRep_Tool::IsImplicit(faces[size_t(a)])
                  && BRep_Tool::IsClosed(TopoDS::Edge(edges.FindKey(id + 1)), faces[size_t(a)]),
                "cylinder_shell_unqualified_self_seam");
      else { adjacent[size_t(a)].push_back(b); adjacent[size_t(b)].push_back(a); }
    }
    for (const auto& link : links)
    {
      checkpoint(s);
      require(!link.empty(), "cylinder_shell_unused_vertex");
      for (const auto& germ : link) require(germ.second.size() == 2, "cylinder_shell_noncyclic_vertex_link");
      std::set<int> visited; std::vector<int> todo {link.begin()->first};
      while (!todo.empty())
      {
        const int id = todo.back(); todo.pop_back();
        if (!visited.insert(id).second) continue;
        for (int next : link.at(id)) todo.push_back(next);
      }
      require(visited.size() == link.size(), "cylinder_shell_pinched_vertex");
    }
    std::set<int> visited; std::vector<int> todo {0};
    while (!todo.empty())
    {
      checkpoint(s);
      const int id = todo.back(); todo.pop_back();
      if (!visited.insert(id).second) continue;
      for (int next : adjacent[size_t(id)]) todo.push_back(next);
    }
    require(visited.size() == faces.size(), "cylinder_shell_disconnected_faces");
    const int euler = int(vertices.Size()) - int(edges.Size()) + int(faces.size());
    require(euler == 2, "cylinder_shell_reference_euler_mismatch");
    shell.Closed(true);
    result.Shell = shell; result.Faces = std::move(faces); result.Euler = euler;
    result.Vertices = vertices.Size(); result.Edges = edges.Size();
    result.TopologyValidated = true; result.Reason = "cylinder_native_shell_topology_validated";
  }
  catch (const Standard_AbortiveTransaction&) { throw; }
  catch (const Standard_Failure& error)
  { result.Reason = error.what() ? error.what() : "cylinder_native_shell_failed"; }
  return result;
}

BRepClass3d_ImplicitCylinder::ShellTolerance BRepClass3d_ImplicitCylinder::CheckShellTolerance(
  const TopoDS_Face& face, const Settings& s, bool negative, size_t maximumNestingTests,
  int maximumDepth, size_t maximumCells)
{
  ShellTolerance result;
  result.Source = BuildNativeShell(face, s, negative, maximumNestingTests, maximumDepth, maximumCells);
  if (!result.Source.TopologyValidated) { result.Reason = result.Source.Reason; return result; }
  try
  {
    const auto& side = result.Source.Source.Reference.Source.Boundary.Source.Boundary;
    result.ContourBound = side.Source.Chart.Graph.Boundary.Section.GeometricBound;
    result.SurfaceDeviationBound = side.SurfaceDeviationBound;
    result.FieldWireBound = result.ContourBound + result.SurfaceDeviationBound;
    const auto within = [&](double value) {
      return std::isfinite(value) && value >= 0 && value <= s.DistanceTolerance;
    };
    require(result.ContourBound > 0 && within(result.FieldWireBound), "cylinder_shell_field_wire_bound");
    for (const auto& fieldFace : result.Source.Source.Faces)
      require(BRep_Tool::Tolerance(fieldFace) >= result.FieldWireBound, "cylinder_shell_field_face_underdeclared");
    for (const auto& edge : side.FieldEdges)
      require(BRep_Tool::Tolerance(edge) >= result.SurfaceDeviationBound,
              "cylinder_shell_side_edge_underdeclared");
    for (TopExp_Explorer it(result.Source.Shell, TopAbs_FACE); it.More(); it.Next())
    {
      checkpoint(s);
      const double tolerance = BRep_Tool::Tolerance(TopoDS::Face(it.Current()));
      require(within(tolerance), "cylinder_shell_face_tolerance_budget");
      result.MaximumFaceTolerance = std::max(result.MaximumFaceTolerance, tolerance);
    }
    for (TopExp_Explorer it(result.Source.Shell, TopAbs_EDGE); it.More(); it.Next())
    {
      checkpoint(s);
      const double tolerance = BRep_Tool::Tolerance(TopoDS::Edge(it.Current()));
      require(within(tolerance), "cylinder_shell_edge_tolerance_budget");
      result.MaximumEdgeTolerance = std::max(result.MaximumEdgeTolerance, tolerance);
    }
    for (TopExp_Explorer it(result.Source.Shell, TopAbs_VERTEX); it.More(); it.Next())
    {
      checkpoint(s);
      const double tolerance = BRep_Tool::Tolerance(TopoDS::Vertex(it.Current()));
      require(within(tolerance), "cylinder_shell_vertex_tolerance_budget");
      result.MaximumVertexTolerance = std::max(result.MaximumVertexTolerance, tolerance);
    }
    result.ToleranceVerified = true; result.Reason = "cylinder_shell_tolerance_verified";
  }
  catch (const Standard_AbortiveTransaction&) { throw; }
  catch (const Standard_Failure& error)
  { result.Reason = error.what() ? error.what() : "cylinder_shell_tolerance_unresolved"; }
  return result;
}

BRepClass3d_ImplicitCylinder::ShellGeometry BRepClass3d_ImplicitCylinder::CheckShellGeometry(
  const TopoDS_Face& face, const Settings& s, bool negative, size_t maximumNestingTests,
  int maximumDepth, size_t maximumCells, int maximumGeometryDepth,
  size_t maximumGeometryCells, size_t maximumSegmentTests)
{
  require(maximumGeometryDepth >= 0 && maximumGeometryDepth <= 32
            && maximumGeometryCells > 0 && maximumGeometryCells <= 2000000
            && maximumSegmentTests > 0 && maximumSegmentTests <= 50000000,
          "invalid_cylinder_geometry_budget");
  ShellGeometry result;
  result.Source = CheckShellTolerance(face, s, negative, maximumNestingTests, maximumDepth, maximumCells);
  if (!result.Source.ToleranceVerified) { result.Reason = result.Source.Reason; return result; }
  try
  {
    const auto& reference = result.Source.Source.Source.Reference;
    const auto& axial = reference.Source;
    const auto& sides = axial.Boundary.Source;
    const auto& chart = sides.Boundary.Source.Chart;
    require(sides.Faces.size() == 1 && axial.MonotonicityVerified
              && reference.MathematicalCorrespondenceVerified, "cylinder_geometry_reference_missing");
    const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
    const auto* field = static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
    const double scale = std::abs((face.Location() * data->Location()).Transformation().ScaleFactor());
    const double lipschitz = field->LipschitzBound() / scale;
    const double precision = 8192 * std::numeric_limits<double>::epsilon()
      * (1 + s.Radius + std::abs(s.HeightMinimum) + std::abs(s.HeightMaximum)
         + s.Frame.Location().Distance(gp_Pnt()))
      * (1 + std::abs(s.AngleMinimum) + std::abs(s.AngleMaximum));
    const double reserve = precision * (1 + lipschitz)
      + 8192 * std::numeric_limits<double>::epsilon() * (1 + std::abs(field->IsoValue()));
    require(std::isfinite(reserve) && precision < s.DistanceTolerance, "cylinder_geometry_precision");
    const double classificationReserve = std::max({precision, BRep_Tool::Tolerance(sides.Faces[0]),
      result.Source.MaximumEdgeTolerance, result.Source.MaximumVertexTolerance});
    struct Segment { gp_Pnt2d A, B; };
    std::vector<Segment> segments;
    const auto scaled = [&](const gp_Pnt2d& p) { return gp_Pnt2d(s.Radius * p.X(), p.Y()); };
    const auto append = [&](const gp_Pnt2d& a, const gp_Pnt2d& b) {
      const auto first = scaled(a), last = scaled(b);
      require(std::isfinite(first.X()) && std::isfinite(first.Y()) && std::isfinite(last.X())
                && std::isfinite(last.Y()) && first.Distance(last) > precision,
              "cylinder_geometry_trim_segment");
      segments.push_back({first, last});
    };
    for (TopExp_Explorer it(sides.Faces[0], TopAbs_EDGE); it.More(); it.Next())
    {
      checkpoint(s);
      double first, last;
      const auto pc = BRep_Tool::CurveOnSurface(TopoDS::Edge(it.Current()), sides.Faces[0], first, last);
      require(!pc.IsNull(), "cylinder_geometry_missing_trim");
      const auto bs = occ::down_cast<Geom2d_BSplineCurve>(pc);
      if (!bs.IsNull())
      {
        require(bs->Degree() == 1, "cylinder_geometry_nonlinear_trim");
        for (int i = 2; i <= bs->NbPoles(); ++i) append(bs->Pole(i - 1), bs->Pole(i));
      }
      else
      {
        require(!occ::down_cast<Geom2d_Line>(pc).IsNull(), "cylinder_geometry_unsupported_trim");
        append(pc->Value(first), pc->Value(last));
      }
    }
    require(!segments.empty(), "cylinder_geometry_empty_trim");
    const auto distance = [&](double u, double v) {
      double nearest = std::numeric_limits<double>::infinity();
      for (const auto& segment : segments)
      {
        if ((result.SegmentTests % 256) == 0) checkpoint(s);
        require(result.SegmentTests < maximumSegmentTests, "cylinder_geometry_segment_budget");
        ++result.SegmentTests;
        const double dx = segment.B.X() - segment.A.X(), dy = segment.B.Y() - segment.A.Y();
        const double length = std::hypot(dx, dy), x = dx / length, y = dy / length;
        const double along = std::clamp((u - segment.A.X()) * x + (v - segment.A.Y()) * y, 0., length);
        nearest = std::min(nearest, std::hypot(u - segment.A.X() - along * x,
                                              v - segment.A.Y() - along * y));
      }
      return nearest;
    };
    const auto rootDistance = [&](double u, double v) {
      double low = s.HeightMinimum, high = s.HeightMaximum;
      const double target = s.DistanceTolerance / 64;
      for (int i = 0; i < 128 && high - low > target; ++i)
      {
        checkpoint(s);
        const double middle = low + (high - low) * .5;
        const double value = axial.Direction * BRep_ImplicitTool::Value(face, chart.Surface->Value(u / s.Radius, middle));
        require(std::isfinite(value), "cylinder_geometry_root_nonfinite");
        if (value > reserve) high = middle;
        else if (value < -reserve) low = middle;
        else
        {
          const double error = (std::abs(value) + reserve) / axial.DerivativeLowerBound;
          const double a = std::max(low, middle - error), b = std::min(high, middle + error);
          require(a <= b && b - a < high - low, "cylinder_geometry_root_precision");
          low = a; high = b;
        }
      }
      require(high - low <= target, "cylinder_geometry_root_unresolved");
      return std::max(std::abs(v - low), std::abs(v - high));
    };
    struct Cell { double U0, U1, V0, V1; int Depth; };
    std::vector<Cell> pending {{s.Radius * s.AngleMinimum, s.Radius * s.AngleMaximum,
                               s.HeightMinimum, s.HeightMaximum, 0}};
    double mathematical = 0, native = 0;
    while (!pending.empty())
    {
      checkpoint(s);
      require(result.Cells < maximumGeometryCells, "cylinder_geometry_cell_budget"); ++result.Cells;
      const auto cell = pending.back(); pending.pop_back();
      const double du = cell.U1 - cell.U0, dv = cell.V1 - cell.V0;
      const double u = cell.U0 + du * .5, v = cell.V0 + dv * .5;
      const double radius = std::hypot(du, dv) * .5 + precision;
      const double q = (negative ? -1 : 1) * BRep_ImplicitTool::Value(face, chart.Surface->Value(u / s.Radius, v));
      require(std::isfinite(q) && std::isfinite(radius), "cylinder_geometry_cell_nonfinite");
      const double uncertainty = lipschitz * radius + reserve;
      const bool mathIn = q > uncertainty, mathOut = q < -uncertainty;
      const double nearest = distance(u, v);
      const auto state = BRepClass_FaceClassifier(sides.Faces[0], gp_Pnt2d(u / s.Radius, v), precision).State();
      const bool nativeIn = state == TopAbs_IN && nearest > radius + classificationReserve;
      const bool nativeOut = state == TopAbs_OUT && nearest > radius + classificationReserve;
      const bool nativeWitness = state == TopAbs_IN && nearest > classificationReserve;
      double toNative = 0, toMath = 0;
      if (!mathOut && !nativeIn) toNative = radius + precision + (nativeWitness ? 0 : nearest);
      if (!nativeOut && !mathIn) toMath = radius + precision + (q > reserve ? 0 : rootDistance(u, v));
      require(std::isfinite(toNative) && std::isfinite(toMath), "cylinder_geometry_witness_nonfinite");
      if (toNative <= s.DistanceTolerance && toMath <= s.DistanceTolerance)
      { mathematical = std::max(mathematical, toNative); native = std::max(native, toMath); continue; }
      require(cell.Depth < maximumGeometryDepth, "cylinder_geometry_depth_budget");
      auto a = cell, b = cell; a.Depth = b.Depth = cell.Depth + 1;
      if (du >= dv) { require(u > cell.U0 && u < cell.U1, "cylinder_geometry_split_precision"); a.U1 = b.U0 = u; }
      else { require(v > cell.V0 && v < cell.V1, "cylinder_geometry_split_precision"); a.V1 = b.V0 = v; }
      pending.push_back(b); pending.push_back(a);
    }
    result.MathematicalToNativeBound = mathematical;
    result.NativeToMathematicalBound = native;
    result.BoundaryBound = std::max({mathematical, native, result.Source.FieldWireBound});
    require(result.BoundaryBound <= s.DistanceTolerance, "cylinder_geometry_combined_bound");
    result.BoundaryVerified = true; result.Reason = "cylinder_shell_boundary_geometry_verified";
  }
  catch (const Standard_AbortiveTransaction&) { throw; }
  catch (const Standard_Failure& error)
  { result.Reason = error.what() ? error.what() : "cylinder_shell_geometry_unresolved"; }
  return result;
}
