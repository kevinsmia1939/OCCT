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
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS_Wire.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_AbortiveTransaction.hxx>
#include <algorithm>
#include <cmath>
#include <limits>

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

Caps BRepClass3d_ImplicitCylinder::BuildCaps(const TopoDS_Face& face, const Settings& s, bool negative,
                                           size_t maximumNestingTests)
{
  Caps result;
  require(maximumNestingTests > 0 && maximumNestingTests <= 10000000, "invalid_cap_nesting_budget");
  result.MaximumNestingTests = maximumNestingTests;
  result.Field = face; result.Parameters = s; result.Parameters.WorkCheckpoint = nullptr;
  result.Negative = negative;
  result.Source = BRep_ImplicitCylinderSection::BuildCapSkeleton(face, s);
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
