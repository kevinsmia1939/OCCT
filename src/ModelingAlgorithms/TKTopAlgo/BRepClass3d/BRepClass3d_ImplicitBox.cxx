// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitBox.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Result.hxx>
#include <sstream>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <ElSLib.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <typeinfo>

namespace
{
using Map = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

void require(bool condition, const char* message)
{
  if (!condition)
    throw Standard_NotImplemented(message);
}

// Conservative Lipschitz exclusion, not a fixed sample sign test. Edges near
// zeros (including tangency) reject rather than inventing box-edge vertices.
void excludeEdge(const Geom_ImplicitField& field,
                 const gp_Pnt&             a,
                 const gp_Pnt&             b,
                 double                    tolerance,
                 int                       depth = 0)
{
  const gp_Pnt middle((a.XYZ() + b.XYZ()) * .5);
  const double length = a.Distance(b);
  if (std::abs(field.Value(middle)) > field.LipschitzBound() * length * .5 + 1.e-12)
    return;
  require(depth < 24 && length > tolerance / 16,
          "Box edge crosses or approaches the gyroid; open trim topology unsupported");
  excludeEdge(field, a, middle, tolerance, depth + 1);
  excludeEdge(field, middle, b, tolerance, depth + 1);
}

TopoDS_Wire wireFor(BRep_Builder& builder, const TopoDS_Edge& edge, bool forward)
{
  TopoDS_Wire wire;
  builder.MakeWire(wire);
  builder.Add(wire, edge.Oriented(forward ? TopAbs_FORWARD : TopAbs_REVERSED));
  wire.Closed(true);
  return wire;
}
} // namespace

BRepClass3d_ImplicitBox::Boundary BRepClass3d_ImplicitBox::Build(const TopoDS_Face&  fieldFace,
                                                                 const TopoDS_Shape& box,
                                                                 Operation           operation,
                                                                 double              tolerance,
                                                                 int                 maximumCells)
{
  if (!std::isfinite(tolerance) || tolerance <= 0)
    throw Standard_DomainError("Mixed Boolean tolerance must be finite and positive");
  if (operation != Operation::Common && operation != Operation::BoxMinusField)
    throw Standard_DomainError("Unknown mixed Boolean operation");
  const auto data = occ::down_cast<BRep_TImplicitFace>(fieldFace.TShape());
  require(!data.IsNull() && typeid(*data) == typeid(BRep_TImplicitFace)
            && typeid(*data->Field()) == typeid(Geom_ImplicitGyroid)
            && !TopoDS_Iterator(fieldFace).More() && fieldFace.Location().IsIdentity()
            && data->Location().IsIdentity() && fieldFace.Orientation() == TopAbs_FORWARD,
          "Mixed box Boolean requires a childless identity-placed forward native gyroid");
  require(!box.IsNull() && box.ShapeType() == TopAbs_SOLID && box.Orientation() == TopAbs_FORWARD
            && BRepCheck_Analyzer(box).IsValid(),
          "Mixed Boolean operand must be one valid forward conventional box solid");
  Map        vertices, edges, faces, shells;
  const auto snapshot = BRepBuilderAPI_Copy(box).Shape();
  TopExp::MapShapes(snapshot, TopAbs_VERTEX, vertices);
  TopExp::MapShapes(snapshot, TopAbs_EDGE, edges);
  TopExp::MapShapes(snapshot, TopAbs_FACE, faces);
  TopExp::MapShapes(snapshot, TopAbs_SHELL, shells);
  require(vertices.Extent() == 8 && edges.Extent() == 12 && faces.Extent() == 6
            && shells.Extent() == 1,
          "Only six-face axis-aligned boxes are supported");
  gp_Pnt low = BRep_Tool::Pnt(TopoDS::Vertex(vertices(1))), high = low;
  for (const auto& vertex : vertices)
  {
    const gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(vertex));
    for (int j = 1; j <= 3; ++j)
    {
      require(std::isfinite(p.Coord(j)), "Non-finite box coordinate");
      low.SetCoord(j, std::min(low.Coord(j), p.Coord(j)));
      high.SetCoord(j, std::max(high.Coord(j), p.Coord(j)));
    }
  }
  for (int j = 1; j <= 3; ++j)
    require(high.Coord(j) - low.Coord(j) > 4 * tolerance && low.Coord(j) >= data->Minimum().Coord(j)
              && high.Coord(j) <= data->Maximum().Coord(j),
            "Box must fit inside the field domain and exceed four tolerances per axis");
  for (const auto& vertex : vertices)
  {
    const gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(vertex));
    for (int j = 1; j <= 3; ++j)
      require(std::min(std::abs(p.Coord(j) - low.Coord(j)), std::abs(p.Coord(j) - high.Coord(j)))
                < 1.e-7,
              "Box vertices must be axis-aligned corners");
  }
  for (const auto& edge : edges)
  {
    const auto e = TopoDS::Edge(edge);
    require(BRepAdaptor_Curve(e).GetType() == GeomAbs_Line,
            "Only straight-edged boxes are supported");
    TopoDS_Vertex a, b;
    TopExp::Vertices(e, a, b);
    excludeEdge(*data->Field(), BRep_Tool::Pnt(a), BRep_Tool::Pnt(b), tolerance);
  }

  const bool   negative  = operation == Operation::Common;
  const bool   keepOuter = (data->Field()->Value(low) < 0) == negative;
  Boundary     result;
  BRep_Builder builder;
  builder.MakeCompound(result.Shape);
  TopoDS_Face implicit;
  builder.MakeImplicitFace(implicit, data->Field(), low, high, TopLoc_Location(), tolerance);
  if (!negative)
    implicit.Reverse();
  std::array<bool, 6> seen{};
  int                 index = 0;
  for (const auto& shape : faces)
  {
    const auto original = TopoDS::Face(shape);
    require(original.Orientation() == TopAbs_FORWARD || original.Orientation() == TopAbs_REVERSED,
            "Unsupported box face orientation");
    const BRepAdaptor_Surface adaptor(original);
    require(adaptor.GetType() == GeomAbs_Plane, "Only planar box faces are supported");
    const gp_Pln plane = adaptor.Plane();
    gp_Vec       normal(plane.Axis().Direction());
    if (original.Orientation() == TopAbs_REVERSED)
      normal.Reverse();
    int axis = 0;
    for (int j = 1; j <= 3; ++j)
      if (std::abs(normal.Coord(j)) > 1 - 1.e-12)
        axis = j;
    require(axis != 0, "Only axis-aligned box planes are supported");
    const bool upper = normal.Coord(axis) > 0;
    const int  side  = 2 * (axis - 1) + int(upper);
    require(
      !seen[size_t(side)]
        && std::abs(plane.Location().Coord(axis) - (upper ? high.Coord(axis) : low.Coord(axis)))
             < 1.e-7,
      "Box planes must be unique outward boundary planes");
    seen[size_t(side)] = true;
    int wireCount      = 0;
    for (TopoDS_Iterator it(original, false); it.More(); it.Next())
      ++wireCount;
    require(wireCount == 1, "Box faces with holes are unsupported");
    BRep_ImplicitPlaneSection::Settings settings;
    settings.Frame = gp_Ax2(plane.Location(), plane.Axis().Direction(), plane.XAxis().Direction());
    const gp_Pln sectionPlane(settings.Frame);
    settings.UMinimum = settings.VMinimum = 1.e100;
    settings.UMaximum = settings.VMaximum = -1.e100;
    for (TopExp_Explorer it(original, TopAbs_VERTEX); it.More(); it.Next())
    {
      double u, v;
      ElSLib::Parameters(sectionPlane, BRep_Tool::Pnt(TopoDS::Vertex(it.Current())), u, v);
      settings.UMinimum = std::min(settings.UMinimum, u);
      settings.UMaximum = std::max(settings.UMaximum, u);
      settings.VMinimum = std::min(settings.VMinimum, v);
      settings.VMaximum = std::max(settings.VMaximum, v);
    }
    settings.DistanceTolerance = tolerance;
    settings.MaximumCells      = maximumCells;
    const auto section         = BRep_ImplicitPlaneSection::Trace(fieldFace, settings);
    if (!section.Resolved)
    {
      std::ostringstream message;
      message << "Box face " << index << " section unresolved: cells=" << section.UnresolvedCells
              << ", critical=" << section.PotentialCriticalCells
              << ", ambiguous=" << section.AmbiguousCells << ", graph=" << section.GraphDefects
              << ", axis=" << axis << ", coordinate=" << plane.Location().Coord(axis)
              << ", u=" << settings.UMinimum << ":" << settings.UMaximum
              << ", v=" << settings.VMinimum << ":" << settings.VMaximum;
      throw Standard_NotImplemented(message.str().c_str());
    }
    require(section.Curves.size() <= 1, "Multiple/nested contours per box face unsupported");
    result.Sections.push_back(section);
    result.GeometricBound = std::max(result.GeometricBound, section.GeometricBound);
    TopoDS_Face     cap;
    TopLoc_Location location;
    const auto      surface = BRep_Tool::Surface(original, location);
    builder.MakeFace(cap, surface, location, tolerance);
    cap.Orientation(original.Orientation());
    if (keepOuter)
      for (TopoDS_Iterator it(original); it.More(); it.Next())
        builder.Add(cap, it.Value());
    for (const auto& loop : section.Curves)
    {
      require(loop.Closed && loop.Points.size() >= 4,
              "Open box-face contours require edge-crossing trim topology");
      const int                    count = int(loop.Points.size());
      NCollection_Array1<gp_Pnt>   poles(1, count);
      NCollection_Array1<gp_Pnt2d> uv(1, count);
      NCollection_Array1<double>   knots(1, count);
      NCollection_Array1<int>      multiplicities(1, count);
      for (int i = 1; i <= count; ++i)
      {
        poles(i) = loop.Points[size_t(i - 1)];
        double u, v;
        ElSLib::Parameters(sectionPlane, poles(i), u, v);
        require(std::min({u - settings.UMinimum,
                          settings.UMaximum - u,
                          v - settings.VMinimum,
                          settings.VMaximum - v})
                  > 2 * tolerance,
                "Trim seam too close to a box edge for the closed-loop prototype");
        ElSLib::Parameters(plane, poles(i), u, v);
        uv(i).SetCoord(u, v);
        knots(i)          = i - 1;
        multiplicities(i) = (i == 1 || i == count) ? 2 : 1;
      }
      const occ::handle<Geom_BSplineCurve> curve =
        new Geom_BSplineCurve(poles, knots, multiplicities, 1);
      TopoDS_Edge edge;
      builder.MakeEdge(edge, curve, tolerance);
      TopoDS_Vertex vertex;
      builder.MakeVertex(vertex, poles(1), tolerance);
      builder.Add(edge, vertex.Oriented(TopAbs_FORWARD));
      builder.Add(edge, vertex.Oriented(TopAbs_REVERSED));
      builder.Range(edge, 0, count - 1);
      builder.UpdateVertex(vertex, 0, edge, tolerance);
      // Both parameter endpoints belong to the same closed-edge vertex.
      builder.UpdateVertex(TopoDS::Vertex(vertex.Reversed()), count - 1, edge, tolerance);
      edge.Closed(true);
      const occ::handle<Geom2d_BSplineCurve> pcurve =
        new Geom2d_BSplineCurve(uv, knots, multiplicities, 1);
      builder.UpdateEdge(edge, pcurve, cap, tolerance);
      builder.Range(edge, cap, 0, count - 1);
      builder.SameRange(edge, true);
      builder.SameParameter(edge, true);
      // Choose the cap's effective world orientation with the retained region
      // on the left. The implicit use is exactly opposite in the compound.
      const gp_Vec tangent(poles(1), poles(2));
      const gp_Vec desired =
        normal.Crossed(BRep_ImplicitTool::Gradient(fieldFace, poles(1))) * (negative ? 1 : -1);
      const double alignment = tangent.Dot(desired);
      require(std::isfinite(alignment) && std::abs(alignment) > 1.e-14,
              "Cannot orient mixed trim seam");
      const bool effectiveForward = alignment > 0;
      builder.Add(cap, wireFor(builder, edge, effectiveForward));
      builder.Add(implicit, wireFor(builder, edge, !effectiveForward));
      ++result.SeamCount;
    }
    if (keepOuter || !section.Curves.empty())
    {
      BRepCheck_Analyzer check(cap);
      if (!check.IsValid())
      {
        std::ostringstream message;
        message << "Constructed planar cap " << index << " failed conventional validation:";
        Map parts;
        TopExp::MapShapes(cap, parts);
        for (const auto& part : parts)
        {
          const auto status = check.Result(part);
          if (status.IsNull())
            continue;
          for (auto value : status->Status())
            if (value != BRepCheck_NoError)
              message << " type=" << part.ShapeType() << " status=" << value;
          for (status->InitContextIterator(); status->MoreShapeInContext();
               status->NextShapeInContext())
            for (auto value : status->StatusOnShape())
              if (value != BRepCheck_NoError)
                message << " context-type=" << part.ShapeType() << " status=" << value;
        }
        throw Standard_NotImplemented(message.str().c_str());
      }
      result.Faces.push_back(cap);
      result.Sources.push_back(index);
      builder.Add(result.Shape, cap);
    }
    ++index;
  }
  result.Faces.push_back(implicit);
  result.Sources.push_back(-1);
  builder.Add(result.Shape, implicit);
  // Necessary closure/orientation checks for this restricted decomposition.
  // This is deliberately NOT a UV-based whole-solid validity certificate.
  Map boundaryEdges;
  TopExp::MapShapes(result.Shape, TopAbs_EDGE, boundaryEdges);
  std::vector<int> uses(size_t(boundaryEdges.Extent()), 0), balance(uses.size(), 0);
  for (TopExp_Explorer it(result.Shape, TopAbs_EDGE); it.More(); it.Next())
  {
    const size_t edge = size_t(boundaryEdges.FindIndex(it.Current()) - 1);
    ++uses[edge];
    require(it.Current().Orientation() == TopAbs_FORWARD
              || it.Current().Orientation() == TopAbs_REVERSED,
            "Unsupported boundary edge orientation");
    balance[edge] += it.Current().Orientation() == TopAbs_FORWARD ? 1 : -1;
  }
  for (size_t i = 0; i < uses.size(); ++i)
    require(uses[i] == 2 && balance[i] == 0,
            "Hybrid boundary has an unpaired or inconsistently oriented edge");
  return result;
}
