// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRep_ImplicitGyroidOps.hxx>
#include <BRep_TImplicitPlanePatch.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_Builder.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <algorithm>
#include <cmath>
#include <typeinfo>

IMPLEMENT_STANDARD_RTTIEXT(BRep_TImplicitPlanePatch, BRep_TImplicitFace)

BRep_TImplicitPlanePatch::BRep_TImplicitPlanePatch(const occ::handle<Geom_ImplicitField>& field,
                                                   const gp_Pnt&                          low,
                                                   const gp_Pnt&                          high,
                                                   const gp_Pln&                          plane,
                                                   double                                 tolerance)
    : BRep_TImplicitFace(field, low, high, tolerance),
      myPlane(plane)
{
  for (int i = 1; i <= 3; ++i)
    if (!std::isfinite(plane.Location().Coord(i))
        || !std::isfinite(plane.Axis().Direction().Coord(i)))
      throw Standard_DomainError("Non-finite native trim plane");
}

occ::handle<TopoDS_TShape> BRep_TImplicitPlanePatch::EmptyCopy() const
{
  occ::handle<BRep_TImplicitPlanePatch> copy =
    new BRep_TImplicitPlanePatch(Field(), Minimum(), Maximum(), myPlane, Tolerance());
  copy->Location(Location());
  return copy;
}

void BRep_TImplicitPlanePatch::DumpJson(Standard_OStream& out, int depth) const {
  OCCT_DUMP_TRANSIENT_CLASS_BEGIN(out) OCCT_DUMP_BASE_CLASS(out, depth, BRep_TImplicitFace)
    OCCT_DUMP_FIELD_VALUES_DUMPED(out, depth, &myPlane)}

BRep_ImplicitGyroidOps::Split BRep_ImplicitGyroidOps::SplitZ(const TopoDS_Face& face,
                                                             double             z,
                                                             double             tolerance,
                                                             int                maximumCells)
{
  const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  if (data.IsNull() || !data->ExcludedDomains().empty() || typeid(*data) != typeid(BRep_TImplicitFace)
      || typeid(*data->Field()) != typeid(Geom_ImplicitGyroid) || TopoDS_Iterator(face).More()
      || face.Orientation() != TopAbs_FORWARD || !face.Location().IsIdentity()
      || !data->Location().IsIdentity())
    throw Standard_NotImplemented(
      "Gyroid SplitZ requires an untrimmed identity-placed forward native gyroid");
  if (!std::isfinite(z) || z <= data->Minimum().Z() || z >= data->Maximum().Z())
    throw Standard_DomainError("Split plane must lie strictly inside the field domain");
  if (!std::isfinite(data->Tolerance()) || data->Tolerance() <= 0)
    throw Standard_DomainError("Invalid source face tolerance");
  BRep_ImplicitPlaneSection::Settings settings;
  settings.Frame             = gp_Ax2(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
  settings.UMinimum          = data->Minimum().X();
  settings.UMaximum          = data->Maximum().X();
  settings.VMinimum          = data->Minimum().Y();
  settings.VMaximum          = data->Maximum().Y();
  settings.DistanceTolerance = tolerance;
  settings.MaximumCells      = maximumCells;
  Split result;
  result.Section = BRep_ImplicitPlaneSection::Trace(face, settings);
  if (!result.Section.Resolved)
    throw Standard_DomainError("Gyroid section is unresolved; no trim topology created");
  for (const auto& curve : result.Section.Curves)
    if (!curve.Closed)
      throw Standard_NotImplemented("Open gyroid section requires domain-side boundary topology");
  const double faceTolerance = std::max(tolerance, data->Tolerance());
  result.Negative.TShape(new BRep_TImplicitPlanePatch(data->Field(),
                                                      data->Minimum(),
                                                      data->Maximum(),
                                                      gp_Pln(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1)),
                                                      faceTolerance));
  result.Positive.TShape(new BRep_TImplicitPlanePatch(data->Field(),
                                                      data->Minimum(),
                                                      data->Maximum(),
                                                      gp_Pln(gp_Pnt(0, 0, z), gp_Dir(0, 0, -1)),
                                                      faceTolerance));
  result.Negative.Orientation(TopAbs_FORWARD);
  result.Positive.Orientation(TopAbs_FORWARD);
  BRep_Builder builder;
  for (const auto& loop : result.Section.Curves)
  {
    const int                  count = int(loop.Points.size());
    NCollection_Array1<gp_Pnt> poles(1, count);
    NCollection_Array1<double> knots(1, count);
    NCollection_Array1<int>    multiplicities(1, count);
    for (int i = 1; i <= count; ++i)
    {
      poles(i)          = loop.Points[size_t(i - 1)];
      knots(i)          = i - 1;
      multiplicities(i) = (i == 1 || i == count) ? 2 : 1;
    }
    const occ::handle<Geom_BSplineCurve> curve =
      new Geom_BSplineCurve(poles, knots, multiplicities, 1);
    TopoDS_Edge edge;
    builder.MakeEdge(edge, curve, faceTolerance);
    builder.Range(edge, 0, count - 1);
    TopoDS_Vertex seam;
    builder.MakeVertex(seam, loop.Points.front(), faceTolerance);
    builder.Add(edge, seam.Oriented(TopAbs_FORWARD));
    builder.Add(edge, seam.Oriented(TopAbs_REVERSED));
    builder.UpdateVertex(seam, 0, edge, faceTolerance);
    edge.Closed(true);
    const gp_Vec tangent(loop.Points[0], loop.Points[1]);
    const gp_Vec desired =
      BRep_ImplicitTool::Gradient(face, loop.Points[0]).Crossed(gp_Vec(0, 0, 1));
    const double alignment = tangent.Dot(desired);
    if (!std::isfinite(alignment) || alignment == 0)
      throw Standard_DomainError("Cannot orient gyroid trim seam");
    TopoDS_Wire wire;
    builder.MakeWire(wire);
    builder.Add(wire, edge.Oriented(alignment > 0 ? TopAbs_FORWARD : TopAbs_REVERSED));
    wire.Closed(true);
    builder.Add(result.Negative, wire.Oriented(TopAbs_FORWARD));
    builder.Add(result.Positive, wire.Oriented(TopAbs_REVERSED));
  }
  return result;
}
