// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRep_TImplicitCylinderPatch.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <TopoDS_Iterator.hxx>
#include <cmath>
#include <typeinfo>

IMPLEMENT_STANDARD_RTTIEXT(BRep_TImplicitCylinderPatch, BRep_TImplicitFace)

BRep_TImplicitCylinderPatch::BRep_TImplicitCylinderPatch(
  const occ::handle<Geom_ImplicitField>& field, const gp_Pnt& low, const gp_Pnt& high,
  const gp_Ax2& frame, double radius, double minimum, double maximum, double tolerance)
  : BRep_TImplicitFace(field, low, high, tolerance), myFrame(frame), myRadius(radius),
    myHeightMinimum(minimum), myHeightMaximum(maximum)
{
  if (!std::isfinite(radius) || radius <= 0 || !std::isfinite(minimum)
      || !std::isfinite(maximum) || minimum >= maximum)
    throw Standard_DomainError("Invalid finite cylinder restriction");
  for (int axis = 1; axis <= 3; ++axis)
    if (!std::isfinite(frame.Location().Coord(axis))
        || !std::isfinite(frame.XDirection().Coord(axis))
        || !std::isfinite(frame.YDirection().Coord(axis))
        || !std::isfinite(frame.Direction().Coord(axis)))
      throw Standard_DomainError("Non-finite cylinder frame");
  for (double x : {-radius, radius}) for (double y : {-radius, radius})
    for (double z : {minimum, maximum})
    {
      const auto p = frame.Location().Translated(gp_Vec(frame.XDirection()) * x
        + gp_Vec(frame.YDirection()) * y + gp_Vec(frame.Direction()) * z);
      for (int axis = 1; axis <= 3; ++axis)
        if (!std::isfinite(p.Coord(axis)) || p.Coord(axis) < low.Coord(axis)
            || p.Coord(axis) > high.Coord(axis))
          throw Standard_DomainError("Cylinder enclosure leaves the field domain");
    }
}

bool BRep_TImplicitCylinderPatch::Contains(const gp_Pnt& point) const
{
  const gp_Vec delta(myFrame.Location(), point);
  const double height = delta.Dot(gp_Vec(myFrame.Direction()));
  const double radial = std::hypot(delta.Dot(gp_Vec(myFrame.XDirection())),
                                   delta.Dot(gp_Vec(myFrame.YDirection())));
  return std::isfinite(height) && std::isfinite(radial) && radial <= myRadius
         && height >= myHeightMinimum && height <= myHeightMaximum;
}

TopoDS_Face BRep_TImplicitCylinderPatch::Restrict(const TopoDS_Face& source, const gp_Ax2& frame,
  double radius, double minimum, double maximum)
{
  const auto data = occ::down_cast<BRep_TImplicitFace>(source.TShape());
  if (data.IsNull() || typeid(*data) != typeid(BRep_TImplicitFace)
      || !data->ExcludedDomains().empty() || TopoDS_Iterator(source).More())
    throw Standard_NotImplemented("Cylinder restriction requires a childless unrestricted native face");
  occ::handle<BRep_TImplicitCylinderPatch> patch = new BRep_TImplicitCylinderPatch(
    data->Field(), data->Minimum(), data->Maximum(), frame, radius, minimum, maximum, data->Tolerance());
  patch->Location(data->Location());
  TopoDS_Face result;
  result.TShape(patch);
  result.Location(source.Location(), false);
  result.Orientation(source.Orientation());
  return result;
}

occ::handle<TopoDS_TShape> BRep_TImplicitCylinderPatch::EmptyCopy() const
{
  occ::handle<BRep_TImplicitCylinderPatch> copy = new BRep_TImplicitCylinderPatch(
    Field(), Minimum(), Maximum(), myFrame, myRadius, myHeightMinimum, myHeightMaximum, Tolerance());
  copy->Location(Location());
  return copy;
}

void BRep_TImplicitCylinderPatch::DumpJson(Standard_OStream& out, int depth) const
{
  OCCT_DUMP_TRANSIENT_CLASS_BEGIN(out)
  OCCT_DUMP_BASE_CLASS(out, depth, BRep_TImplicitFace)
  OCCT_DUMP_FIELD_VALUES_DUMPED(out, depth, &myFrame)
  OCCT_DUMP_FIELD_VALUE_NUMERICAL(out, myRadius)
  OCCT_DUMP_FIELD_VALUE_NUMERICAL(out, myHeightMinimum)
  OCCT_DUMP_FIELD_VALUE_NUMERICAL(out, myHeightMaximum)
}
