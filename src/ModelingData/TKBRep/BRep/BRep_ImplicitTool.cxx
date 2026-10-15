// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRep_ImplicitTool.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_Tool.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_TypeMismatch.hxx>
#include <algorithm>
#include <cmath>

namespace
{
const BRep_TImplicitFace& faceData(const TopoDS_Face& f)
{
  if (!BRep_Tool::IsImplicit(f))
    throw Standard_TypeMismatch("BRep_ImplicitTool requires a native implicit face");
  return *static_cast<const BRep_TImplicitFace*>(f.TShape().get());
}

double checked(double v)
{
  if (!std::isfinite(v))
    throw Standard_DomainError("Non-finite implicit face query");
  return v;
}

void checkPoint(const gp_Pnt& p)
{
  checked(p.X());
  checked(p.Y());
  checked(p.Z());
}

gp_Trsf placement(const TopoDS_Face& f)
{
  return (f.Location() * faceData(f).Location()).Transformation();
}
} // namespace

double BRep_ImplicitTool::Value(const TopoDS_Face& f, const gp_Pnt& p)
{
  checkPoint(p);
  gp_Pnt local = p.Transformed(placement(f).Inverted());
  checkPoint(local);
  return checked(faceData(f).Field()->Value(local));
}

gp_Vec BRep_ImplicitTool::Gradient(const TopoDS_Face& f, const gp_Pnt& p)
{
  checkPoint(p);
  gp_Trsf t     = placement(f);
  gp_Pnt  local = p.Transformed(t.Inverted());
  checkPoint(local);
  gp_Vec g = faceData(f).Field()->Gradient(local);
  // Inverse transpose for the uniform-scale orthogonal linear part.
  g.Transform(t);
  g /= t.ScaleFactor() * t.ScaleFactor();
  checked(g.X());
  checked(g.Y());
  checked(g.Z());
  return g;
}

gp_Vec BRep_ImplicitTool::Normal(const TopoDS_Face& f, const gp_Pnt& p)
{
  gp_Vec n      = Gradient(f, p);
  double length = checked(std::hypot(n.X(), n.Y(), n.Z()));
  if (length == 0 || (f.Orientation() != TopAbs_FORWARD && f.Orientation() != TopAbs_REVERSED))
    throw Standard_DomainError("Implicit face normal undefined at critical point/orientation");
  n /= length;
  if (f.Orientation() == TopAbs_REVERSED)
    n.Reverse();
  return n;
}

Bnd_Box BRep_ImplicitTool::Bounds(const TopoDS_Face& f)
{
  const auto& data = faceData(f);
  gp_Trsf     t    = placement(f);
  Bnd_Box     bounds;
  for (int x = 0; x < 2; ++x)
    for (int y = 0; y < 2; ++y)
      for (int z = 0; z < 2; ++z)
      {
        gp_Pnt p(x ? data.Maximum().X() : data.Minimum().X(),
                 y ? data.Maximum().Y() : data.Minimum().Y(),
                 z ? data.Maximum().Z() : data.Minimum().Z());
        p.Transform(t);
        checkPoint(p);
        bounds.Add(p);
      }
  return bounds;
}

bool BRep_ImplicitTool::IsInDomain(const TopoDS_Face& f, const gp_Pnt& p)
{
  checkPoint(p);
  const auto& data  = faceData(f);
  gp_Pnt      local = p.Transformed(placement(f).Inverted());
  checkPoint(local);
  for (int axis = 1; axis <= 3; ++axis)
    if (local.Coord(axis) < data.Minimum().Coord(axis)
        || local.Coord(axis) > data.Maximum().Coord(axis))
      return false;
  return true;
}

std::vector<BRep_ImplicitTool::RayHit> BRep_ImplicitTool::IntersectRay(const TopoDS_Face& f,
                                                                       const gp_Pnt&      origin,
                                                                       const gp_Vec&      direction,
                                                                       double             minimum,
                                                                       double             maximum,
                                                                       int                samples,
                                                                       double             tolerance)
{
  const auto& data = faceData(f);
  checkPoint(origin);
  const double length = checked(std::hypot(direction.X(), direction.Y(), direction.Z()));
  if (tolerance == 0)
    tolerance = BRep_Tool::Tolerance(f);
  if (length <= 0 || !std::isfinite(minimum) || !std::isfinite(maximum) || minimum > maximum
      || !std::isfinite(maximum - minimum) || samples < 1 || samples > 65536
      || !std::isfinite(tolerance) || tolerance <= 0)
    throw Standard_DomainError(
      "Invalid implicit ray interval/direction/samples/geometric tolerance");
  gp_Vec  unit           = direction / length;
  gp_Trsf inverse        = placement(f).Inverted();
  gp_Pnt  localOrigin    = origin.Transformed(inverse);
  gp_Vec  localDirection = unit.Transformed(inverse);
  checkPoint(localOrigin);
  checked(localDirection.X());
  checked(localDirection.Y());
  checked(localDirection.Z());
  for (int axis = 1; axis <= 3; ++axis)
  {
    double d = localDirection.Coord(axis);
    if (d == 0)
    {
      if (localOrigin.Coord(axis) < data.Minimum().Coord(axis)
          || localOrigin.Coord(axis) > data.Maximum().Coord(axis))
        return {};
      continue;
    }
    double a = (data.Minimum().Coord(axis) - localOrigin.Coord(axis)) / d;
    double b = (data.Maximum().Coord(axis) - localOrigin.Coord(axis)) / d;
    if (a > b)
      std::swap(a, b);
    minimum = std::max(minimum, a);
    maximum = std::min(maximum, b);
    if (minimum > maximum)
      return {};
  }
  std::vector<RayHit> hits;
  auto                point  = [&](double t) { return origin.Translated(unit * t); };
  auto                append = [&](double t, double width) {
    if (!hits.empty() && t - hits.back().Distance <= tolerance)
      return;
    gp_Pnt p = point(t);
    hits.push_back({t, p, std::abs(Value(f, p)), width});
  };
  double previous = minimum, previousValue = Value(f, point(previous));
  if (previousValue == 0)
    append(previous, 0);
  if (minimum == maximum)
    return hits;
  for (int i = 1; i <= samples; ++i)
  {
    double current = i == samples ? maximum : minimum + (maximum - minimum) * (double(i) / samples);
    double value   = Value(f, point(current));
    if (value == 0 && previousValue == 0)
      throw Standard_DomainError(
        "Implicit ray adjacent zero samples: non-isolated/ambiguous roots");
    if (value == 0)
      append(current, 0);
    else if (previousValue != 0 && std::signbit(value) != std::signbit(previousValue))
    {
      double low = previous, high = current, lowValue = previousValue;
      for (int k = 0; k < 100 && high - low > tolerance; ++k)
      {
        double mid = low + (high - low) * 0.5;
        if (mid == low || mid == high)
          break;
        double v = Value(f, point(mid));
        if (v == 0)
          low = high = mid;
        else if (std::signbit(v) != std::signbit(lowValue))
          high = mid;
        else
        {
          low      = mid;
          lowValue = v;
        }
      }
      if (high - low > tolerance)
        throw Standard_DomainError("Implicit ray cannot reach requested geometric tolerance");
      append(low + (high - low) * 0.5, high - low);
    }
    previous      = current;
    previousValue = value;
  }
  return hits;
}
