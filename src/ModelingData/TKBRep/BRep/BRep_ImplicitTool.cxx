// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRep_ImplicitTool.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_TImplicitSpherePatch.hxx>
#include <BRep_TImplicitPlanePatch.hxx>
#include <BRep_Tool.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_TypeMismatch.hxx>
#include <algorithm>
#include <cmath>
#include <array>
#include <limits>

namespace
{
const BRep_TImplicitFace& faceData(const TopoDS_Face& f)
{
  if (!BRep_Tool::IsImplicit(f))
    throw Standard_TypeMismatch("BRep_ImplicitTool requires a native implicit face");
  const auto& data = *static_cast<const BRep_TImplicitFace*>(f.TShape().get());
  if ((dynamic_cast<const BRep_TImplicitSpherePatch*>(&data)
       || dynamic_cast<const BRep_TImplicitPlanePatch*>(&data))
      && !data.Location().IsIdentity())
    throw Standard_DomainError(
      "Patch geometry location must remain identity; move its face instance");
  return data;
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
  for (const auto& box : data.ExcludedDomains())
    if (local.X() >= box[0].X() && local.X() <= box[1].X()
        && local.Y() >= box[0].Y() && local.Y() <= box[1].Y()
        && local.Z() >= box[0].Z() && local.Z() <= box[1].Z())
      return false;
  if (const auto* patch = dynamic_cast<const BRep_TImplicitSpherePatch*>(&data))
    return gp_Vec(patch->Plane().Location(), local).Dot(gp_Vec(patch->Plane().Axis().Direction()))
           <= 0;
  if (const auto* patch = dynamic_cast<const BRep_TImplicitPlanePatch*>(&data))
    return gp_Vec(patch->Plane().Location(), local).Dot(gp_Vec(patch->Plane().Axis().Direction()))
           <= 0;
  return true;
}

BRep_ImplicitTool::Projection BRep_ImplicitTool::Project(
  const TopoDS_Face& f, const gp_Pnt& seed, double tolerance, double maximumDistance,
  int maximumIterations)
{
  const auto& data = faceData(f);
  checkPoint(seed);
  if (!std::isfinite(tolerance) || tolerance <= 0 || !std::isfinite(maximumDistance)
      || maximumDistance <= 0 || maximumIterations < 1 || maximumIterations > 64)
    throw Standard_DomainError("Invalid native projection tolerance/distance/iteration limit");
  Projection result;
  result.Point = seed;
  const auto bounds = Bounds(f);
  double xmin, ymin, zmin, xmax, ymax, zmax;
  bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  const double magnitude = std::max({1., std::abs(xmin), std::abs(ymin), std::abs(zmin),
    std::abs(xmax), std::abs(ymax), std::abs(zmax), std::abs(seed.X()),
    std::abs(seed.Y()), std::abs(seed.Z())});
  const double reserve = 512 * std::numeric_limits<double>::epsilon() * magnitude;
  if (!std::isfinite(reserve) || reserve >= tolerance / 65536)
  {
    result.Reason = "coordinate_precision";
    return result;
  }
  const double step = std::max(16 * reserve, std::min(tolerance / 32,
    std::cbrt(std::numeric_limits<double>::epsilon()) * magnitude));
  using Hessian = std::array<std::array<double, 3>, 3>;
  const auto hessian = [&](const gp_Pnt& p)
  {
    Hessian h{};
    for (int axis = 0; axis < 3; ++axis)
    {
      gp_Vec offset;
      offset.SetCoord(axis + 1, step);
      const auto column = (Gradient(f, p.Translated(offset))
                           - Gradient(f, p.Translated(-offset))) / (2 * step);
      for (int row = 0; row < 3; ++row)
        h[size_t(row)][size_t(axis)] = column.Coord(row + 1);
    }
    for (int i = 0; i < 3; ++i)
      for (int j = i + 1; j < 3; ++j)
        h[size_t(i)][size_t(j)] = h[size_t(j)][size_t(i)] =
          (h[size_t(i)][size_t(j)] + h[size_t(j)][size_t(i)]) * .5;
    return h;
  };
  const auto product = [](const Hessian& h, const gp_Vec& v)
  {
    gp_Vec answer;
    for (int row = 0; row < 3; ++row)
      for (int col = 0; col < 3; ++col)
        answer.SetCoord(row + 1, answer.Coord(row + 1) + h[size_t(row)][size_t(col)] * v.Coord(col + 1));
    return answer;
  };
  const auto admissibleBracket = [&](const gp_Pnt& a, const gp_Pnt& b)
  {
    if (!IsInDomain(f, a) || !IsInDomain(f, b))
      return false;
    // Endpoint checks suffice for the convex box and patch halfspaces. Closed
    // exclusion boxes also require a segment test, so a thin removed region
    // cannot hide between admissible endpoints.
    const auto inverse = placement(f).Inverted();
    const auto localA = a.Transformed(inverse), localB = b.Transformed(inverse);
    for (const auto& box : data.ExcludedDomains())
    {
      double low = 0, high = 1;
      for (int axis = 1; axis <= 3 && low <= high; ++axis)
      {
        const double d = localB.Coord(axis) - localA.Coord(axis);
        if (d == 0)
        {
          if (localA.Coord(axis) < box[0].Coord(axis) || localA.Coord(axis) > box[1].Coord(axis))
            high = -1;
        }
        else
        {
          double first = (box[0].Coord(axis) - localA.Coord(axis)) / d;
          double last = (box[1].Coord(axis) - localA.Coord(axis)) / d;
          if (first > last) std::swap(first, last);
          low = std::max(low, first); high = std::min(high, last);
        }
      }
      if (low <= high) return false;
    }
    return true;
  };
  double lambda = 0;
  for (int iteration = 0; iteration <= maximumIterations; ++iteration)
  {
    const auto g = Gradient(f, result.Point);
    const double norm = checked(g.Magnitude()), value = Value(f, result.Point);
    result.Iterations = iteration;
    result.Residual = std::abs(value);
    result.Distance = seed.Distance(result.Point);
    const gp_Vec stationary = gp_Vec(seed, result.Point) + g * lambda;
    result.Stationarity = stationary.Magnitude();
    if (norm <= std::numeric_limits<double>::min())
    {
      result.Reason = "critical_point";
      return result;
    }
    const auto h = hessian(result.Point);
    if (result.Residual / norm <= tolerance / 4096 && result.Stationarity <= tolerance / 64)
    {
      const gp_Vec normal = g / norm;
      const gp_Vec reference = std::abs(normal.X()) < .8 ? gp_Vec(1, 0, 0) : gp_Vec(0, 1, 0);
      const gp_Vec u = reference.Crossed(normal).Normalized(), v = normal.Crossed(u);
      const double a = u.Dot(u + product(h, u) * lambda);
      const double b = u.Dot(v + product(h, v) * lambda);
      const double c = v.Dot(v + product(h, v) * lambda);
      if (!std::isfinite(a + b + c) || (a + c - std::hypot(a - c, 2 * b)) * .5 <= 1.e-8)
      {
        result.Reason = "non_minimum";
        return result;
      }
      gp_Pnt first = result.Point.Translated(normal * (-tolerance / 1024));
      gp_Pnt last = result.Point.Translated(normal * (tolerance / 1024));
      double firstValue = Value(f, first), lastValue = Value(f, last);
      const double fieldReserve = 64 * std::numeric_limits<double>::epsilon()
                                  * (1 + result.Residual + norm * magnitude);
      if (!(firstValue < -fieldReserve && lastValue > fieldReserve))
      {
        result.Reason = "unresolved_bracket";
        return result;
      }
      if (!admissibleBracket(first, last))
      {
        result.Reason = "outside_domain";
        return result;
      }
      for (int refine = 0; refine < 64 && first.Distance(last) > tolerance / 65536; ++refine)
      {
        const gp_Pnt mid((first.XYZ() + last.XYZ()) * .5);
        const double midValue = Value(f, mid);
        if (std::abs(midValue) <= fieldReserve) break;
        if (midValue < 0) first = mid;
        else last = mid;
      }
      result.First = first; result.Last = last;
      result.Point = gp_Pnt((first.XYZ() + last.XYZ()) * .5);
      result.BracketWidth = first.Distance(last);
      result.Residual = std::abs(Value(f, result.Point));
      result.Stationarity = (gp_Vec(seed, result.Point) + Gradient(f, result.Point) * lambda).Magnitude();
      result.Distance = seed.Distance(result.Point);
      if (result.Distance > maximumDistance || result.Stationarity > tolerance / 32
          || result.BracketWidth > tolerance / 256 || !IsInDomain(f, result.Point))
      {
        result.Reason = "unresolved_refinement";
        return result;
      }
      result.Converged = true;
      result.Reason = "converged";
      return result;
    }
    if (iteration == maximumIterations) break;
    std::array<std::array<double, 5>, 4> matrix{};
    for (int row = 0; row < 3; ++row)
    {
      for (int col = 0; col < 3; ++col)
        matrix[size_t(row)][size_t(col)] = (row == col ? 1. : 0.) + lambda * h[size_t(row)][size_t(col)];
      matrix[size_t(row)][3] = g.Coord(row + 1);
      matrix[size_t(row)][4] = -stationary.Coord(row + 1);
      matrix[3][size_t(row)] = g.Coord(row + 1) / norm;
    }
    matrix[3][4] = -value / norm;
    for (auto& row : matrix)
    {
      double scale = 0;
      for (int col = 0; col < 4; ++col) scale = std::max(scale, std::abs(row[size_t(col)]));
      if (!std::isfinite(scale) || scale == 0) { result.Reason = "singular_system"; return result; }
      for (auto& entry : row) entry /= scale;
    }
    for (int col = 0; col < 4; ++col)
    {
      int pivot = col;
      for (int row = col + 1; row < 4; ++row)
        if (std::abs(matrix[size_t(row)][size_t(col)]) > std::abs(matrix[size_t(pivot)][size_t(col)])) pivot = row;
      if (std::abs(matrix[size_t(pivot)][size_t(col)]) < 1.e-14)
      { result.Reason = "singular_system"; return result; }
      std::swap(matrix[size_t(pivot)], matrix[size_t(col)]);
      const double divisor = matrix[size_t(col)][size_t(col)];
      for (int j = col; j < 5; ++j) matrix[size_t(col)][size_t(j)] /= divisor;
      for (int row = 0; row < 4; ++row)
        if (row != col)
        {
          const double factor = matrix[size_t(row)][size_t(col)];
          for (int j = col; j < 5; ++j) matrix[size_t(row)][size_t(j)] -= factor * matrix[size_t(col)][size_t(j)];
        }
    }
    const gp_Vec delta(matrix[0][4], matrix[1][4], matrix[2][4]);
    const double merit = result.Residual / norm + result.Stationarity;
    bool accepted = false;
    for (double fraction = 1; fraction >= 1. / 1024; fraction *= .5)
    {
      const auto trial = result.Point.Translated(delta * fraction);
      const double trialLambda = lambda + matrix[3][4] * fraction;
      const auto trialGradient = Gradient(f, trial);
      const double trialNorm = checked(trialGradient.Magnitude());
      if (trialNorm <= std::numeric_limits<double>::min() || seed.Distance(trial) > maximumDistance) continue;
      const double trialMerit = std::abs(Value(f, trial)) / trialNorm
                               + (gp_Vec(seed, trial) + trialGradient * trialLambda).Magnitude();
      if (std::isfinite(trialMerit) && trialMerit < merit)
      { result.Point = trial; lambda = trialLambda; accepted = true; break; }
    }
    if (!accepted) { result.Reason = "step_limit"; return result; }
  }
  return result;
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
    if (!IsInDomain(f, p))
      return;
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
