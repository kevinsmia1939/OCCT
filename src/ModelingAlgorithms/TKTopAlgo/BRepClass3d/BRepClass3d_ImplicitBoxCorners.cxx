// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitBox.hxx>
#include <BRep_TImplicitFace.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <typeinfo>

namespace
{
using Ops = BRepClass3d_ImplicitBox;
bool identity(const TopLoc_Location& location)
{
  const auto transform = location.Transformation();
  for (int row = 1; row <= 3; ++row)
    for (int column = 1; column <= 4; ++column)
      if (transform.Value(row, column) != (row == column ? 1. : 0.))
        return false;
  return true;
}
const BRep_TImplicitFace& geometry(const TopoDS_Face& face)
{
  const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  if (data.IsNull() || typeid(*data) != typeid(BRep_TImplicitFace)
      || typeid(*data->Field()) != typeid(Geom_ImplicitGyroid) || !identity(face.Location())
      || !identity(data->Location())
      || (face.Orientation() != TopAbs_FORWARD && face.Orientation() != TopAbs_REVERSED))
    throw Standard_NotImplemented("Corner germs require an identity-placed native gyroid box domain");
  return *data;
}
}

Ops::CornerModel Ops::CheckCornerGerms(const TopoDS_Face& face, double maximumRadius)
{
  if (!std::isfinite(maximumRadius) || maximumRadius <= 0)
    throw Standard_DomainError("Corner neighborhood radius must be finite and positive");
  const auto& data = geometry(face);
  const auto& field = *static_cast<const Geom_ImplicitGyroid*>(data.Field().get());
  const auto low = data.Minimum(), high = data.Maximum();
  const double k = 2 * std::acos(-1.) / field.Period(), h = 4 * k * k;
  const double l = field.LipschitzBound(), eps = std::numeric_limits<double>::epsilon();
  double magnitude = 1, radius = maximumRadius;
  for (int axis = 1; axis <= 3; ++axis)
  {
    magnitude = std::max({magnitude, std::abs(low.Coord(axis)), std::abs(high.Coord(axis))});
    radius = std::min(radius, (high.Coord(axis) - low.Coord(axis)) / 8);
    for (double coordinate : {low.Coord(axis), high.Coord(axis)})
      if (!std::isfinite(k * coordinate) || std::abs(k * coordinate) >= 1.e6)
        throw Standard_DomainError("Corner trigonometric precision exhausted");
  }
  if (magnitude > 1.e50 || !std::isfinite(h) || h <= 0 || !std::isfinite(l)
      || !std::isfinite(radius) || radius <= 0)
    throw Standard_DomainError("Corner derivative/domain precision exhausted");
  CornerModel result;
  result.CoordinateReserve = 4096 * eps * magnitude;
  result.HessianBound = h;
  const double reserve = result.CoordinateReserve;
  const double fieldReserve = 1024 * eps * (1 + std::abs(field.IsoValue())) + l * reserve;
  const double gradientReserve = h * reserve + 1024 * eps * l;
  result.Resolved = true;
  for (int index = 0; index < 8; ++index)
  {
    auto& corner = result.Corners[size_t(index)];
    corner.Anchor = gp_Pnt(index & 1 ? high.X() : low.X(),
                           index & 2 ? high.Y() : low.Y(),
                           index & 4 ? high.Z() : low.Z());
    corner.Minimum = corner.Maximum = corner.Anchor;
    double value = field.Value(corner.Anchor), exactValue = 0;
    gp_Vec exactGradient;
    corner.ExactZero = field.LatticeJet(corner.Anchor, exactValue, exactGradient) && exactValue == 0;
    gp_Vec gradient = corner.ExactZero ? exactGradient : field.Gradient(corner.Anchor);
    const double derivative = std::min({std::abs(gradient.X()), std::abs(gradient.Y()),
                                        std::abs(gradient.Z())});
    double edge = radius;
    if (corner.ExactZero)
    {
      ++result.ExactZeroCount;
      edge = std::min(edge, (derivative - 4 * gradientReserve) / (8 * std::sqrt(3.) * h));
    }
    else if (std::abs(value) > 4 * fieldReserve)
      edge = std::min(edge, (std::abs(value) - 4 * fieldReserve) / (8 * std::sqrt(3.) * l));
    else
    {
      corner.Reason = "corner_zero_identity_unresolved";
      result.Resolved = false;
      continue;
    }
    if (!std::isfinite(edge) || edge <= 16 * reserve)
    {
      corner.Reason = "corner_derivative_or_coordinate_precision_unresolved";
      result.Resolved = false;
      continue;
    }
    for (int axis = 1; axis <= 3; ++axis)
    {
      double end = corner.Anchor.Coord(axis) + (index & (1 << (axis - 1)) ? -edge : edge);
      if (std::abs(static_cast<long double>(end) - corner.Anchor.Coord(axis)) > edge)
        end = std::nextafter(end, corner.Anchor.Coord(axis));
      corner.Minimum.SetCoord(axis, std::min(end, corner.Anchor.Coord(axis)));
      corner.Maximum.SetCoord(axis, std::max(end, corner.Anchor.Coord(axis)));
    }
    const gp_Vec diagonal(corner.Minimum, corner.Maximum);
    const double extent = std::hypot(diagonal.X(), diagonal.Y(), diagonal.Z()) + 4 * reserve;
    const double variation = h * extent + gradientReserve;
    int positive = 0, negative = 0;
    for (int axis = 0; axis < 3; ++axis)
    {
      const double inward = gradient.Coord(axis + 1) * (index & (1 << axis) ? -1 : 1);
      corner.InwardDerivativeLower[size_t(axis)] = inward - variation;
      corner.InwardDerivativeUpper[size_t(axis)] = inward + variation;
      positive += inward - variation > gradientReserve;
      negative += inward + variation < -gradientReserve;
    }
    if (corner.ExactZero && positive + negative == 3)
    {
      corner.Kind = positive == 3 ? CornerKind::PositiveCone
                  : negative == 3 ? CornerKind::NegativeCone : CornerKind::SplitCone;
      result.SplitCount += corner.Kind == CornerKind::SplitCone;
    }
    else if (!corner.ExactZero && std::abs(value) > l * extent + fieldReserve)
      corner.Kind = value > 0 ? CornerKind::Positive : CornerKind::Negative;
    if (corner.Kind == CornerKind::Unresolved)
    {
      corner.Reason = "corner_neighborhood_bound_unresolved";
      result.Resolved = false;
      continue;
    }
    corner.Radius = edge;
    corner.Reason = corner.ExactZero ? "exact_lattice_zero_monotone_corner" : "zero_free_corner";
  }
  return result;
}

std::vector<std::array<gp_Pnt, 2>> Ops::SurfaceCornerExclusions(const TopoDS_Face& face)
{
  const auto& data = geometry(face);
  const auto model = CheckCornerGerms(face, data.Tolerance());
  if (!model.Resolved)
    throw Standard_NotImplemented("Implicit surface corner exclusions unresolved");
  std::vector<std::array<gp_Pnt, 2>> result;
  for (const auto& corner : model.Corners)
    if (corner.ExactZero && (corner.Kind == CornerKind::PositiveCone
                            || corner.Kind == CornerKind::NegativeCone))
      result.push_back({corner.Minimum, corner.Maximum});
  return result;
}

Ops::CornerMembership Ops::ClassifyCornerGerm(const TopoDS_Face& face, int index,
                                            const gp_Pnt& point, Operation operation,
                                            double maximumRadius)
{
  if (index < 0 || index >= 8 || !std::isfinite(point.X()) || !std::isfinite(point.Y())
      || !std::isfinite(point.Z())
      || (operation != Operation::Common && operation != Operation::BoxMinusField))
    throw Standard_DomainError("Invalid corner membership query");
  const auto model = CheckCornerGerms(face, maximumRadius);
  const auto& corner = model.Corners[size_t(index)];
  CornerMembership result;
  if (corner.Kind == CornerKind::Unresolved)
  {
    result.Reason = corner.Reason;
    return result;
  }
  bool anchor = true;
  for (int axis = 1; axis <= 3; ++axis)
  {
    if (point.Coord(axis) < corner.Minimum.Coord(axis)
        || point.Coord(axis) > corner.Maximum.Coord(axis))
    {
      result.Reason = "outside_corner_neighborhood";
      return result;
    }
    anchor &= point.Coord(axis) == corner.Anchor.Coord(axis);
  }
  const bool common = operation == Operation::Common;
  if (anchor && corner.ExactZero)
  {
    result.OnField = result.ClosedRetained = true;
    result.RegularizedRetained = corner.Kind == CornerKind::SplitCone
                              || (common ? corner.Kind == CornerKind::NegativeCone
                                         : corner.Kind == CornerKind::PositiveCone);
  }
  else
  {
    bool positive = corner.Kind == CornerKind::Positive;
    bool negative = corner.Kind == CornerKind::Negative;
    if (corner.ExactZero)
    {
      long double lower = 0, upper = 0, magnitude = 0;
      for (int axis = 0; axis < 3; ++axis)
      {
        const long double displacement =
          (static_cast<long double>(point.Coord(axis + 1)) - corner.Anchor.Coord(axis + 1))
          * (index & (1 << axis) ? -1 : 1);
        const auto a = displacement * corner.InwardDerivativeLower[size_t(axis)];
        const auto b = displacement * corner.InwardDerivativeUpper[size_t(axis)];
        lower += a;
        upper += b;
        magnitude += std::abs(a) + std::abs(b);
      }
      const long double error = 128 * std::numeric_limits<double>::epsilon() * magnitude;
      positive = lower > error && lower > 1024 * std::numeric_limits<double>::min();
      negative = upper < -error && upper < -1024 * std::numeric_limits<double>::min();
    }
    if (!positive && !negative)
    {
      result.Reason = "corner_field_sign_unresolved";
      return result;
    }
    result.ClosedRetained = result.RegularizedRetained = common ? negative : positive;
  }
  result.Resolved = true;
  result.Reason = "local_corner_membership_not_distance_classification";
  return result;
}
