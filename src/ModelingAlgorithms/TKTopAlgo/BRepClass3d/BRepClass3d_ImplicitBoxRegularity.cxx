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
struct Cell
{
  gp_Pnt Low, High;
};

struct Checker
{
  const Geom_ImplicitGyroid&          Field;
  double                              L, H, Reserve, FieldReserve, GradientReserve;
  int                                 MaximumDepth;
  size_t                              MaximumCells, Visited = 0;
  BRepClass3d_ImplicitBox::Regularity Result;

  bool check(const Cell& cell, int dimension, int depth)
  {
    Result.UnresolvedDimension = dimension;
    Result.UnresolvedMinimum   = cell.Low;
    Result.UnresolvedMaximum   = cell.High;
    if (Visited >= MaximumCells)
      return false;
    ++Visited;
    ++Result.Visited[size_t(dimension)];
    const gp_Pnt p((cell.Low.XYZ() + cell.High.XYZ()) * .5);
    auto         gradient = Field.Gradient(p);
    for (int axis = 1; axis <= 3; ++axis)
      if (cell.Low.Coord(axis) == cell.High.Coord(axis))
        gradient.SetCoord(axis, 0);
    const double norm = std::hypot(gradient.X(), gradient.Y(), gradient.Z());
    const gp_Vec size(cell.Low, cell.High);
    const double r         = std::hypot(size.X(), size.Y(), size.Z()) * .5 + 4 * Reserve;
    const double variation = std::min(L * r, (norm + GradientReserve) * r + H * r * r * .5)
                               * (1 + 64 * std::numeric_limits<double>::epsilon())
                             + FieldReserve;
    if (std::abs(Field.Value(p)) > variation)
    {
      ++Result.Excluded[size_t(dimension)];
      return true;
    }
    const double lower = norm - H * r - GradientReserve;
    // Keep a useful quantitative bound, not an arbitrarily tiny positive
    // rounding margin. Projection to each stratum is essential at tangencies.
    if (dimension > 0 && std::isfinite(lower) && lower > GradientReserve && lower >= .5 * norm)
    {
      ++Result.RegularCells[size_t(dimension)];
      auto& minimum = Result.MinimumDerivative[size_t(dimension)];
      minimum       = minimum == 0 ? lower : std::min(minimum, lower);
      return true;
    }
    if (dimension == 0 || depth >= MaximumDepth)
      return false;
    std::vector<Cell> children{cell};
    for (int axis = 1; axis <= 3; ++axis)
    {
      if (cell.Low.Coord(axis) == cell.High.Coord(axis))
        continue;
      const double mid = p.Coord(axis);
      if (!(mid > cell.Low.Coord(axis) && mid < cell.High.Coord(axis)))
        return false;
      const size_t count = children.size();
      for (size_t i = 0; i < count; ++i)
      {
        auto other = children[i];
        children[i].High.SetCoord(axis, mid);
        other.Low.SetCoord(axis, mid);
        children.push_back(other);
      }
    }
    for (const auto& child : children)
      if (!check(child, dimension, depth + 1))
        return false;
    return true;
  }
};
} // namespace

BRepClass3d_ImplicitBox::Regularity BRepClass3d_ImplicitBox::CheckRegularity(
  const TopoDS_Face& face,
  int                maximumDepth,
  size_t             maximumCells)
{
  if (maximumDepth < 1 || maximumDepth > 20 || !maximumCells || maximumCells > 10000000)
    throw Standard_DomainError("Invalid native box regularity budgets");
  const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  if (data.IsNull() || typeid(*data) != typeid(BRep_TImplicitFace)
      || typeid(*data->Field()) != typeid(Geom_ImplicitGyroid) || !face.Location().IsIdentity()
      || !data->Location().IsIdentity()
      || (face.Orientation() != TopAbs_FORWARD && face.Orientation() != TopAbs_REVERSED))
    throw Standard_NotImplemented(
      "Regularity requires an identity-placed native gyroid box domain");
  const auto&  field = *static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  const auto   low = data->Minimum(), high = data->Maximum();
  const double k = 2 * std::acos(-1.) / field.Period(), h = 4 * k * k, l = field.LipschitzBound();
  double       magnitude = 1;
  for (int axis = 1; axis <= 3; ++axis)
  {
    magnitude = std::max({magnitude, std::abs(low.Coord(axis)), std::abs(high.Coord(axis))});
    for (double value : {low.Coord(axis), high.Coord(axis)})
      if (!std::isfinite(k * value) || std::abs(k * value) >= 1.e6)
        throw Standard_DomainError("Regularity trigonometric precision exhausted");
  }
  if (magnitude > 1.e50 || !std::isfinite(h) || h <= 0 || !std::isfinite(l))
    throw Standard_DomainError("Regularity derivative/coordinate precision exhausted");
  const double reserve = 4096 * std::numeric_limits<double>::epsilon() * magnitude;
  Checker checker{field,
                  l,
                  h,
                  reserve,
                  1024 * std::numeric_limits<double>::epsilon() * (1 + std::abs(field.IsoValue()))
                    + l * reserve,
                  h * reserve + 1024 * std::numeric_limits<double>::epsilon() * l,
                  maximumDepth,
                  maximumCells,
                  0,
                  {}};
  checker.Result.HessianBound = h;
  // Dimensions are checked independently, highest first. No field-sign or
  // orientation choice changes the underlying zero-set transversality.
  for (int dimension = 3; dimension >= 0; --dimension)
  {
    for (int mask = 0; mask < 8; ++mask)
    {
      int free = 0;
      for (int axis = 0; axis < 3; ++axis)
        free += (mask >> axis) & 1;
      if (free != dimension)
        continue;
      for (int side = 0; side < 8; ++side)
      {
        if (side & mask)
          continue; // free-coordinate side bits would duplicate strata
        Cell cell{low, high};
        for (int axis = 0; axis < 3; ++axis)
          if (!(mask & (1 << axis)))
          {
            const double value = side & (1 << axis) ? high.Coord(axis + 1) : low.Coord(axis + 1);
            cell.Low.SetCoord(axis + 1, value);
            cell.High.SetCoord(axis + 1, value);
          }
        if (!checker.check(cell, dimension, 0))
          return checker.Result;
      }
    }
    checker.Result.Complete[size_t(dimension)] = true;
  }
  checker.Result.Resolved            = true;
  checker.Result.UnresolvedDimension = -1;
  return checker.Result;
}
