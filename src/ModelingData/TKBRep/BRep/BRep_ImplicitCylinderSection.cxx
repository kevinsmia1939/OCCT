// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRep_ImplicitCylinderSection.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_TImplicitFace.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <TopoDS_Iterator.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <typeinfo>

namespace
{
using Section = BRep_ImplicitCylinderSection;
void require(bool test, const char* message)
{
  if (!test) throw Standard_DomainError(message);
}
void checkpoint(const Section::Settings& s)
{
  if (s.WorkCheckpoint) s.WorkCheckpoint();
}
bool cells(int n) { return n >= 4 && n <= 512 && (n & (n - 1)) == 0; }
bool opposite(double a, double b) { return std::signbit(a) != std::signbit(b); }
struct Chart
{
  const Section::Settings& s;
  gp_Pnt point(double u, double v) const
  {
    const double angle = u / s.Radius;
    return s.Frame.Location().Translated(
      gp_Vec(s.Frame.XDirection()) * (s.Radius * std::cos(angle))
      + gp_Vec(s.Frame.YDirection()) * (s.Radius * std::sin(angle))
      + gp_Vec(s.Frame.Direction()) * v);
  }
  gp_Vec tangent(double u) const
  {
    const double angle = u / s.Radius;
    return gp_Vec(s.Frame.XDirection()) * -std::sin(angle)
           + gp_Vec(s.Frame.YDirection()) * std::cos(angle);
  }
};
struct Edge { int state = 0; size_t vertex = 0; }; // empty, unique root, unresolved

Section::Result pass(const TopoDS_Face& face, const Section::Settings& s,
                     int n, double lipschitz, double hessian, double reserve, double precision)
{
  Section::Result result;
  result.Resolution = n;
  const Chart chart {s};
  const double u0 = s.Radius * s.AngleMinimum;
  const double du = s.Radius * (s.AngleMaximum - s.AngleMinimum) / n;
  const double dv = (s.HeightMaximum - s.HeightMinimum) / n;
  const double diagonal = std::hypot(du, dv);
  const size_t width = size_t(n) + 1;
  std::vector<double> values(width * width);
  for (int j = 0; j <= n; ++j)
    for (int i = 0; i <= n; ++i)
    {
      if ((size_t(j) * width + i) % 256 == 0) checkpoint(s);
      values[size_t(j) * width + i] = BRep_ImplicitTool::Value(
        face, chart.point(u0 + du * i, s.HeightMinimum + dv * j));
    }
  std::vector<gp_Pnt> points;
  std::vector<bool> boundary;
  std::vector<Edge> horizontal(size_t(n) * width), vertical(size_t(n) * width);
  const auto edge = [&](int i, int j, bool alongU) {
    if ((size_t(j) * width + i) % 256 == 0) checkpoint(s);
    const double u = u0 + du * i, v = s.HeightMinimum + dv * j;
    const double length = alongU ? du : dv;
    const size_t a = size_t(j) * width + i, b = a + (alongU ? 1 : width);
    const double middleU = u + (alongU ? du * .5 : 0);
    const double middleV = v + (alongU ? 0 : dv * .5);
    const auto middle = chart.point(middleU, middleV);
    const double derivative = BRep_ImplicitTool::Gradient(face, middle).Dot(
      alongU ? chart.tangent(middleU) : gp_Vec(s.Frame.Direction()));
    const bool monotone = std::abs(derivative) > hessian * length * .5 + reserve;
    Edge e;
    if (std::abs(values[a]) <= reserve || std::abs(values[b]) <= reserve)
    { e.state = 2; return e; }
    if (!opposite(values[a], values[b]))
    {
      if (!monotone && std::min(std::abs(values[a]), std::abs(values[b]))
                           <= hessian * length * length / 8 + reserve)
        e.state = 2;
      return e;
    }
    if (!monotone) { e.state = 2; return e; }
    double low = 0, high = 1, first = values[a];
    for (int iteration = 0; iteration < 100 && (high - low) * length > s.DistanceTolerance / 32;
         ++iteration)
    {
      const double t = low + (high - low) * .5;
      const double f = BRep_ImplicitTool::Value(face, chart.point(
        u + (alongU ? du * t : 0), v + (alongU ? 0 : dv * t)));
      if (f == 0) { low = high = t; break; }
      if (opposite(first, f)) high = t;
      else { low = t; first = f; }
    }
    const double bracket = (high - low) * length;
    if (bracket > s.DistanceTolerance / 32) { e.state = 2; return e; }
    const double t = low + (high - low) * .5;
    const auto p = chart.point(u + (alongU ? du * t : 0), v + (alongU ? 0 : dv * t));
    e.state = 1; e.vertex = points.size(); points.push_back(p);
    boundary.push_back(alongU ? j == 0 || j == n : i == 0 || i == n);
    result.MaximumResidual = std::max(result.MaximumResidual, std::abs(BRep_ImplicitTool::Value(face, p)));
    result.MaximumRootBracket = std::max(result.MaximumRootBracket, bracket);
    return e;
  };
  for (int j = 0; j <= n; ++j)
    for (int i = 0; i < n; ++i) horizontal[size_t(j) * n + i] = edge(i, j, true);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i <= n; ++i) vertical[size_t(j) * width + i] = edge(i, j, false);
  std::vector<std::pair<size_t, size_t>> segments;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
    {
      if ((size_t(j) * n + i) % 256 == 0) checkpoint(s);
      const double u = u0 + du * (i + .5);
      const auto center = chart.point(u, s.HeightMinimum + dv * (j + .5));
      const double f = BRep_ImplicitTool::Value(face, center);
      if (std::abs(f) > lipschitz * diagonal * .5 + reserve) continue;
      const auto g = BRep_ImplicitTool::Gradient(face, center);
      const double uncertainty = hessian * diagonal * .5 + reserve;
      const bool regular = std::abs(g.Dot(chart.tangent(u))) > uncertainty
                        || std::abs(g.Dot(gp_Vec(s.Frame.Direction()))) > uncertainty;
      std::array<Edge, 4> edges {{horizontal[size_t(j) * n + i], horizontal[size_t(j + 1) * n + i],
        vertical[size_t(j) * width + i], vertical[size_t(j) * width + i + 1]}};
      std::vector<size_t> roots;
      bool unknown = false;
      for (const auto& e : edges) { unknown |= e.state == 2; if (e.state == 1) roots.push_back(e.vertex); }
      if (regular && !unknown && roots.empty()) continue;
      if (regular && !unknown && roots.size() == 2 && diagonal + precision <= s.DistanceTolerance)
      {
        segments.emplace_back(roots[0], roots[1]);
        result.GeometricBound = std::max(result.GeometricBound, diagonal + precision);
      }
      else
      {
        ++result.UnresolvedCells;
        if (!regular) ++result.PotentialCriticalCells;
        if (unknown || roots.size() != 2) ++result.AmbiguousCells;
      }
    }
  std::vector<std::vector<size_t>> adjacent(points.size());
  for (size_t index = 0; index < segments.size(); ++index)
  { adjacent[segments[index].first].push_back(index); adjacent[segments[index].second].push_back(index); }
  for (size_t index = 0; index < adjacent.size(); ++index)
    if (!adjacent[index].empty() && adjacent[index].size() != 2
        && !(boundary[index] && adjacent[index].size() == 1)) ++result.GraphDefects;
  std::vector<bool> used(segments.size(), false);
  const auto walk = [&](size_t start, size_t next) {
    BRep_ImplicitPlaneSection::Curve curve;
    curve.Points.push_back(points[start]);
    size_t current = start;
    while (!used[next])
    {
      used[next] = true;
      current = segments[next].first == current ? segments[next].second : segments[next].first;
      curve.Points.push_back(points[current]);
      if (current == start) { curve.Closed = true; break; }
      bool found = false;
      for (auto candidate : adjacent[current]) if (!used[candidate]) { next = candidate; found = true; break; }
      if (!found) break;
    }
    result.Curves.push_back(std::move(curve));
  };
  for (size_t i = 0; i < adjacent.size(); ++i)
    if (adjacent[i].size() == 1 && !used[adjacent[i][0]]) walk(i, adjacent[i][0]);
  for (size_t i = 0; i < segments.size(); ++i) if (!used[i]) walk(segments[i].first, i);
  result.Resolved = result.UnresolvedCells == 0 && result.GraphDefects == 0;
  return result;
}
} // namespace

Section::Result Section::Trace(const TopoDS_Face& face, const Settings& s)
{
  const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  if (data.IsNull() || !data->ExcludedDomains().empty() || typeid(*data) != typeid(BRep_TImplicitFace)
      || typeid(*data->Field()) != typeid(Geom_ImplicitGyroid) || TopoDS_Iterator(face).More())
    throw Standard_NotImplemented("Cylinder section requires a childless exact native gyroid face");
  require(cells(s.MinimumCells) && cells(s.MaximumCells) && s.MinimumCells <= s.MaximumCells,
          "Invalid cylinder section resolution budget");
  for (double value : {s.Radius, s.AngleMinimum, s.AngleMaximum, s.HeightMinimum, s.HeightMaximum,
                       s.DistanceTolerance}) require(std::isfinite(value), "Non-finite cylinder section setting");
  require(s.Radius > 0 && s.AngleMaximum > s.AngleMinimum
          && s.AngleMaximum - s.AngleMinimum < 2 * std::acos(-1.)
          && s.HeightMaximum > s.HeightMinimum && s.DistanceTolerance > 0,
          "Invalid finite open cylinder section");
  require(std::abs(s.AngleMinimum) < 1.e6 && std::abs(s.AngleMaximum) < 1.e6,
          "Cylinder angles exceed trigonometric precision policy");
  const auto transform = (face.Location() * data->Location()).Transformation();
  const double scale = std::abs(transform.ScaleFactor());
  const auto* field = static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  const double k = 2 * std::acos(-1.) / field->Period();
  const double lipschitz = field->LipschitzBound() / scale;
  const double hessian = 4 * (k / scale) * (k / scale) + lipschitz / s.Radius;
  double magnitude = 1;
  // Conservative oriented square prism encloses every chart point, including
  // angular extrema. Checking its corners in the field's affine box is sufficient.
  for (double x : {-s.Radius, s.Radius}) for (double y : {-s.Radius, s.Radius})
    for (double z : {s.HeightMinimum, s.HeightMaximum})
    {
      const auto p = s.Frame.Location().Translated(gp_Vec(s.Frame.XDirection()) * x
        + gp_Vec(s.Frame.YDirection()) * y + gp_Vec(s.Frame.Direction()) * z);
      require(BRep_ImplicitTool::IsInDomain(face, p), "Cylinder enclosure outside field domain");
      const auto local = p.Transformed(transform.Inverted());
      for (int axis = 1; axis <= 3; ++axis)
      {
        require(std::isfinite(local.Coord(axis) * k) && std::abs(local.Coord(axis) * k) < 1.e6,
                "Cylinder coordinates exceed trigonometric precision policy");
        magnitude = std::max(magnitude, std::abs(p.Coord(axis)));
      }
    }
  const double epsilon = std::numeric_limits<double>::epsilon();
  const double precision = 512 * epsilon * magnitude * (1 + std::abs(s.AngleMinimum) + std::abs(s.AngleMaximum));
  const double reserve = 512 * epsilon * (1 + lipschitz + hessian + std::abs(field->IsoValue()))
                         + (lipschitz + hessian) * precision;
  require(std::isfinite(hessian) && std::isfinite(reserve) && s.DistanceTolerance > precision
          && std::isfinite(s.Radius * s.AngleMinimum) && std::isfinite(s.Radius * s.AngleMaximum)
          && std::isfinite(s.Radius * (s.AngleMaximum - s.AngleMinimum))
          && std::isfinite(s.HeightMaximum - s.HeightMinimum), "Cylinder precision/span overflow");
  Result result;
  for (int n = s.MinimumCells; n <= s.MaximumCells; n *= 2)
  {
    result = pass(face, s, n, lipschitz, hessian, reserve, precision);
    if (result.Resolved) return result;
  }
  return result;
}
