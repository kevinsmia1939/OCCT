// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRep_ImplicitPlaneSection.hxx>
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
#include <map>
#include <typeinfo>

namespace
{
void require(bool value, const char* message)
{
  if (!value)
    throw Standard_DomainError(message);
}

bool cells(int n)
{
  return n >= 4 && n <= 512 && (n & (n - 1)) == 0;
}

bool opposite(double a, double b)
{
  return std::signbit(a) != std::signbit(b);
}

struct Edge
{
  int    State  = 0;
  size_t Vertex = 0;
}; // 0 empty, 1 unique root, 2 unresolved

BRep_ImplicitPlaneSection::Result pass(const TopoDS_Face&                         face,
                                       const BRep_ImplicitPlaneSection::Settings& settings,
                                       int                                        n,
                                       double                                     lipschitz,
                                       double                                     hessian,
                                       double                                     fieldReserve,
                                       double                                     gradientReserve,
                                       double                                     coordinateReserve)
{
  using Section = BRep_ImplicitPlaneSection;
  Section::Result result;
  result.Resolution            = n;
  const size_t        width    = size_t(n) + 1;
  const double        du       = (settings.UMaximum - settings.UMinimum) / n;
  const double        dv       = (settings.VMaximum - settings.VMinimum) / n;
  const double        diagonal = std::hypot(du, dv);
  const gp_Vec        u(settings.Frame.XDirection()), v(settings.Frame.YDirection());
  std::vector<gp_Pnt> grid(width * width);
  std::vector<double> values(width * width);
  for (int j = 0; j <= n; ++j)
    for (int i = 0; i <= n; ++i)
    {
      const size_t index = size_t(j) * width + i;
      grid[index]        = settings.Frame.Location().Translated(u * (settings.UMinimum + du * i)
                                                                + v * (settings.VMinimum + dv * j));
      values[index]      = BRep_ImplicitTool::Value(face, grid[index]);
    }
  std::vector<gp_Pnt>                       vertices;
  std::vector<bool>                         boundary;
  std::vector<std::pair<size_t, size_t>>    segments;
  std::map<std::pair<size_t, size_t>, Edge> edgeCache;
  auto                                      edge = [&](size_t a, size_t b) {
    const auto key   = std::minmax(a, b);
    const auto found = edgeCache.find(key);
    if (found != edgeCache.end())
      return found->second;
    Edge         answer;
    const double fa = values[a], fb = values[b];
    const gp_Pnt mid((grid[a].XYZ() + grid[b].XYZ()) * .5);
    const gp_Vec delta(grid[a], grid[b]);
    const double length = delta.Magnitude();
    if (!std::isfinite(length) || length <= 0)
    {
      answer.State = 2; // Coordinate collapse is unresolved, never a zero-length division.
      edgeCache.emplace(key, answer);
      return answer;
    }
    const double derivative = BRep_ImplicitTool::Gradient(face, mid).Dot(delta / length);
    const bool   monotone   = std::abs(derivative) > hessian * length * .5 + gradientReserve;
    if (std::abs(fa) <= fieldReserve || std::abs(fb) <= fieldReserve)
      answer.State = 2; // Grid-aligned zero/near-zero: do not guess its connectivity.
    else if (!opposite(fa, fb)
             && (monotone
                 || std::abs(BRep_ImplicitTool::Value(face, mid))
                      > (std::abs(derivative) + gradientReserve) * length * .5
                          + hessian * length * length * .125 + fieldReserve))
      answer.State = 0;
    else if (!monotone || !opposite(fa, fb))
      answer.State = 2;
    else
    {
      double low = 0, high = 1, fLow = fa;
      for (int iteration = 0;
           iteration < 100 && (high - low) * length > settings.DistanceTolerance / 16;
           ++iteration)
      {
        const double t = low + (high - low) * .5;
        if (t == low || t == high)
          break;
        const double f = BRep_ImplicitTool::Value(face, grid[a].Translated(delta * t));
        if (f == 0)
          low = high = t;
        else if (opposite(fLow, f))
          high = t;
        else
        {
          low  = t;
          fLow = f;
        }
      }
      const double bracket = (high - low) * length;
      if (bracket > settings.DistanceTolerance / 16)
        answer.State = 2;
      else
      {
        answer.State  = 1;
        answer.Vertex = vertices.size();
        const auto p  = grid[a].Translated(delta * (low + (high - low) * .5));
        vertices.push_back(p);
        const bool horizontal = a / width == b / width;
        boundary.push_back(horizontal ? (a / width == 0 || a / width == size_t(n))
                                      : (a % width == 0 || a % width == size_t(n)));
        result.MaximumResidual =
          std::max(result.MaximumResidual, std::abs(BRep_ImplicitTool::Value(face, p)));
        result.MaximumRootBracket = std::max(result.MaximumRootBracket, bracket);
      }
    }
    edgeCache.emplace(key, answer);
    return answer;
  };
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
    {
      const size_t                a       = size_t(j) * width + i;
      const std::array<size_t, 4> corners = {a, a + 1, a + 1 + width, a + width};
      const gp_Pnt                center((grid[a].XYZ() + grid[a + width + 1].XYZ()) * .5);
      const double                f = BRep_ImplicitTool::Value(face, center);
      if (std::abs(f) > lipschitz * diagonal * .5 + fieldReserve)
        continue;
      const auto   g           = BRep_ImplicitTool::Gradient(face, center);
      const double uncertainty = hessian * diagonal * .5 + gradientReserve;
      const bool   regular = std::abs(g.Dot(u)) > uncertainty || std::abs(g.Dot(v)) > uncertainty;
      std::array<size_t, 4> roots;
      int                   count     = 0;
      bool                  ambiguous = false;
      for (int k = 0; k < 4; ++k)
      {
        const auto crossing = edge(corners[k], corners[(k + 1) % 4]);
        ambiguous |= crossing.State == 2;
        if (crossing.State == 1)
          roots[count++] = crossing.Vertex;
      }
      // A partial derivative of fixed nonzero sign rules out closed internal
      // components. Unique edge roots + two crossings then imply one arc.
      if (regular && !ambiguous && count == 0)
        continue;
      if (regular && !ambiguous && count == 2
          && diagonal + coordinateReserve <= settings.DistanceTolerance)
      {
        segments.emplace_back(roots[0], roots[1]);
        result.GeometricBound = std::max(result.GeometricBound, diagonal + coordinateReserve);
      }
      else
      {
        ++result.UnresolvedCells;
        if (!regular)
          ++result.PotentialCriticalCells;
        if (ambiguous || count > 2)
          ++result.AmbiguousCells;
      }
    }
  std::vector<std::vector<size_t>> adjacent(vertices.size());
  for (size_t i = 0; i < segments.size(); ++i)
  {
    adjacent[segments[i].first].push_back(i);
    adjacent[segments[i].second].push_back(i);
  }
  for (size_t i = 0; i < adjacent.size(); ++i)
    if (adjacent[i].size() > 2 || (adjacent[i].size() == 1 && !boundary[i]))
      ++result.GraphDefects;
  std::vector<bool> used(segments.size(), false);
  auto              walk = [&](size_t start, size_t next) {
    Section::Curve curve;
    size_t         current = start;
    curve.Points.push_back(vertices[start]);
    while (!used[next])
    {
      used[next]          = true;
      const auto& segment = segments[next];
      current             = segment.first == current ? segment.second : segment.first;
      curve.Points.push_back(vertices[current]);
      if (current == start)
      {
        curve.Closed = true;
        break;
      }
      bool found = false;
      for (auto candidate : adjacent[current])
        if (!used[candidate])
        {
          next  = candidate;
          found = true;
          break;
        }
      if (!found)
        break;
    }
    result.Curves.push_back(std::move(curve));
  };
  for (size_t i = 0; i < adjacent.size(); ++i)
    if (adjacent[i].size() == 1 && !used[adjacent[i][0]])
      walk(i, adjacent[i][0]);
  for (size_t i = 0; i < segments.size(); ++i)
    if (!used[i])
      walk(segments[i].first, i);
  result.Resolved = result.UnresolvedCells == 0 && result.GraphDefects == 0;
  return result;
}
} // namespace

BRep_ImplicitPlaneSection::Result BRep_ImplicitPlaneSection::Trace(const TopoDS_Face& face,
                                                                   const Settings&    settings)
{
  const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  if (data.IsNull() || typeid(*data) != typeid(BRep_TImplicitFace)
      || typeid(*data->Field()) != typeid(Geom_ImplicitGyroid) || TopoDS_Iterator(face).More())
    throw Standard_NotImplemented("Plane section requires a childless exact native gyroid face");
  require(cells(settings.MinimumCells) && cells(settings.MaximumCells)
            && settings.MinimumCells <= settings.MaximumCells,
          "Invalid section resolution budget");
  for (double value : {settings.UMinimum,
                       settings.UMaximum,
                       settings.VMinimum,
                       settings.VMaximum,
                       settings.DistanceTolerance})
    require(std::isfinite(value), "Non-finite section bounds/tolerance");
  require(settings.UMinimum < settings.UMaximum && settings.VMinimum < settings.VMaximum
            && settings.DistanceTolerance > 0,
          "Invalid section bounds/tolerance");
  require(std::isfinite(settings.UMaximum - settings.UMinimum)
            && std::isfinite(settings.VMaximum - settings.VMinimum),
          "Section span overflow");
  const gp_Trsf transform = (face.Location() * data->Location()).Transformation();
  const double  scale     = std::abs(transform.ScaleFactor());
  require(std::isfinite(scale) && scale > 0, "Invalid section placement");
  const auto*  field = static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  const double k     = 2 * std::acos(-1.0) / field->Period();
  // Hessian diagonal bounds 2*k^2, off-diagonal bounds k^2: symmetric
  // maximum row sum bounds its operator norm by 4*k^2, before placement.
  const double lipschitz = field->LipschitzBound() / scale, hessian = 4 * (k / scale) * (k / scale);
  require(std::isfinite(lipschitz) && std::isfinite(hessian) && hessian > 0,
          "Unrepresentable derivative bounds");
  double magnitude = 1;
  for (double u : {settings.UMinimum, settings.UMaximum})
    for (double v : {settings.VMinimum, settings.VMaximum})
    {
      const auto p = settings.Frame.Location().Translated(
        gp_Vec(settings.Frame.XDirection()) * u + gp_Vec(settings.Frame.YDirection()) * v);
      require(BRep_ImplicitTool::IsInDomain(face, p),
              "Section rectangle is outside the field domain");
      const auto local = p.Transformed(transform.Inverted());
      for (int axis = 1; axis <= 3; ++axis)
      {
        require(std::isfinite(local.Coord(axis) * k) && std::abs(local.Coord(axis) * k) < 1.e6,
                "Section coordinates exceed the trigonometric precision policy");
        magnitude = std::max(magnitude, std::abs(p.Coord(axis)));
      }
    }
  const double epsilon           = std::numeric_limits<double>::epsilon();
  const double coordinateReserve = 128 * epsilon * magnitude;
  require(settings.DistanceTolerance > coordinateReserve,
          "Section tolerance below coordinate precision");
  const double fieldReserve =
    128 * epsilon * (1 + std::abs(field->IsoValue())) + lipschitz * coordinateReserve;
  const double gradientReserve = 128 * epsilon * lipschitz + hessian * coordinateReserve;
  require(std::isfinite(fieldReserve) && std::isfinite(gradientReserve),
          "Section reserve overflow");
  Result result;
  for (int n = settings.MinimumCells; n <= settings.MaximumCells; n *= 2)
  {
    result =
      pass(face, settings, n, lipschitz, hessian, fieldReserve, gradientReserve, coordinateReserve);
    if (result.Resolved)
      break;
  }
  return result;
}
