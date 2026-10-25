// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRep_ImplicitCylinderSection.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Line.hxx>
#include <Geom2d_Line.hxx>
#include <BRep_TImplicitFace.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopExp.hxx>
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
                     int n, double lipschitz, double hessian, double reserve, double precision,
                     bool periodic, Section::BoundaryResult* details = nullptr,
                     Section::PeriodicGraph* graph = nullptr)
{
  if (graph)
  {
    graph->Vertices.clear(); graph->Arcs.clear(); graph->CurveVertices.clear(); graph->CurveUV.clear();
    graph->Verified = false;
  }
  if (details)
  {
    details->Roots.clear(); details->CurveRoots.clear();
    details->BoundaryIncidenceVerified = false;
  }
  Section::Result result;
  result.Resolution = n;
  const Chart chart {s};
  const double u0 = s.Radius * s.AngleMinimum;
  const double du = s.Radius * (periodic ? 2 * std::acos(-1.) : s.AngleMaximum - s.AngleMinimum) / n;
  const double dv = (s.HeightMaximum - s.HeightMinimum) / n;
  const double diagonal = std::hypot(du, dv);
  const size_t width = size_t(n) + 1;
  std::vector<double> values(width * width);
  for (int j = 0; j <= n; ++j)
    for (int i = 0; i <= n; ++i)
    {
      if ((size_t(j) * width + i) % 256 == 0) checkpoint(s);
      if (periodic && i == n)
      {
        values[size_t(j) * width + i] = values[size_t(j) * width];
        continue;
      }
      values[size_t(j) * width + i] = BRep_ImplicitTool::Value(
        face, chart.point(u0 + du * i, s.HeightMinimum + dv * j));
    }
  std::vector<gp_Pnt> points;
  std::vector<bool> boundary;
  std::vector<std::pair<size_t, Section::BoundaryRoot>> capRoots;
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
    const double rootTolerance = s.DistanceTolerance
      / (details && alongU && (j == 0 || j == n) ? 65536 : 32);
    for (int iteration = 0; iteration < 100 && (high - low) * length > rootTolerance;
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
    if (bracket > rootTolerance) { e.state = 2; return e; }
    const double t = low + (high - low) * .5;
    const auto p = chart.point(u + (alongU ? du * t : 0), v + (alongU ? 0 : dv * t));
    e.state = 1; e.vertex = points.size(); points.push_back(p);
    if (graph)
    {
      Section::ChartVertex node;
      node.Point = p;
      node.UV.SetCoord((u + (alongU ? du * t : 0)) / s.Radius,
                       v + (alongU ? 0 : dv * t));
      node.Seam = periodic && !alongU && i == 0;
      if (node.Seam) BRep_Builder().MakeVertex(node.Vertex, p, bracket * .5 + precision);
      graph->Vertices.push_back(std::move(node));
    }
    boundary.push_back(alongU ? j == 0 || j == n : !periodic && (i == 0 || i == n));
    if (details && alongU && (j == 0 || j == n))
    {
      Section::BoundaryRoot root;
      root.AngleFirst = (u + du * low) / s.Radius;
      root.AngleLast = (u + du * high) / s.Radius;
      root.First = chart.point(u + du * low, v);
      root.Last = chart.point(u + du * high, v);
      root.HeightBoundary = j == 0 ? 0 : 1;
      BRep_Builder().MakeVertex(root.Vertex, p, bracket * .5 + precision);
      if (graph) graph->Vertices[e.vertex].Vertex = root.Vertex;
      capRoots.emplace_back(e.vertex, std::move(root));
    }
    result.MaximumResidual = std::max(result.MaximumResidual, std::abs(BRep_ImplicitTool::Value(face, p)));
    result.MaximumRootBracket = std::max(result.MaximumRootBracket, bracket);
    return e;
  };
  for (int j = 0; j <= n; ++j)
    for (int i = 0; i < n; ++i) horizontal[size_t(j) * n + i] = edge(i, j, true);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i <= n; ++i)
      vertical[size_t(j) * width + i] = periodic && i == n
        ? vertical[size_t(j) * width] : edge(i, j, false);
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
        if (graph) graph->Arcs.push_back({roots[0], roots[1], i, j});
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
  std::vector<std::array<size_t, 2>> curveEnds;
  const auto walk = [&](size_t start, size_t next) {
    BRep_ImplicitPlaneSection::Curve curve;
    curve.Points.push_back(points[start]);
    std::vector<size_t> nodes {start};
    std::vector<gp_Pnt2d> uv;
    size_t current = start;
    while (!used[next])
    {
      if (graph)
      {
        const auto& arc = graph->Arcs[next];
        const size_t destination = segments[next].first == current
          ? segments[next].second : segments[next].first;
        const double turn = 2 * std::acos(-1.);
        const auto local = [&](size_t id) {
          auto p = graph->Vertices[id].UV;
          if (graph->Vertices[id].Seam && arc.U == n - 1) p.SetX(p.X() + turn);
          return p;
        };
        auto a = local(current), b = local(destination);
        if (uv.empty()) uv.push_back(a);
        // Branch translation uses the same logical node, never point proximity.
        const double shift = turn * std::round((uv.back().X() - a.X()) / turn);
        b.SetX(b.X() + shift);
        uv.push_back(b);
      }
      used[next] = true;
      current = segments[next].first == current ? segments[next].second : segments[next].first;
      curve.Points.push_back(points[current]);
      nodes.push_back(current);
      if (current == start) { curve.Closed = true; break; }
      bool found = false;
      for (auto candidate : adjacent[current]) if (!used[candidate]) { next = candidate; found = true; break; }
      if (!found) break;
    }
    result.Curves.push_back(std::move(curve));
    curveEnds.push_back({start, current});
    if (graph)
    {
      graph->CurveVertices.push_back(std::move(nodes));
      graph->CurveUV.push_back(std::move(uv));
    }
  };
  for (size_t i = 0; i < adjacent.size(); ++i)
    if (adjacent[i].size() == 1 && !used[adjacent[i][0]]) walk(i, adjacent[i][0]);
  for (size_t i = 0; i < segments.size(); ++i) if (!used[i]) walk(segments[i].first, i);
  result.Resolved = result.UnresolvedCells == 0 && result.GraphDefects == 0;
  if (details && result.Resolved)
  {
    std::sort(capRoots.begin(), capRoots.end(), [](const auto& a, const auto& b) {
      return a.second.HeightBoundary != b.second.HeightBoundary
        ? a.second.HeightBoundary < b.second.HeightBoundary
        : a.second.AngleFirst < b.second.AngleFirst;
    });
    std::vector<size_t> ids(points.size(), 0), uses(capRoots.size(), 0);
    for (size_t index = 0; index < capRoots.size(); ++index)
    {
      ids[capRoots[index].first] = index + 1;
      details->Roots.push_back(capRoots[index].second);
    }
    bool verified = true;
    for (size_t index = 0; index < curveEnds.size(); ++index)
    {
      const auto& ends = curveEnds[index];
      const std::array<size_t, 2> roots {{ids[ends[0]], ids[ends[1]]}};
      if (result.Curves[index].Closed)
        verified &= roots[0] == 0 && roots[1] == 0;
      else
        for (auto id : roots)
        {
          verified &= id != 0;
          if (id) ++uses[id - 1];
        }
      details->CurveRoots.push_back(roots);
    }
    for (auto count : uses) verified &= count == 1;
    if (!verified)
    {
      result.Resolved = false; ++result.GraphDefects;
      details->Roots.clear(); details->CurveRoots.clear();
    }
    details->BoundaryIncidenceVerified = verified;
  }
  if (graph)
  {
    graph->Verified = result.Resolved && details && details->BoundaryIncidenceVerified;
    if (!graph->Verified)
    {
      graph->Vertices.clear(); graph->Arcs.clear(); graph->CurveVertices.clear(); graph->CurveUV.clear();
    }
  }
  return result;
}
} // namespace

static Section::Result traceCylinder(const TopoDS_Face& face, const Section::Settings& s, bool periodic,
                                     Section::BoundaryResult* details = nullptr,
                                     Section::PeriodicGraph* graph = nullptr)
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
          && (periodic ? s.AngleMaximum == s.AngleMinimum + 2 * std::acos(-1.)
                       : s.AngleMaximum - s.AngleMinimum < 2 * std::acos(-1.))
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
  Section::Result result;
  for (int n = s.MinimumCells; n <= s.MaximumCells; n *= 2)
  {
    result = pass(face, s, n, lipschitz, hessian, reserve, precision, periodic, details, graph);
    if (result.Resolved) return result;
  }
  if (periodic)
    for (int numerator : {6, 5, 7})
    {
      const int n = s.MaximumCells * numerator / 8;
      if (n < s.MinimumCells) continue;
      auto candidate = pass(face, s, n, lipschitz, hessian, reserve, precision, true, details, graph);
      if (candidate.Resolved) return candidate;
    }
  return result;
}

Section::Result Section::Trace(const TopoDS_Face& face, const Settings& settings)
{
  return traceCylinder(face, settings, false);
}

Section::Result Section::TracePeriodic(const TopoDS_Face& face, const Settings& settings)
{
  return traceCylinder(face, settings, true);
}

Section::BoundaryResult Section::TracePeriodicBoundary(const TopoDS_Face& face, const Settings& settings)
{
  BoundaryResult result;
  result.Section = traceCylinder(face, settings, true, &result);
  if (!result.Section.Resolved)
  {
    result.Roots.clear(); result.CurveRoots.clear();
    result.BoundaryIncidenceVerified = false;
  }
  return result;
}

Section::PeriodicGraph Section::TracePeriodicGraph(const TopoDS_Face& face, const Settings& settings)
{
  PeriodicGraph graph;
  graph.Boundary.Section = traceCylinder(face, settings, true, &graph.Boundary, &graph);
  if (!graph.Boundary.Section.Resolved)
  {
    graph.Boundary.Roots.clear(); graph.Boundary.CurveRoots.clear();
    graph.Boundary.BoundaryIncidenceVerified = false;
    graph.Vertices.clear(); graph.Arcs.clear(); graph.CurveVertices.clear(); graph.CurveUV.clear();
    graph.Verified = false;
  }
  return graph;
}

static Section::RingResult buildBoundaryRings(const TopoDS_Face& face,
  const Section::Settings& s, const Section::BoundaryResult& boundary)
{
  Section::RingResult result;
  result.Boundary = boundary;
  if (!result.Boundary.BoundaryIncidenceVerified) return result;
  const Chart chart {s};
  const double turn = 2 * std::acos(-1.);
  const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  const auto* field = static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  const double scale = std::abs((face.Location() * data->Location()).Transformation().ScaleFactor());
  const double lipschitz = field->LipschitzBound() / scale;
  BRep_Builder builder;
  bool verified = true;
  for (int cap = 0; cap < 2; ++cap)
  {
    checkpoint(s);
    const double height = cap == 0 ? s.HeightMinimum : s.HeightMaximum;
    const gp_Pnt origin = s.Frame.Location().Translated(gp_Vec(s.Frame.Direction()) * height);
    const occ::handle<Geom_Circle> circle = new Geom_Circle(
      gp_Ax2(origin, s.Frame.Direction(), s.Frame.XDirection()), s.Radius);
    std::vector<size_t> roots;
    for (size_t i = 0; i < result.Boundary.Roots.size(); ++i)
      if (result.Boundary.Roots[i].HeightBoundary == cap) roots.push_back(i + 1);
    if (roots.size() % 2 != 0) { verified = false; break; }
    std::vector<bool> signs;
    const size_t count = roots.empty() ? 2 : roots.size();
    std::array<TopoDS_Vertex, 2> synthetic;
    const double magnitude = 1 + std::max({std::abs(origin.X()), std::abs(origin.Y()),
                                           std::abs(origin.Z()), s.Radius});
    const double precision = 1024 * std::numeric_limits<double>::epsilon() * magnitude
                             * (1 + std::abs(s.AngleMinimum) + std::abs(s.AngleMaximum));
    const double reserve = 1024 * std::numeric_limits<double>::epsilon()
                             * (1 + std::abs(field->IsoValue())) + lipschitz * precision;
    if (!std::isfinite(precision) || !std::isfinite(reserve) || precision > s.DistanceTolerance)
    { verified = false; break; }
    if (roots.empty())
      for (int i = 0; i < 2; ++i)
        builder.MakeVertex(synthetic[i], circle->Value(s.AngleMinimum + i * turn / 2), precision);
    for (size_t i = 0; i < count; ++i)
    {
      checkpoint(s);
      const size_t next = (i + 1) % count;
      Section::RingArc arc;
      double first, last, low, high;
      TopoDS_Vertex a, b;
      if (roots.empty())
      {
        first = low = s.AngleMinimum + i * turn / 2;
        last = high = first + turn / 2;
        a = synthetic[i]; b = synthetic[next];
      }
      else
      {
        arc.LowerRoot = roots[i]; arc.UpperRoot = roots[next];
        const auto& lower = result.Boundary.Roots[arc.LowerRoot - 1];
        const auto& upper = result.Boundary.Roots[arc.UpperRoot - 1];
        first = (lower.AngleFirst + lower.AngleLast) * .5;
        last = (upper.AngleFirst + upper.AngleLast) * .5 + (next == 0 ? turn : 0);
        low = lower.AngleLast;
        high = upper.AngleFirst + (next == 0 ? turn : 0);
        a = lower.Vertex; b = upper.Vertex;
      }
      if (!(high > low) || s.Radius * (high - low) <= precision)
      { verified = false; break; }
      const double sample = low + (high - low) * .5;
      const double value = BRep_ImplicitTool::Value(face, chart.point(s.Radius * sample, height));
      if (!std::isfinite(value) || std::abs(value) <= reserve)
      { verified = false; break; }
      arc.Negative = value < 0; signs.push_back(arc.Negative);
      const double tolerance = std::max({precision, BRep_Tool::Tolerance(a), BRep_Tool::Tolerance(b)});
      if (tolerance > s.DistanceTolerance) { verified = false; break; }
      builder.MakeEdge(arc.Edge, circle, tolerance);
      a.Orientation(TopAbs_FORWARD); b.Orientation(TopAbs_REVERSED);
      builder.Add(arc.Edge, a); builder.Add(arc.Edge, b);
      builder.UpdateVertex(a, first, arc.Edge, BRep_Tool::Tolerance(a));
      builder.UpdateVertex(b, last, arc.Edge, BRep_Tool::Tolerance(b));
      builder.Range(arc.Edge, first, last);
      if (cap == 0) arc.Edge.Reverse();
      result.Rings[cap].push_back(std::move(arc));
    }
    if (!verified) break;
    for (size_t i = 0; i < signs.size(); ++i)
      verified &= roots.empty() ? signs[i] == signs[(i + 1) % signs.size()]
                                : signs[i] != signs[(i + 1) % signs.size()];
    if (!verified) break;
  }
  if (!verified) for (auto& ring : result.Rings) ring.clear();
  result.Verified = verified;
  return result;
}

Section::RingResult Section::BuildBoundaryRings(const TopoDS_Face& face, const Settings& s)
{
  return buildBoundaryRings(face, s, TracePeriodicBoundary(face, s));
}

Section::ChartSkeleton Section::BuildChartSkeleton(const TopoDS_Face& face, const Settings& s)
{
  ChartSkeleton result;
  result.Graph = TracePeriodicGraph(face, s);
  if (!result.Graph.Verified) return result;
  result.Rings = buildBoundaryRings(face, s, result.Graph.Boundary);
  if (!result.Rings.Verified) return result;
  BRep_Builder builder;
  for (int cap = 0; cap < 2; ++cap)
  {
    checkpoint(s);
    auto& arcs = result.Rings.Rings[cap];
    // Root-free rings already have two half circles meeting at the chart seam.
    if (arcs.front().LowerRoot == 0)
    {
      TopoDS_Vertex other;
      TopExp::Vertices(arcs.front().Edge, result.SeamEnds[cap], other);
      continue;
    }
    const auto wrap = arcs.back();
    double first, last;
    const auto circle = BRep_Tool::Curve(wrap.Edge, first, last);
    TopoDS_Vertex a, b;
    TopExp::Vertices(wrap.Edge, a, b);
    const double tolerance = BRep_Tool::Tolerance(wrap.Edge);
    builder.MakeVertex(result.SeamEnds[cap], circle->Value(s.AngleMaximum), tolerance);
    const auto make = [&](TopoDS_Vertex start, TopoDS_Vertex end,
                          double low, double high, size_t lower, size_t upper) {
      RingArc arc;
      arc.LowerRoot = lower; arc.UpperRoot = upper; arc.Negative = wrap.Negative;
      builder.MakeEdge(arc.Edge, circle, tolerance);
      start.Orientation(TopAbs_FORWARD); end.Orientation(TopAbs_REVERSED);
      builder.Add(arc.Edge, start); builder.Add(arc.Edge, end);
      builder.UpdateVertex(start, low, arc.Edge, BRep_Tool::Tolerance(start));
      builder.UpdateVertex(end, high, arc.Edge, BRep_Tool::Tolerance(end));
      builder.Range(arc.Edge, low, high);
      if (cap == 0) arc.Edge.Reverse();
      return arc;
    };
    arcs.back() = make(a, result.SeamEnds[cap], first, s.AngleMaximum, wrap.LowerRoot, 0);
    arcs.push_back(make(result.SeamEnds[cap], b, s.AngleMinimum,
                        last - 2 * std::acos(-1.), 0, wrap.UpperRoot));
  }
  result.Surface = new Geom_CylindricalSurface(s.Frame, s.Radius);
  for (int cap = 0; cap < 2; ++cap)
    for (const auto& arc : result.Rings.Rings[cap])
    {
      double first, last; BRep_Tool::Range(arc.Edge, first, last);
      const occ::handle<Geom2d_Line> pc = new Geom2d_Line(
        gp_Pnt2d(0, cap == 0 ? s.HeightMinimum : s.HeightMaximum), gp_Dir2d(1, 0));
      builder.UpdateEdge(arc.Edge, pc, result.Surface, TopLoc_Location(), BRep_Tool::Tolerance(arc.Edge));
      builder.Range(arc.Edge, result.Surface, TopLoc_Location(), first, last);
      builder.SameRange(arc.Edge, true); builder.SameParameter(arc.Edge, true);
    }
  std::vector<std::pair<double, TopoDS_Vertex>> seam {{s.HeightMinimum, result.SeamEnds[0]},
                                                    {s.HeightMaximum, result.SeamEnds[1]}};
  for (const auto& node : result.Graph.Vertices)
    if (node.Seam) seam.emplace_back(node.UV.Y(), node.Vertex);
  std::sort(seam.begin(), seam.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  const auto origin = Chart{s}.point(s.Radius * s.AngleMinimum, 0);
  const occ::handle<Geom_Line> line = new Geom_Line(origin, s.Frame.Direction());
  const occ::handle<Geom2d_Line> low = new Geom2d_Line(gp_Pnt2d(s.AngleMinimum, 0), gp_Dir2d(0, 1));
  const occ::handle<Geom2d_Line> high = new Geom2d_Line(gp_Pnt2d(s.AngleMaximum, 0), gp_Dir2d(0, 1));
  for (size_t i = 1; i < seam.size(); ++i)
  {
    checkpoint(s);
    auto a = seam[i - 1].second, b = seam[i].second;
    const double tolerance = std::max(BRep_Tool::Tolerance(a), BRep_Tool::Tolerance(b));
    TopoDS_Edge edge;
    builder.MakeEdge(edge, line, tolerance);
    a.Orientation(TopAbs_FORWARD); b.Orientation(TopAbs_REVERSED);
    builder.Add(edge, a); builder.Add(edge, b);
    builder.UpdateVertex(a, seam[i - 1].first, edge, BRep_Tool::Tolerance(a));
    builder.UpdateVertex(b, seam[i].first, edge, BRep_Tool::Tolerance(b));
    builder.Range(edge, seam[i - 1].first, seam[i].first);
    builder.UpdateEdge(edge, low, high, result.Surface, TopLoc_Location(), tolerance);
    builder.Range(edge, result.Surface, TopLoc_Location(), seam[i - 1].first, seam[i].first);
    builder.SameRange(edge, true); builder.SameParameter(edge, true);
    result.SeamEdges.push_back(edge);
  }
  result.Verified = true;
  return result;
}

Section::CapSkeleton Section::BuildCapSkeleton(const TopoDS_Face& face, const Settings& s)
{
  CapSkeleton result;
  result.Disks = TraceEndCapDisks(face, s);
  if (!result.Disks.Verified) return result;
  result.Rings = buildBoundaryRings(face, s, result.Disks.Boundary);
  result.Verified = result.Rings.Verified;
  return result;
}
