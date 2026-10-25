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
void checkpoint(const BRep_ImplicitPlaneSection::Settings& settings)
{
  if (settings.WorkCheckpoint)
    settings.WorkCheckpoint();
}

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

bool exactSum(double a, double b, double sum)
{
  // TwoSum residual: an affine map rounded ONTO a lattice point is not an
  // exact lattice identity of the requested parameter-space point.
  const double part = sum - a;
  return (a - (sum - part)) + (b - part) == 0;
}

struct Edge
{
  int    State  = 0;
  size_t Vertex = 0;
}; // 0 empty, 1 unique root, 2 unresolved

BRep_ImplicitPlaneSection::Result pass(const TopoDS_Face&                         face,
                                       const BRep_ImplicitPlaneSection::Settings& settings,
                                       int                                        n,
                                       const Geom_ImplicitGyroid&                 field,
                                       double                                     lipschitz,
                                       double                                     hessian,
                                       double                                     fieldReserve,
                                       double                                     gradientReserve,
                                       double                                     coordinateReserve,
                                       BRep_ImplicitPlaneSection::GraphResult* graph = nullptr)
{
  if (graph) { graph->Vertices.clear(); graph->Arcs.clear(); graph->VertexRadius = 0; }
  using Section = BRep_ImplicitPlaneSection;
  checkpoint(settings);
  Section::Result result;
  result.Resolution            = n;
  const size_t        width    = size_t(n) + 1;
  const double        du       = (settings.UMaximum - settings.UMinimum) / n;
  const double        dv       = (settings.VMaximum - settings.VMinimum) / n;
  const double        diagonal = std::hypot(du, dv);
  const gp_Vec        u(settings.Frame.XDirection()), v(settings.Frame.YDirection());
  const auto          data       = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  const auto          inverse    = (face.Location() * data->Location()).Transformation().Inverted();
  const double        waveNumber = 2 * std::acos(-1.) / field.Period();
  std::vector<gp_Pnt> grid(width * width);
  std::vector<double> values(width * width);
  std::vector<bool> exactZero(width * width, false);
  for (int j = 0; j <= n; ++j)
    for (int i = 0; i <= n; ++i)
    {
      const size_t index = size_t(j) * width + i;
      if (index % 256 == 0)
        checkpoint(settings);
      const double x = i == n ? settings.UMaximum : settings.UMinimum + du * i;
      const double y = j == n ? settings.VMaximum : settings.VMinimum + dv * j;
      grid[index]        = settings.Frame.Location().Translated(u * x + v * y);
      values[index]      = BRep_ImplicitTool::Value(face, grid[index]);
      if (settings.AllowExactLatticeNodes)
      {
        double value;
        gp_Vec gradient;
        bool exactMap = true;
        for (int axis = 1; axis <= 3; ++axis)
        {
          const double offset = u.Coord(axis) != 0 ? u.Coord(axis) * x : v.Coord(axis) * y;
          exactMap &= exactSum(settings.Frame.Location().Coord(axis), offset, grid[index].Coord(axis));
        }
        exactZero[index] = exactMap && field.LatticeJet(grid[index], value, gradient) && value == 0;
      }
    }
  std::vector<gp_Pnt>                       vertices;
  std::vector<bool>                         boundary;
  std::vector<std::pair<size_t, size_t>>    segments;
  std::map<std::pair<size_t, size_t>, Edge> edgeCache;
  const size_t invalidVertex = std::numeric_limits<size_t>::max();
  std::vector<size_t> nodeVertex(grid.size(), invalidVertex);
  std::vector<bool> isolatedContact(grid.size(), false);
  if (settings.AllowExactLatticeNodes)
    for (size_t index = 0; index < grid.size(); ++index)
      if (exactZero[index])
      {
        nodeVertex[index] = vertices.size();
        vertices.push_back(grid[index]);
        boundary.push_back(index / width == 0 || index / width == size_t(n)
                           || index % width == 0 || index % width == size_t(n));
        result.MaximumResidual = std::max(result.MaximumResidual, std::abs(values[index]));
        values[index] = 0; // Analytic identity only, never a near-zero residual snap.
        ++result.ExactLatticeNodeCount;
      }
  auto                                      edge = [&](size_t a, size_t b) {
    checkpoint(settings);
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
    const double derivative     = BRep_ImplicitTool::Gradient(face, mid).Dot(delta / length);
    const bool   monotone       = std::abs(derivative) > hessian * length * .5 + gradientReserve;
    bool         harmonicUnique = false, harmonicEmpty = false;
    const auto   localDelta  = delta.Transformed(inverse);
    int          varyingAxis = -1, components = 0;
    for (int axis = 1; axis <= 3; ++axis)
      if (localDelta.Coord(axis) != 0)
      {
        varyingAxis = axis;
        ++components;
      }
    if (components == 1 && waveNumber * localDelta.Magnitude() < std::acos(-1.))
    {
      // On a field-coordinate axis the gyroid is exactly A sin(theta) +
      // B cos(theta) + C. An interval shorter than pi has at most two roots,
      // so opposite endpoint signs imply exactly one transverse root even
      // if the derivative changes sign elsewhere in that interval.
      harmonicUnique = opposite(fa, fb);
      if (!harmonicUnique)
      {
        const auto   p = grid[a].Transformed(inverse), q = grid[b].Transformed(inverse);
        const int    next = varyingAxis % 3 + 1, prev = (varyingAxis + 1) % 3 + 1;
        const double A = std::cos(waveNumber * p.Coord(next));
        const double B = std::sin(waveNumber * p.Coord(prev));
        const double C = std::sin(waveNumber * p.Coord(next)) * std::cos(waveNumber * p.Coord(prev))
                         - field.IsoValue();
        const double pi = std::acos(-1.), base = std::atan2(A, B);
        const double low  = waveNumber * std::min(p.Coord(varyingAxis), q.Coord(varyingAxis));
        const double high = waveNumber * std::max(p.Coord(varyingAxis), q.Coord(varyingAxis));
        harmonicEmpty     = true;
        for (int index = int(std::ceil((low - base) / pi));
             index <= int(std::floor((high - base) / pi));
             ++index)
        {
          const double theta = base + index * pi;
          const double value = A * std::sin(theta) + B * std::cos(theta) + C;
          if (std::abs(value) <= 4 * fieldReserve || opposite(value, fa))
            harmonicEmpty = false;
        }
      }
    }
    if (exactZero[a] != exactZero[b])
    {
      const bool first = exactZero[a];
      const double other = first ? fb : fa;
      // A unique endpoint root needs monotonicity on the ENTIRE grid edge.
      // Harmonic opposite-sign uniqueness is not an endpoint/tangency proof.
      if (monotone && std::abs(other) > fieldReserve
          && std::signbit(other) == std::signbit(first ? derivative : -derivative))
      {
        answer.State = 1;
        answer.Vertex = nodeVertex[first ? a : b];
      }
      else
        answer.State = 2;
    }
    else if (std::abs(fa) <= fieldReserve || std::abs(fb) <= fieldReserve)
      answer.State = 2; // Grid-aligned zero/near-zero: do not guess its connectivity.
    else if (!opposite(fa, fb)
             && (monotone || harmonicEmpty
                 || std::abs(BRep_ImplicitTool::Value(face, mid))
                      > (std::abs(derivative) + gradientReserve) * length * .5
                          + hessian * length * length * .125 + fieldReserve))
      answer.State = 0;
    else if ((!monotone && !harmonicUnique) || !opposite(fa, fb))
      answer.State = 2;
    else
    {
      double low = 0, high = 1, fLow = fa;
      const double rootTolerance = settings.DistanceTolerance / (graph ? 4096 : 16);
      for (int iteration = 0;
           iteration < 100 && (high - low) * length > rootTolerance;
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
      if (bracket > rootTolerance)
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
      if (a % 256 == 0)
        checkpoint(settings);
      const std::array<size_t, 4> corners = {a, a + 1, a + 1 + width, a + width};
      const gp_Pnt                center((grid[a].XYZ() + grid[a + width + 1].XYZ()) * .5);
      const double                f = BRep_ImplicitTool::Value(face, center);
      if (std::abs(f) > lipschitz * diagonal * .5 + fieldReserve)
        continue;
      const auto g = BRep_ImplicitTool::Gradient(face, center);
      // Taylor exclusion uses the in-plane gradient and the existing global
      // Hessian operator bound. Unlike L*r alone it can certify empty cells
      // near a stationary point whose field value is bounded away from zero.
      // Never infer emptiness merely from equal signs at the corners.
      const double radius     = diagonal * .5;
      const double localBound = (std::hypot(g.Dot(u), g.Dot(v)) + gradientReserve) * radius
                                + hessian * radius * radius * .5 + fieldReserve;
      if (std::abs(f) > localBound)
        continue;
      const double uncertainty = hessian * diagonal * .5 + gradientReserve;
      const bool   regular = std::abs(g.Dot(u)) > uncertainty || std::abs(g.Dot(v)) > uncertainty;
      if (settings.AllowExactLatticeNodes
          && std::abs(g.Dot(u)) > uncertainty && std::abs(g.Dot(v)) > uncertainty)
      {
        int zero = -1, count = 0;
        for (int k = 0; k < 4; ++k)
          if (exactZero[corners[k]])
          {
            zero = k;
            ++count;
          }
        if (count == 1)
        {
          const double gu = g.Dot(u) * (zero == 1 || zero == 2 ? -1 : 1);
          const double gv = g.Dot(v) * (zero >= 2 ? -1 : 1);
          if (std::signbit(gu) == std::signbit(gv))
          {
            // Both inward derivatives keep this sign throughout the cell.
            // Thus only its exact zero corner meets F=0: no one-dimensional
            // arc enters this quadrant. Other incident cells may have an arc.
            const size_t node = corners[size_t(zero)];
            const bool rectangleCorner = (node / width == 0 || node / width == size_t(n))
                                      && (node % width == 0 || node % width == size_t(n));
            isolatedContact[node] = rectangleCorner;
            continue;
          }
        }
      }
      std::array<size_t, 4> roots;
      int                   count     = 0;
      bool                  ambiguous = false;
      for (int k = 0; k < 4; ++k)
      {
        const auto crossing = edge(corners[k], corners[(k + 1) % 4]);
        ambiguous |= crossing.State == 2;
        if (crossing.State == 1)
          if (std::find(roots.begin(), roots.begin() + count, crossing.Vertex) == roots.begin() + count)
            roots[count++] = crossing.Vertex; // A node shared by two cell edges is ONE root.
      }
      // A partial derivative of fixed nonzero sign rules out closed internal
      // components. Unique edge roots + two crossings then imply one arc.
      if (regular && !ambiguous && count == 0)
        continue;
      if (regular && !ambiguous && count == 2
          && diagonal + coordinateReserve <= settings.DistanceTolerance)
      {
        segments.emplace_back(roots[0], roots[1]);
        if (graph) graph->Arcs.push_back({i, j, roots[0], roots[1]});
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
  if (settings.AllowExactLatticeNodes)
    for (size_t node = 0; node < nodeVertex.size(); ++node)
      if (nodeVertex[node] != invalidVertex && adjacent[nodeVertex[node]].empty())
      {
        if (isolatedContact[node])
          result.PointContacts.push_back(grid[node]);
        else
          ++result.GraphDefects; // An unaccounted exact zero cannot disappear.
      }
  std::vector<bool> used(segments.size(), false);
  auto              walk = [&](size_t start, size_t next) {
    Section::Curve curve;
    size_t         current = start;
    curve.Points.push_back(vertices[start]);
    while (!used[next])
    {
      checkpoint(settings);
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
  if (graph && result.Resolved)
  {
    graph->Vertices = vertices;
    graph->VertexRadius = result.MaximumRootBracket * .5 + coordinateReserve;
  }
  return result;
}
} // namespace

static BRep_ImplicitPlaneSection::Result tracePlane(const TopoDS_Face& face,
  const BRep_ImplicitPlaneSection::Settings& settings,
  BRep_ImplicitPlaneSection::GraphResult* graph = nullptr)
{
  const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  if (data.IsNull() || !data->ExcludedDomains().empty() || typeid(*data) != typeid(BRep_TImplicitFace)
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
  if (settings.AllowExactLatticeNodes)
  {
    for (int row = 1; row <= 3; ++row)
      for (int column = 1; column <= 4; ++column)
        if (transform.Value(row, column) != (row == column ? 1. : 0.))
          throw Standard_NotImplemented("Exact lattice section nodes require identity field placement");
    for (const auto& direction : {settings.Frame.XDirection(), settings.Frame.YDirection()})
      for (int axis = 1; axis <= 3; ++axis)
        if (direction.Coord(axis) != 0 && std::abs(direction.Coord(axis)) != 1)
          throw Standard_NotImplemented("Exact lattice section nodes require a coordinate-plane frame");
  }
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
  BRep_ImplicitPlaneSection::Result result;
  for (int n = settings.MinimumCells; n <= settings.MaximumCells; n *= 2)
  {
    result = pass(face,
                  settings,
                  n,
                  *field,
                  lipschitz,
                  hessian,
                  fieldReserve,
                  gradientReserve,
                  coordinateReserve, graph);
    if (result.Resolved)
      return result;
    if (settings.AllowExactLatticeNodes && n - 1 >= settings.MinimumCells)
    {
      // Odd grids avoid repeated unrecognized eighth-period samples without
      // moving the rectangle, changing the field or enlarging the axis budget.
      auto candidate = pass(face, settings, n - 1, *field, lipschitz, hessian,
                            fieldReserve, gradientReserve, coordinateReserve, graph);
      if (candidate.Resolved)
        return candidate;
    }
  }
  // A dyadic sequence can repeatedly place a regular arc almost tangent to a
  // grid edge. Try bounded alternative uniform grids, with identical acceptance
  // proofs and no increase of the per-axis budget. Never perturb the geometry
  // or reinterpret unresolved cells as empty to force a successful result.
  for (int numerator : {6, 5, 7})
  {
    const int n = settings.MaximumCells * numerator / 8;
    if (n < settings.MinimumCells)
      continue;
    auto candidate = pass(face,
                          settings,
                          n,
                          *field,
                          lipschitz,
                          hessian,
                          fieldReserve,
                          gradientReserve,
                          coordinateReserve, graph);
    if (candidate.Resolved)
      return candidate;
  }
  return result;
}

BRep_ImplicitPlaneSection::Result BRep_ImplicitPlaneSection::Trace(const TopoDS_Face& face,
  const Settings& settings)
{
  return tracePlane(face, settings);
}

BRep_ImplicitPlaneSection::GraphResult BRep_ImplicitPlaneSection::TraceGraph(
  const TopoDS_Face& face, const Settings& settings)
{
  GraphResult result;
  result.Section = tracePlane(face, settings, &result);
  if (!result.Section.Resolved)
  { result.Vertices.clear(); result.Arcs.clear(); result.VertexRadius = 0; }
  return result;
}
