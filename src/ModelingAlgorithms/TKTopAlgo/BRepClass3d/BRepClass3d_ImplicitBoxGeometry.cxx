// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitBox.hxx>
#include <BRepClass3d_ImplicitBoxIO.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom_Plane.hxx>
#include <ElSLib.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <typeinfo>

namespace
{
using Box    = std::array<gp_Pnt, 2>;
using Report = BRepClass3d_ImplicitBox::CapGeometry;

void require(bool condition, const char* reason)
{
  if (!condition)
    throw Standard_NotImplemented(reason);
}

struct Polygon
{
  struct Node
  {
    Box    Bounds;
    size_t Begin = 0, End = 0, Left = 0, Right = 0;
  };

  gp_Pln                                 Plane;
  int                                    Side;
  Box                                    Bounds;
  std::vector<std::pair<gp_Pnt, gp_Pnt>> Segments;
  std::vector<Node>                      Nodes;

  size_t index(size_t begin, size_t end)
  {
    Node node;
    node.Begin  = begin;
    node.End    = end;
    node.Bounds = {Segments[begin].first, Segments[begin].first};
    for (size_t i = begin; i < end; ++i)
      for (const auto& p : {Segments[i].first, Segments[i].second})
        for (int a = 1; a <= 3; ++a)
        {
          node.Bounds[0].SetCoord(a, std::min(node.Bounds[0].Coord(a), p.Coord(a)));
          node.Bounds[1].SetCoord(a, std::max(node.Bounds[1].Coord(a), p.Coord(a)));
        }
    const size_t id = Nodes.size();
    Nodes.push_back(node);
    if (end - begin > 8)
    {
      int       axis  = (Side / 2 + 1) % 3 + 1;
      const int other = (Side / 2 + 2) % 3 + 1;
      if (node.Bounds[1].Coord(other) - node.Bounds[0].Coord(other)
          > node.Bounds[1].Coord(axis) - node.Bounds[0].Coord(axis))
        axis = other;
      const size_t middle = begin + (end - begin) / 2;
      std::nth_element(Segments.begin() + begin,
                       Segments.begin() + middle,
                       Segments.begin() + end,
                       [axis](const auto& a, const auto& b) {
                         return .5 * a.first.Coord(axis) + .5 * a.second.Coord(axis)
                                < .5 * b.first.Coord(axis) + .5 * b.second.Coord(axis);
                       });
      const size_t left = index(begin, middle), right = index(middle, end);
      Nodes[id].Left  = left;
      Nodes[id].Right = right;
    }
    return id;
  }

  gp_Pnt lift(double u, double v) const
  {
    const int  a = Side / 2, b = (a + 1) % 3, c = (a + 2) % 3;
    auto       p = Plane.Location();
    const auto n = Plane.Axis().Direction();
    p.SetCoord(b + 1, u);
    p.SetCoord(c + 1, v);
    p.SetCoord(a + 1,
               Plane.Location().Coord(a + 1)
                 - (n.Coord(b + 1) * (u - Plane.Location().Coord(b + 1))
                    + n.Coord(c + 1) * (v - Plane.Location().Coord(c + 1)))
                     / n.Coord(a + 1));
    return p;
  }
};

struct Mask
{
  bool   Inside = false;
  double Lower = 0, Upper = std::numeric_limits<double>::infinity();
};

struct Checker
{
  const Geom_ImplicitGyroid& Field;
  gp_Pnt                     Low, High;
  std::vector<Polygon>       Caps;
  Report&                    Result;
  double                     Tolerance, Reserve, FieldReserve, GradientReserve, H, L, Direction;
  int                        MaximumDepth;
  size_t                     MaximumCells, MaximumSegments;

  double variation(const gp_Pnt& p, double r, int axis) const
  {
    auto g = Field.Gradient(p);
    g.SetCoord(axis + 1, 0);
    return std::min(L * r, (g.Magnitude() + GradientReserve) * r + .5 * H * r * r)
             * (1 + 128 * std::numeric_limits<double>::epsilon())
           + FieldReserve;
  }

  Mask mask(const Polygon& cap, const gp_Pnt& p)
  {
    const int u = (cap.Side / 2 + 1) % 3 + 1, v = (cap.Side / 2 + 2) % 3 + 1;
    bool      inside   = false;
    double    distance = std::numeric_limits<double>::infinity(), upper = distance;
    // Each skipped node must be irrelevant to BOTH nearest-distance bounds
    // and the half-open parity ray. A distance-only BVH would corrupt holes.
    auto lower = [&](const Box& bounds, bool spatial) {
      double squared = 0;
      for (int a = 1; a <= 3; ++a)
        if (spatial || a == u || a == v)
        {
          const double delta =
            std::max({bounds[0].Coord(a) - p.Coord(a), p.Coord(a) - bounds[1].Coord(a), 0.});
          squared += delta * delta;
        }
      return std::max(0., std::sqrt(squared) - 8 * Reserve);
    };
    auto visit = [&](auto&& self, size_t id) -> void {
      const auto& node = cap.Nodes[id];
      const bool  ray  = p.Coord(v) >= node.Bounds[0].Coord(v) - 8 * Reserve
                         && p.Coord(v) <= node.Bounds[1].Coord(v) + 8 * Reserve
                         && p.Coord(u) <= node.Bounds[1].Coord(u) + 8 * Reserve;
      if (!ray && lower(node.Bounds, false) > distance && lower(node.Bounds, true) > upper)
        return;
      if (node.Left != 0)
      {
        // Near-first traversal tightens the upper bounds before pruning.
        const bool leftFirst =
          lower(cap.Nodes[node.Left].Bounds, false) <= lower(cap.Nodes[node.Right].Bounds, false);
        self(self, leftFirst ? node.Left : node.Right);
        self(self, leftFirst ? node.Right : node.Left);
        return;
      }
      for (size_t i = node.Begin; i < node.End; ++i)
      {
        const auto& edge = cap.Segments[i];
        require(Result.SegmentTests < MaximumSegments, "Cap geometry segment budget exhausted");
        ++Result.SegmentTests;
        const auto&  a  = edge.first;
        const auto&  b  = edge.second;
        const double du = b.Coord(u) - a.Coord(u), dv = b.Coord(v) - a.Coord(v),
                     length = du * du + dv * dv;
        const double t =
          length > 0
            ? std::clamp(((p.Coord(u) - a.Coord(u)) * du + (p.Coord(v) - a.Coord(v)) * dv) / length,
                         0.,
                         1.)
            : 0;
        distance =
          std::min(distance,
                   std::hypot(p.Coord(u) - a.Coord(u) - t * du, p.Coord(v) - a.Coord(v) - t * dv));
        const gp_Vec ab(a, b);
        const double s = ab.SquareMagnitude() > 0
                           ? std::clamp(gp_Vec(a, p).Dot(ab) / ab.SquareMagnitude(), 0., 1.)
                           : 0;
        upper          = std::min(upper, p.Distance(a.Translated(ab * s)) + 8 * Reserve);
        if ((a.Coord(v) > p.Coord(v)) != (b.Coord(v) > p.Coord(v)))
        {
          const double x = a.Coord(u) + (p.Coord(v) - a.Coord(v)) * du / dv;
          if (p.Coord(u) < x)
            inside = !inside;
        }
      }
    };
    visit(visit, 0);
    Mask result;
    result.Lower  = std::max(0., distance - 8 * Reserve);
    result.Inside = inside && result.Lower > 8 * Reserve;
    result.Upper =
      result.Inside ? p.Distance(cap.lift(p.Coord(u), p.Coord(v))) + 8 * Reserve : upper;
    return result;
  }

  gp_Pnt projected(const gp_Pnt& p, int side) const
  {
    auto q = p;
    for (int a = 1; a <= 3; ++a)
      q.SetCoord(a, std::clamp(q.Coord(a), Low.Coord(a), High.Coord(a)));
    q.SetCoord(side / 2 + 1, side % 2 ? High.Coord(side / 2 + 1) : Low.Coord(side / 2 + 1));
    return q;
  }

  double witness(const gp_Pnt& p, int side, double available) const
  {
    if (!(available > 16 * Reserve))
      return std::numeric_limits<double>::infinity();
    const auto base = projected(p, side);
    double     best = std::numeric_limits<double>::infinity();
    auto       test = [&](const gp_Pnt& q) {
      const double distance = p.Distance(q) + 8 * Reserve;
      if (distance <= available && Direction * Field.Value(q) < -FieldReserve)
        best = std::min(best, distance);
    };
    test(base);
    auto gradient = Field.Gradient(base);
    gradient.SetCoord(side / 2 + 1, 0);
    std::vector<gp_Vec> directions;
    if (gradient.Magnitude() > GradientReserve)
      directions.push_back(gradient / gradient.Magnitude());
    for (int a = 0; a < 3; ++a)
      if (a != side / 2)
      {
        gp_Vec d(0, 0, 0);
        d.SetCoord(a + 1, 1);
        directions.push_back(d);
      }
    for (const auto& d : directions)
      for (double factor : {-.25, .25, -.5, .5, -.9, .9})
        test(projected(base.Translated(d * (available * factor)), side));
    return best;
  }

  bool visit(const Box& box, int side, const Polygon* native, int depth)
  {
    const int direction = native ? 1 : 0, axis = side / 2, u = (axis + 1) % 3, v = (axis + 2) % 3;
    Result.UnresolvedSide      = side;
    Result.UnresolvedDirection = direction;
    Result.UnresolvedMinimum   = box[0];
    Result.UnresolvedMaximum   = box[1];
    if (Result.Cells[0] + Result.Cells[1] >= MaximumCells)
    {
      Result.Reason = "Cap geometry cell budget exhausted";
      return false;
    }
    ++Result.Cells[size_t(direction)];
    const double cu = (box[0].Coord(u + 1) + box[1].Coord(u + 1)) * .5,
                 cv = (box[0].Coord(v + 1) + box[1].Coord(v + 1)) * .5;
    gp_Pnt       p;
    p = native ? native->lift(cu, cv) : projected(gp_Pnt((box[0].XYZ() + box[1].XYZ()) * .5), side);
    double r = 0, projection = 0, nativePlaneDistance = 0;
    for (int a : {0, 1})
      for (int b : {0, 1})
      {
        const double x = box[size_t(a)].Coord(u + 1), y = box[size_t(b)].Coord(v + 1);
        gp_Pnt       q = native ? native->lift(x, y) : p;
        if (!native)
        {
          q.SetCoord(u + 1, x);
          q.SetCoord(v + 1, y);
        }
        r          = std::max(r, p.Distance(q));
        projection = std::max(projection, q.Distance(projected(q, side)));
      }
    r += 8 * Reserve;
    if (native)
    {
      const auto   state        = mask(*native, p);
      const double planarRadius = std::hypot(box[1].Coord(u + 1) - box[0].Coord(u + 1),
                                             box[1].Coord(v + 1) - box[0].Coord(v + 1))
                                    * .5
                                  + 8 * Reserve;
      if (!state.Inside && state.Lower > planarRadius)
        return true;
      const auto base = projected(p, side);
      if (Direction * Field.Value(base) + variation(base, r, axis) < 0
          && projection + 16 * Reserve <= Tolerance)
      {
        Result.CapsToMathematicalBound =
          std::max(Result.CapsToMathematicalBound, projection + 16 * Reserve);
        return true;
      }
      const double upper = witness(p, side, Tolerance - r - 16 * Reserve) + r + 16 * Reserve;
      if (upper <= Tolerance)
      {
        Result.CapsToMathematicalBound = std::max(Result.CapsToMathematicalBound, upper);
        return true;
      }
    }
    else
    {
      if (Direction * Field.Value(p) > variation(p, r, axis))
        return true;
      double upper = std::numeric_limits<double>::infinity();
      for (const auto& cap : Caps)
        if (cap.Side == side)
        {
          const auto state = mask(cap, p);
          upper            = std::min(upper, state.Upper + r + 16 * Reserve);
          if (state.Inside && state.Lower > r)
          {
            nativePlaneDistance = 0;
            for (int a : {0, 1})
              for (int b : {0, 1})
              {
                auto q = p;
                q.SetCoord(u + 1, box[size_t(a)].Coord(u + 1));
                q.SetCoord(v + 1, box[size_t(b)].Coord(v + 1));
                nativePlaneDistance =
                  std::max(nativePlaneDistance,
                           q.Distance(cap.lift(q.Coord(u + 1), q.Coord(v + 1))));
              }
            upper = std::min(upper, nativePlaneDistance + 16 * Reserve);
          }
        }
      if (upper <= Tolerance)
      {
        Result.MathematicalToCapsBound = std::max(Result.MathematicalToCapsBound, upper);
        return true;
      }
    }
    if (depth >= MaximumDepth)
    {
      Result.Reason = "Cap geometry depth exhausted";
      return false;
    }
    const int split =
      box[1].Coord(u + 1) - box[0].Coord(u + 1) >= box[1].Coord(v + 1) - box[0].Coord(v + 1) ? u
                                                                                             : v;
    const double middle = (box[0].Coord(split + 1) + box[1].Coord(split + 1)) * .5;
    if (!(middle > box[0].Coord(split + 1) && middle < box[1].Coord(split + 1)))
    {
      Result.Reason = "Cap geometry coordinate precision exhausted";
      return false;
    }
    auto left = box, right = box;
    left[1].SetCoord(split + 1, middle);
    right[0].SetCoord(split + 1, middle);
    return visit(left, side, native, depth + 1) && visit(right, side, native, depth + 1);
  }
};
} // namespace

BRepClass3d_ImplicitBox::CapGeometry BRepClass3d_ImplicitBox::CheckCapGeometry(
  const Boundary&                 source,
  const std::vector<TopoDS_Face>& candidates,
  double                          tolerance,
  int                             maximumDepth,
  size_t                          maximumCells,
  size_t                          maximumSegments)
{
  require(source.Policy == CornerPolicy::Reject,
          "Corner-aware cap geometry validation is not implemented");
  if (!std::isfinite(tolerance) || tolerance <= 0 || maximumDepth < 1 || maximumDepth > 24
      || !maximumCells || maximumCells > 10000000 || !maximumSegments || maximumSegments > 50000000
      || candidates.size() > 128)
    throw Standard_DomainError("Invalid native cap geometry settings");
  std::ostringstream record;
  BRepClass3d_ImplicitBoxIO::Write(source, record);
  std::istringstream input(record.str());
  Boundary           copy;
  BRepClass3d_ImplicitBoxIO::Read(input, copy);
  const auto   data  = occ::down_cast<BRep_TImplicitFace>(copy.Faces.back().TShape());
  const auto&  field = *static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  const auto   low = data->Minimum(), high = data->Maximum();
  const double k = 2 * std::acos(-1.) / field.Period(), h = 4 * k * k, l = field.LipschitzBound();
  double       magnitude = 1;
  for (int a = 1; a <= 3; ++a)
  {
    magnitude = std::max({magnitude, std::abs(low.Coord(a)), std::abs(high.Coord(a))});
    require(std::abs(k * low.Coord(a)) < 1.e6 && std::abs(k * high.Coord(a)) < 1.e6,
            "Cap geometry trigonometric precision exhausted");
  }
  Report               report;
  std::vector<Polygon> polygons;
  size_t               poles = 0;
  for (const auto& candidate : candidates)
  {
    require(!candidate.IsNull() && !BRep_Tool::IsImplicit(candidate)
              && candidate.Location().IsIdentity(),
            "Expected conventional planar cap");
    const auto face = TopoDS::Face(BRepBuilderAPI_Copy(candidate, true, false).Shape());
    require(BRepCheck_Analyzer(face).IsValid(), "Invalid conventional cap geometry");
    TopLoc_Location surfaceLocation;
    const auto      surfaceGeometry = BRep_Tool::Surface(face, surfaceLocation);
    require(surfaceLocation.IsIdentity() && !surfaceGeometry.IsNull()
              && typeid(*surfaceGeometry) == typeid(Geom_Plane),
            "Unsupported cap surface representation/location");
    BRepAdaptor_Surface surface(face);
    require(surface.GetType() == GeomAbs_Plane, "Cap geometry requires planes");
    Polygon cap;
    cap.Plane = surface.Plane();
    for (int d = 1; d <= 3; ++d)
      magnitude = std::max(magnitude, std::abs(cap.Plane.Location().Coord(d)));
    const auto normal = cap.Plane.Axis().Direction();
    int        axis   = 0;
    for (int a = 1; a < 3; ++a)
      if (std::abs(normal.Coord(a + 1)) > std::abs(normal.Coord(axis + 1)))
        axis = a;
    require(std::abs(normal.Coord(axis + 1)) > 1 - 1.e-12, "Cap plane is not box-aligned");
    cap.Side              = 2 * axis + int(normal.Coord(axis + 1) > 0);
    const double position = cap.Side % 2 ? high.Coord(axis + 1) : low.Coord(axis + 1);
    require(std::abs(cap.Plane.Location().Coord(axis + 1) - position)
              <= std::max(1.e-7, 8192 * std::numeric_limits<double>::epsilon() * magnitude),
            "Cap plane is not a box side");
    bool firstPoint = true;
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
    {
      const auto      edge = TopoDS::Edge(it.Current());
      TopLoc_Location location;
      double          a, b, ua, ub;
      const auto curve = occ::down_cast<Geom_BSplineCurve>(BRep_Tool::Curve(edge, location, a, b));
      const auto uv =
        occ::down_cast<Geom2d_BSplineCurve>(BRep_Tool::CurveOnSurface(edge, face, ua, ub));
      require(location.IsIdentity() && !curve.IsNull() && !uv.IsNull() && curve->Degree() == 1
                && uv->Degree() == 1 && typeid(*curve) == typeid(Geom_BSplineCurve)
                && typeid(*uv) == typeid(Geom2d_BSplineCurve) && !curve->IsRational()
                && !uv->IsRational() && !curve->IsPeriodic() && !uv->IsPeriodic()
                && curve->NbPoles() == uv->NbPoles() && curve->NbKnots() == curve->NbPoles()
                && uv->NbKnots() == curve->NbKnots() && a == ua && b == ub
                && a == curve->FirstParameter() && b == curve->LastParameter()
                && ua == uv->FirstParameter() && ub == uv->LastParameter(),
              "Unsupported cap curve parameterization");
      require((poles += size_t(uv->NbPoles())) <= 1000000,
              "Cap geometry input pole budget exhausted");
      gp_Pnt previous;
      for (int i = 1; i <= uv->NbPoles(); ++i)
      {
        require(curve->Knot(i) == uv->Knot(i) && curve->Multiplicity(i) == uv->Multiplicity(i)
                  && curve->Multiplicity(i) == ((i == 1 || i == curve->NbPoles()) ? 2 : 1),
                "Cap curve parameters differ");
        const auto q = uv->Pole(i);
        const auto p = ElSLib::Value(q.X(), q.Y(), cap.Plane);
        require(std::isfinite(q.X()) && std::isfinite(q.Y()), "Nonfinite cap parameter pole");
        magnitude = std::max({magnitude, std::abs(q.X()), std::abs(q.Y())});
        for (int d = 1; d <= 3; ++d)
        {
          require(std::isfinite(p.Coord(d)) && std::isfinite(curve->Pole(i).Coord(d)),
                  "Nonfinite cap pole");
          magnitude =
            std::max({magnitude, std::abs(p.Coord(d)), std::abs(curve->Pole(i).Coord(d))});
          if (firstPoint)
            cap.Bounds = {p, p};
          cap.Bounds[0].SetCoord(d, std::min(cap.Bounds[0].Coord(d), p.Coord(d)));
          cap.Bounds[1].SetCoord(d, std::max(cap.Bounds[1].Coord(d), p.Coord(d)));
        }
        firstPoint = false;
        report.EdgeVertexDeviation =
          std::max(report.EdgeVertexDeviation, p.Distance(curve->Pole(i)));
        if (i > 1)
          cap.Segments.emplace_back(previous, p);
        previous = p;
      }
      TopoDS_Vertex first, last;
      TopExp::Vertices(TopoDS::Edge(edge.Oriented(TopAbs_FORWARD)), first, last);
      for (const auto& vertex : {first, last})
        for (int d = 1; d <= 3; ++d)
        {
          const double coordinate = BRep_Tool::Pnt(vertex).Coord(d);
          require(std::isfinite(coordinate), "Nonfinite cap vertex");
          magnitude = std::max(magnitude, std::abs(coordinate));
        }
      report.EdgeVertexDeviation =
        std::max({report.EdgeVertexDeviation,
                  BRep_Tool::Pnt(first).Distance(curve->Pole(1)),
                  BRep_Tool::Pnt(last).Distance(curve->Pole(curve->NbPoles()))});
    }
    require(!cap.Segments.empty(), "Cap has no polygon boundary");
    cap.index(0, cap.Segments.size());
    polygons.push_back(std::move(cap));
  }
  const double reserve = 8192 * std::numeric_limits<double>::epsilon() * magnitude;
  if (magnitude > 1.e50 || !std::isfinite(h) || h <= 0 || !std::isfinite(l)
      || tolerance <= 4096 * reserve)
    throw Standard_DomainError("Cap geometry coordinate/derivative precision exhausted");
  // Two possible deviations add: vertex->3D edge and 3D edge->lifted p-curve.
  report.EdgeVertexDeviation = 2 * report.EdgeVertexDeviation + 16 * reserve;
  const double available     = tolerance - report.EdgeVertexDeviation - 32 * reserve;
  if (!(available > 4096 * reserve))
  {
    report.Reason = "Edge/vertex deviation exhausts geometry tolerance";
    return report;
  }
  Checker checker{field,
                  low,
                  high,
                  std::move(polygons),
                  report,
                  available,
                  reserve,
                  4096 * std::numeric_limits<double>::epsilon() * (1 + std::abs(field.IsoValue()))
                    + l * reserve,
                  h * reserve + 2048 * std::numeric_limits<double>::epsilon() * l,
                  h,
                  l,
                  copy.Faces.back().Orientation() == TopAbs_FORWARD ? 1. : -1.,
                  maximumDepth,
                  maximumCells,
                  maximumSegments};
  try
  {
    for (int side = 0; side < 6; ++side)
    {
      Box          domain{low, high};
      const int    axis     = side / 2 + 1;
      const double position = side % 2 ? high.Coord(axis) : low.Coord(axis);
      domain[0].SetCoord(axis, position);
      domain[1].SetCoord(axis, position);
      if (!checker.visit(domain, side, nullptr, 0))
        return report;
    }
    for (const auto& cap : checker.Caps)
      if (!checker.visit(cap.Bounds, cap.Side, &cap, 0))
        return report;
    report.BoundaryBound = std::max(report.MathematicalToCapsBound,
                                    report.CapsToMathematicalBound + report.EdgeVertexDeviation)
                           + 16 * reserve;
    report.Resolved      = report.BoundaryBound <= tolerance;
    report.Reason =
      report.Resolved ? "native_cap_geometry_resolved" : "Native boundary tolerance exhausted";
    report.UnresolvedSide = report.UnresolvedDirection = -1;
  }
  catch (const Standard_NotImplemented& error)
  {
    report.Reason = error.what();
  }
  return report;
}

BRepClass3d_ImplicitBox::NativeSolids BRepClass3d_ImplicitBox::BuildNativeSolids(
  const Boundary& source,
  double          tolerance,
  int             resolution,
  int             maximumDepth,
  size_t          maximumCells,
  size_t          maximumTests,
  size_t          maximumSeparatorCells,
  size_t          maximumGeometryCells,
  size_t          maximumSegments)
{
  if (!std::isfinite(tolerance) || tolerance <= 0 || !maximumGeometryCells
      || maximumGeometryCells > 10000000 || !maximumSegments || maximumSegments > 50000000)
    throw Standard_DomainError("Invalid native solid geometry settings");
  NativeSolids result;
  result.Tolerance = tolerance;
  result.Boundary  = BuildNativeShells(source,
                                       resolution,
                                       maximumDepth,
                                       maximumCells,
                                       maximumTests,
                                       maximumSeparatorCells);
  if (!result.Boundary.ShellsConstructed || !result.Boundary.TopologyValidated)
  {
    result.Reason = "native_shells_unresolved: " + result.Boundary.Reason;
    return result;
  }
  std::vector<TopoDS_Face> caps;
  for (const auto& face : result.Boundary.Faces)
  {
    result.MaximumStoredTolerance =
      std::max(result.MaximumStoredTolerance, BRep_Tool::Tolerance(face));
    if (!BRep_Tool::IsImplicit(face))
      caps.push_back(face);
  }
  for (TopExp_Explorer it(result.Boundary.Shape, TopAbs_EDGE); it.More(); it.Next())
    result.MaximumStoredTolerance =
      std::max(result.MaximumStoredTolerance, BRep_Tool::Tolerance(TopoDS::Edge(it.Current())));
  for (TopExp_Explorer it(result.Boundary.Shape, TopAbs_VERTEX); it.More(); it.Next())
    result.MaximumStoredTolerance =
      std::max(result.MaximumStoredTolerance, BRep_Tool::Tolerance(TopoDS::Vertex(it.Current())));
  if (result.MaximumStoredTolerance > tolerance)
  {
    result.Reason = "Stored native tolerances exceed requested solid tolerance";
    return result;
  }
  result.Geometry = CheckCapGeometry(result.Boundary.Source,
                                     caps,
                                     tolerance,
                                     maximumDepth,
                                     maximumGeometryCells,
                                     maximumSegments);
  if (!result.Geometry.Resolved)
  {
    result.Reason = "native_geometry_unresolved: " + result.Geometry.Reason;
    return result;
  }
  int regions = 0;
  for (int region : result.Boundary.ShellRegions)
    regions = std::max(regions, region + 1);
  BRep_Builder              builder;
  std::vector<TopoDS_Solid> solids;
  for (int region = 0; region < regions; ++region)
  {
    TopoDS_Solid solid;
    builder.MakeSolid(solid);
    bool present = false;
    for (size_t i = 0; i < result.Boundary.Shells.size(); ++i)
      if (result.Boundary.ShellRegions[i] == region)
      {
        builder.Add(solid, result.Boundary.Shells[i]);
        present = true;
      }
    if (!present)
    {
      result.Reason = "Native material region has no shells";
      return result;
    }
    solid.Closed(true);
    solids.push_back(solid);
  }
  builder.MakeCompound(result.Shape);
  for (const auto& solid : solids)
    builder.Add(result.Shape, solid);
  result.Solids      = std::move(solids);
  result.Constructed = result.Validated = true;
  result.Reason                         = "native_solids_validated";
  return result;
}
