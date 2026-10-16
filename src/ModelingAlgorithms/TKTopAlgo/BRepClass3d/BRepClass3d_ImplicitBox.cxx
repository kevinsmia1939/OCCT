// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitBox.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom_Plane.hxx>
#include <ElSLib.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <typeinfo>

namespace
{
using Map = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

void require(bool condition, const char* message)
{
  if (!condition)
    throw Standard_NotImplemented(message);
}

struct Node
{
  gp_Pnt        Point;
  TopoDS_Vertex Vertex;
  int           RootEdge = -1;
};

struct Arc
{
  int                 First, Last;
  std::vector<gp_Pnt> Points;
  TopoDS_Edge         Edge;
};

struct Use
{
  int  Index;
  bool Forward;
};

struct Loop
{
  TopoDS_Wire           Wire;
  std::vector<Use>      Uses;
  std::vector<gp_Pnt2d> Polygon;
  double                Area = 0;
};

struct BoxEdge
{
  int              First, Last;
  std::vector<int> Roots, Arcs;
};

struct Side
{
  int                                 Axis, Sign, Source;
  gp_Pln                              Plane;
  BRep_ImplicitPlaneSection::Settings Settings;
  std::vector<Use>                    Uses;
};

gp_Vec unit(int axis, int sign = 1)
{
  gp_Vec v(0, 0, 0);
  v.SetCoord(axis + 1, sign);
  return v;
}

int addNode(std::vector<Node>& nodes, const gp_Pnt& point, double precision, int rootEdge = -1)
{
  Node node;
  node.Point    = point;
  node.RootEdge = rootEdge;
  BRep_Builder().MakeVertex(node.Vertex, point, precision);
  nodes.push_back(node);
  return int(nodes.size()) - 1;
}

// Recursively exclude empty intervals using derivative bounds. An interval
// contributes a root only after uniqueness is established by monotonicity.
// Tangencies, near-zero subdivision nodes and resource exhaustion reject.
void roots(const Geom_ImplicitGyroid& field,
           const gp_Pnt&              first,
           const gp_Pnt&              last,
           double                     hessian,
           double                     bracket,
           double                     reserve,
           std::vector<gp_Pnt>&       out,
           int&                       budget,
           int                        depth = 0)
{
  require(--budget >= 0 && depth <= 32, "Box-edge root isolation budget exhausted");
  const gp_Pnt mid((first.XYZ() + last.XYZ()) * .5);
  const double length = first.Distance(last), fm = field.Value(mid);
  require(std::isfinite(length) && length > 0, "Box-edge coordinate collapse");
  if (std::abs(fm) > field.LipschitzBound() * length * .5 + reserve)
    return;
  const double fa = field.Value(first), fb = field.Value(last);
  require(std::abs(fa) > reserve && std::abs(fb) > reserve,
          "Box corner/subdivision node lies on the field; ambiguous edge topology");
  const double derivative = field.Gradient(mid).Dot(gp_Vec(first, last) / length);
  const bool   monotone   = std::abs(derivative) > hessian * length * .5 + reserve;
  if (monotone)
  {
    if (std::signbit(fa) == std::signbit(fb))
      return;
    gp_Pnt a = first, b = last;
    double va = fa;
    for (int n = 0; n < 100 && a.Distance(b) > bracket; ++n)
    {
      const gp_Pnt p((a.XYZ() + b.XYZ()) * .5);
      const double value = field.Value(p);
      if (value == 0)
      {
        a = b = p;
        break;
      }
      if (std::signbit(va) != std::signbit(value))
        b = p;
      else
      {
        a  = p;
        va = value;
      }
    }
    require(a.Distance(b) <= bracket, "Box-edge root failed its distance bracket");
    out.emplace_back((a.XYZ() + b.XYZ()) * .5);
    return;
  }
  require(length > bracket, "Tangent/critical box-edge root unresolved");
  roots(field, first, mid, hessian, bracket, reserve, out, budget, depth + 1);
  roots(field, mid, last, hessian, bracket, reserve, out, budget, depth + 1);
}

int addArc(std::vector<Arc>&        arcs,
           const std::vector<Node>& nodes,
           int                      first,
           int                      last,
           std::vector<gp_Pnt>      points,
           double                   tolerance)
{
  require(points.size() >= 2, "Empty trim arc");
  points.front()                   = nodes[size_t(first)].Point;
  points.back()                    = nodes[size_t(last)].Point;
  const int                  count = int(points.size());
  NCollection_Array1<gp_Pnt> poles(1, count);
  NCollection_Array1<double> knots(1, count);
  NCollection_Array1<int>    multiplicities(1, count);
  for (int i = 1; i <= count; ++i)
  {
    poles(i)          = points[size_t(i - 1)];
    knots(i)          = i - 1;
    multiplicities(i) = i == 1 || i == count ? 2 : 1;
    if (i > 1)
      require(poles(i).Distance(poles(i - 1)) > 0, "Collapsed trim segment");
  }
  const occ::handle<Geom_BSplineCurve> curve =
    new Geom_BSplineCurve(poles, knots, multiplicities, 1);
  Arc          arc{first, last, std::move(points), {}};
  BRep_Builder builder;
  builder.MakeEdge(arc.Edge, curve, tolerance);
  builder.Add(arc.Edge, nodes[size_t(first)].Vertex.Oriented(TopAbs_FORWARD));
  builder.Add(arc.Edge, nodes[size_t(last)].Vertex.Oriented(TopAbs_REVERSED));
  builder.Range(arc.Edge, 0, count - 1);
  builder.UpdateVertex(TopoDS::Vertex(nodes[size_t(first)].Vertex.Oriented(TopAbs_FORWARD)),
                       0,
                       arc.Edge,
                       BRep_Tool::Tolerance(nodes[size_t(first)].Vertex));
  builder.UpdateVertex(TopoDS::Vertex(nodes[size_t(last)].Vertex.Oriented(TopAbs_REVERSED)),
                       count - 1,
                       arc.Edge,
                       BRep_Tool::Tolerance(nodes[size_t(last)].Vertex));
  arc.Edge.Closed(first == last);
  arcs.push_back(std::move(arc));
  return int(arcs.size()) - 1;
}

std::vector<Loop> assemble(const std::vector<Use>& uses,
                           const std::vector<Arc>& arcs,
                           const gp_Pln*           plane)
{
  std::map<int, size_t> outgoing;
  std::map<int, int>    incoming;
  for (size_t i = 0; i < uses.size(); ++i)
  {
    const auto& use   = uses[i];
    const auto& arc   = arcs[size_t(use.Index)];
    const int   first = use.Forward ? arc.First : arc.Last;
    const int   last  = use.Forward ? arc.Last : arc.First;
    require(outgoing.emplace(first, i).second && ++incoming[last] == 1,
            "Branching trim graph; no boundary published");
  }
  for (const auto& entry : outgoing)
    require(incoming[entry.first] == 1, "Open trim graph; no boundary published");
  std::vector<bool> visited(uses.size(), false);
  std::vector<Loop> result;
  BRep_Builder      builder;
  for (size_t i = 0; i < uses.size(); ++i)
  {
    if (visited[i])
      continue;
    Loop loop;
    builder.MakeWire(loop.Wire);
    size_t current = i;
    do
    {
      require(!visited[current], "Trim traversal entered another loop");
      visited[current] = true;
      const auto  use  = uses[current];
      const auto& arc  = arcs[size_t(use.Index)];
      builder.Add(loop.Wire, arc.Edge.Oriented(use.Forward ? TopAbs_FORWARD : TopAbs_REVERSED));
      loop.Uses.push_back(use);
      if (plane)
        for (size_t p = 0; p + 1 < arc.Points.size(); ++p)
        {
          const auto& point = arc.Points[use.Forward ? p : arc.Points.size() - 1 - p];
          double      u, v;
          ElSLib::Parameters(*plane, point, u, v);
          loop.Polygon.emplace_back(u, v);
        }
      const int  end  = use.Forward ? arc.Last : arc.First;
      const auto next = outgoing.find(end);
      require(next != outgoing.end(), "Missing trim continuation");
      current = next->second;
    } while (current != i);
    loop.Wire.Closed(true);
    for (size_t p = 0; p < loop.Polygon.size(); ++p)
    {
      const auto& a = loop.Polygon[p];
      const auto& b = loop.Polygon[(p + 1) % loop.Polygon.size()];
      loop.Area += (a.X() * b.Y() - a.Y() * b.X()) * .5;
    }
    result.push_back(std::move(loop));
  }
  return result;
}

} // namespace

BRepClass3d_ImplicitBox::PointClassification BRepClass3d_ImplicitBox::Classify(
  const Boundary& boundary,
  const gp_Pnt&   point,
  double          tolerance)
{
  if (!std::isfinite(tolerance) || tolerance <= 0)
    throw Standard_DomainError("Classification distance tolerance must be finite and positive");
  require(!boundary.Faces.empty(), "Classification requires a constructed native mixed boundary");
  const auto& face = boundary.Faces.back();
  const auto  data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  require(!data.IsNull() && typeid(*data) == typeid(BRep_TImplicitFace)
            && typeid(*data->Field()) == typeid(Geom_ImplicitGyroid) && face.Location().IsIdentity()
            && data->Location().IsIdentity()
            && (face.Orientation() == TopAbs_FORWARD || face.Orientation() == TopAbs_REVERSED),
          "Classification requires an unchanged identity-placed mixed gyroid/box boundary");
  const auto&  field = *static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  const auto   low = data->Minimum(), high = data->Maximum();
  double       magnitude = 1;
  const double k         = 2 * std::acos(-1.) / field.Period();
  for (int axis = 1; axis <= 3; ++axis)
  {
    if (!std::isfinite(point.Coord(axis)) || !std::isfinite(k * point.Coord(axis))
        || std::abs(k * point.Coord(axis)) >= 1.e6)
      throw Standard_DomainError(
        "Classification coordinate exceeds the finite trigonometric precision policy");
    magnitude = std::max({magnitude,
                          std::abs(point.Coord(axis)),
                          std::abs(low.Coord(axis)),
                          std::abs(high.Coord(axis))});
  }
  const double coordinateReserve = 512 * std::numeric_limits<double>::epsilon() * magnitude;
  if (tolerance <= 4096 * coordinateReserve)
    throw Standard_DomainError("Classification tolerance below coordinate precision policy");
  const double lipschitz = field.LipschitzBound();
  const double fieldReserve =
    512 * std::numeric_limits<double>::epsilon() * (1 + std::abs(field.IsoValue()))
    + lipschitz * coordinateReserve;
  const double value    = field.Value(point);
  const auto   gradient = field.Gradient(point);
  const double norm = std::hypot(gradient.X(), gradient.Y(), gradient.Z()), hessian = 4 * k * k;
  require(std::isfinite(norm) && std::isfinite(hessian) && hessian > 0,
          "Unrepresentable classification derivative bounds");
  const bool negative          = face.Orientation() == TopAbs_FORWARD;
  const bool retained          = (value < 0) == negative;
  const bool signResolved      = std::abs(value) > fieldReserve;
  gp_Pnt     clipped           = point;
  double     boxInsideDistance = std::numeric_limits<double>::infinity();
  bool       insideBox         = true;
  for (int axis = 1; axis <= 3; ++axis)
  {
    clipped.SetCoord(axis, std::clamp(point.Coord(axis), low.Coord(axis), high.Coord(axis)));
    insideBox &= point.Coord(axis) >= low.Coord(axis) && point.Coord(axis) <= high.Coord(axis);
    boxInsideDistance = std::min({boxInsideDistance,
                                  point.Coord(axis) - low.Coord(axis),
                                  high.Coord(axis) - point.Coord(axis)});
  }
  PointClassification result;
  const auto          distance = [](const gp_Pnt& a, const gp_Pnt& b) {
    return std::hypot(a.X() - b.X(), a.Y() - b.Y(), a.Z() - b.Z());
  };
  const double margin = std::max(0., std::abs(value) - fieldReserve);
  const double gradientReserve =
    hessian * coordinateReserve + 512 * std::numeric_limits<double>::epsilon() * lipschitz;
  const double derivative = norm + gradientReserve;
  // Solve margin = |grad|*r + H*r^2/2 without subtractive cancellation.
  // This is another conservative distance lower bound, useful near stationary
  // points; it is NOT an approximation to signed distance.
  const double taylorDistance =
    margin == 0
      ? 0
      : 2 * margin / (derivative + std::hypot(derivative, std::sqrt(2 * hessian * margin)));
  const double fieldDistance =
    std::max(0., std::max(margin / lipschitz, taylorDistance) - coordinateReserve);
  const double boxOutsideDistance = std::max(0., distance(point, clipped) - coordinateReserve);
  // The region is the intersection of the box and one field-sign region.
  // Inside: distance to either possible boundary is at least their minimum.
  // Outside: each violated constraint gives a lower bound; take the maximum.
  if (insideBox && signResolved && retained)
    result.DistanceLowerBound =
      std::max(0., std::min(boxInsideDistance - coordinateReserve, fieldDistance));
  else
    result.DistanceLowerBound =
      std::max(boxOutsideDistance, signResolved && !retained ? fieldDistance : 0.);
  if (result.DistanceLowerBound > tolerance + coordinateReserve)
  {
    result.State = insideBox && retained && signResolved ? PointState::Inside : PointState::Outside;
    return result;
  }
  // A retained point on a box side is an actual boundary witness even when a
  // gyroid root is absent nearby. Test all six clipped orthogonal cap points.
  for (int axis = 1; axis <= 3; ++axis)
    for (double coordinate : {low.Coord(axis), high.Coord(axis)})
    {
      auto witness = clipped;
      witness.SetCoord(axis, coordinate);
      const double witnessDistance = distance(point, witness) + coordinateReserve;
      if (witnessDistance > tolerance)
        continue;
      const double f = field.Value(witness);
      if (std::abs(f) <= fieldReserve || ((f < 0) != negative))
        continue;
      result.State              = PointState::NearBoundary;
      result.HasBoundaryWitness = true;
      result.WitnessFirst = result.WitnessLast = witness;
      result.WitnessDistanceUpperBound         = witnessDistance;
      return result;
    }
  // Search along the local gradient, clipping the entire search segment to
  // the box. A sign bracket proves existence by continuity, not proximity
  // inferred from a dimensionless residual or a Newton convergence flag.
  auto direction = gradient;
  if (!std::isfinite(norm) || norm <= gradientReserve)
    return result;
  direction /= norm;
  double first = -(tolerance - 4 * coordinateReserve), last = -first;
  for (int axis = 1; axis <= 3; ++axis)
  {
    const double component = direction.Coord(axis);
    if (component == 0)
    {
      if (point.Coord(axis) < low.Coord(axis) || point.Coord(axis) > high.Coord(axis))
        return result;
      continue;
    }
    const double a = (low.Coord(axis) - point.Coord(axis)) / component;
    const double b = (high.Coord(axis) - point.Coord(axis)) / component;
    first          = std::max(first, std::min(a, b));
    last           = std::min(last, std::max(a, b));
  }
  // Shrink away from box clipping roundoff rather than evaluating outside it.
  first += 2 * coordinateReserve;
  last -= 2 * coordinateReserve;
  if (first >= last)
    return result;
  auto   at = [&](double t) { return point.Translated(direction * t); };
  double fa = field.Value(at(first)), fb = field.Value(at(last));
  if (std::abs(fa) <= fieldReserve || std::abs(fb) <= fieldReserve
      || std::signbit(fa) == std::signbit(fb))
    return result;
  for (int iteration = 0; iteration < 100 && last - first > tolerance / 256; ++iteration)
  {
    const double mid = first + (last - first) * .5;
    if (mid == first || mid == last)
      break;
    const double fm = field.Value(at(mid));
    // Keep a valid bracket when the midpoint's sign is within roundoff.
    // Test quarter points; never collapse an uncertain midpoint to an exact root.
    if (std::abs(fm) <= fieldReserve)
    {
      const double a = first + (mid - first) * .5, b = mid + (last - mid) * .5;
      const double va = field.Value(at(a)), vb = field.Value(at(b));
      if (std::abs(va) <= fieldReserve || std::abs(vb) <= fieldReserve
          || std::signbit(va) == std::signbit(vb))
        break;
      first = a;
      last  = b;
      fa    = va;
    }
    else if (std::signbit(fa) != std::signbit(fm))
      last = mid;
    else
    {
      first = mid;
      fa    = fm;
    }
  }
  const auto a = at(first), b = at(last);
  for (int axis = 1; axis <= 3; ++axis)
    if (a.Coord(axis) < low.Coord(axis) || a.Coord(axis) > high.Coord(axis)
        || b.Coord(axis) < low.Coord(axis) || b.Coord(axis) > high.Coord(axis))
      return result;
  const double upper = std::max(distance(point, a), distance(point, b)) + coordinateReserve;
  if (upper > tolerance || distance(a, b) > tolerance / 256)
    return result;
  result.State                     = PointState::NearBoundary;
  result.HasBoundaryWitness        = true;
  result.WitnessFirst              = a;
  result.WitnessLast               = b;
  result.WitnessBracketWidth       = distance(a, b);
  result.WitnessDistanceUpperBound = upper;
  return result;
}

namespace
{
bool contains(const Loop& loop, const gp_Pnt2d& point)
{
  bool inside = false;
  for (size_t i = 0, j = loop.Polygon.size() - 1; i < loop.Polygon.size(); j = i++)
  {
    const auto& a = loop.Polygon[i];
    const auto& b = loop.Polygon[j];
    if ((a.Y() > point.Y()) != (b.Y() > point.Y())
        && point.X() < (b.X() - a.X()) * (point.Y() - a.Y()) / (b.Y() - a.Y()) + a.X())
      inside = !inside;
  }
  return inside;
}

void pcurve(const Arc& arc, const TopoDS_Face& cap, const gp_Pln& plane, double tolerance)
{
  const int                    count = int(arc.Points.size());
  NCollection_Array1<gp_Pnt2d> uv(1, count);
  NCollection_Array1<double>   knots(1, count);
  NCollection_Array1<int>      multiplicities(1, count);
  for (int i = 1; i <= count; ++i)
  {
    double u, v;
    ElSLib::Parameters(plane, arc.Points[size_t(i - 1)], u, v);
    uv(i).SetCoord(u, v);
    knots(i)          = i - 1;
    multiplicities(i) = i == 1 || i == count ? 2 : 1;
  }
  const occ::handle<Geom2d_BSplineCurve> curve =
    new Geom2d_BSplineCurve(uv, knots, multiplicities, 1);
  BRep_Builder builder;
  builder.UpdateEdge(arc.Edge, curve, cap, tolerance);
  builder.Range(arc.Edge, cap, 0, count - 1);
  builder.SameRange(arc.Edge, true);
  builder.SameParameter(arc.Edge, true);
}
} // namespace

BRepClass3d_ImplicitBox::Boundary BRepClass3d_ImplicitBox::Build(const TopoDS_Face&  fieldFace,
                                                                 const TopoDS_Shape& box,
                                                                 Operation           operation,
                                                                 double              tolerance,
                                                                 int                 maximumCells)
{
  if (!std::isfinite(tolerance) || tolerance <= 0)
    throw Standard_DomainError("Mixed Boolean tolerance must be finite and positive");
  if (operation != Operation::Common && operation != Operation::BoxMinusField)
    throw Standard_DomainError("Unknown mixed Boolean operation");
  const auto data = occ::down_cast<BRep_TImplicitFace>(fieldFace.TShape());
  require(!data.IsNull() && typeid(*data) == typeid(BRep_TImplicitFace)
            && typeid(*data->Field()) == typeid(Geom_ImplicitGyroid)
            && !TopoDS_Iterator(fieldFace).More() && fieldFace.Location().IsIdentity()
            && data->Location().IsIdentity() && fieldFace.Orientation() == TopAbs_FORWARD,
          "Mixed box Boolean requires a childless identity-placed forward native gyroid");
  require(!box.IsNull() && box.ShapeType() == TopAbs_SOLID && box.Orientation() == TopAbs_FORWARD
            && BRepCheck_Analyzer(box).IsValid(),
          "Operand must be one valid forward conventional box solid");
  Map vertices, edges, faces, shells;
  TopExp::MapShapes(box, TopAbs_VERTEX, vertices);
  TopExp::MapShapes(box, TopAbs_EDGE, edges);
  TopExp::MapShapes(box, TopAbs_FACE, faces);
  TopExp::MapShapes(box, TopAbs_SHELL, shells);
  require(vertices.Extent() == 8 && edges.Extent() == 12 && faces.Extent() == 6
            && shells.Extent() == 1,
          "Only six-face axis-aligned boxes are supported");
  gp_Pnt low = BRep_Tool::Pnt(TopoDS::Vertex(vertices(1))), high = low;
  double magnitude = 1;
  for (const auto& vertex : vertices)
  {
    const auto p = BRep_Tool::Pnt(TopoDS::Vertex(vertex));
    for (int j = 1; j <= 3; ++j)
    {
      require(std::isfinite(p.Coord(j)), "Non-finite box coordinate");
      low.SetCoord(j, std::min(low.Coord(j), p.Coord(j)));
      high.SetCoord(j, std::max(high.Coord(j), p.Coord(j)));
      magnitude = std::max(magnitude, std::abs(p.Coord(j)));
    }
  }
  const double reserve     = 512 * std::numeric_limits<double>::epsilon() * magnitude;
  const double recognition = std::min(1.e-7, tolerance / 64);
  require(tolerance > 4096 * reserve, "Mixed tolerance below coordinate precision policy");
  for (int j = 1; j <= 3; ++j)
    require(high.Coord(j) - low.Coord(j) > 4 * tolerance && low.Coord(j) >= data->Minimum().Coord(j)
              && high.Coord(j) <= data->Maximum().Coord(j),
            "Box must fit inside the field domain and exceed four tolerances per axis");
  std::array<bool, 8> cornerSeen{};
  for (const auto& vertex : vertices)
  {
    const auto p      = BRep_Tool::Pnt(TopoDS::Vertex(vertex));
    int        corner = 0;
    for (int j = 0; j < 3; ++j)
    {
      require(std::min(std::abs(p.Coord(j + 1) - low.Coord(j + 1)),
                       std::abs(p.Coord(j + 1) - high.Coord(j + 1)))
                <= recognition,
              "Box vertices must be axis-aligned corners");
      if (std::abs(p.Coord(j + 1) - high.Coord(j + 1)) <= recognition)
        corner |= 1 << j;
    }
    require(!cornerSeen[size_t(corner)], "Duplicate box corner");
    cornerSeen[size_t(corner)] = true;
  }
  for (const auto& edge : edges)
  {
    const auto              e = TopoDS::Edge(edge);
    const BRepAdaptor_Curve curve(e);
    require(curve.GetType() == GeomAbs_Line, "Only straight-edged boxes are supported");
    TopoDS_Vertex a, b;
    TopExp::Vertices(e, a, b);
    const auto p = BRep_Tool::Pnt(a), q = BRep_Tool::Pnt(b);
    const auto first = curve.Value(curve.FirstParameter()),
               last  = curve.Value(curve.LastParameter());
    require(std::min(std::max(first.Distance(p), last.Distance(q)),
                     std::max(first.Distance(q), last.Distance(p)))
              <= recognition,
            "Box edge geometry does not match its corner vertices");
    int varying = 0;
    for (int j = 1; j <= 3; ++j)
      varying += std::abs(p.Coord(j) - q.Coord(j)) > recognition;
    require(varying == 1, "Box edges must connect adjacent axis-aligned corners");
  }
  std::vector<Side>   sides;
  std::array<bool, 6> seen{};
  int                 source = 0;
  for (const auto& shape : faces)
  {
    const auto face = TopoDS::Face(shape);
    require(face.Orientation() == TopAbs_FORWARD || face.Orientation() == TopAbs_REVERSED,
            "Unsupported box face orientation");
    BRepAdaptor_Surface adaptor(face);
    require(adaptor.GetType() == GeomAbs_Plane, "Only planar box faces are supported");
    const auto plane = adaptor.Plane();
    gp_Vec     normal(plane.Axis().Direction());
    if (face.Orientation() == TopAbs_REVERSED)
      normal.Reverse();
    int axis = -1;
    for (int j = 0; j < 3; ++j)
      if (std::abs(normal.Coord(j + 1)) > 1 - 1.e-12)
        axis = j;
    require(axis >= 0, "Only axis-aligned box planes are supported");
    const int sign  = normal.Coord(axis + 1) > 0 ? 1 : -1;
    const int index = 2 * axis + int(sign > 0);
    require(!seen[size_t(index)], "Duplicate box side");
    seen[size_t(index)]     = true;
    const double coordinate = sign > 0 ? high.Coord(axis + 1) : low.Coord(axis + 1);
    require(std::abs(plane.Location().Coord(axis + 1) - coordinate) <= recognition,
            "Box plane does not match an outward boundary");
    for (TopExp_Explorer it(face, TopAbs_VERTEX); it.More(); it.Next())
    {
      const auto point = BRep_Tool::Pnt(TopoDS::Vertex(it.Current()));
      require(std::abs(BRep_Tool::Pnt(TopoDS::Vertex(it.Current())).Coord(axis + 1) - coordinate)
                <= recognition,
              "Box face vertices do not lie on the same boundary side");
      require(plane.Distance(point) <= recognition, "Box plane does not match its corner geometry");
    }
    Map sideEdges;
    TopExp::MapShapes(face, TopAbs_EDGE, sideEdges);
    require(sideEdges.Extent() == 4, "Box faces must be four-edged rectangles");
    int count = 0;
    for (TopoDS_Iterator it(face); it.More(); it.Next())
      ++count;
    require(count == 1, "Box faces with holes are unsupported");
    gp_Pnt origin(0, 0, 0);
    origin.SetCoord(axis + 1, coordinate);
    Side side;
    side.Axis           = axis;
    side.Sign           = sign;
    side.Source         = source++;
    side.Settings.Frame = gp_Ax2(origin, gp_Dir(unit(axis, sign)), gp_Dir(unit((axis + 1) % 3)));
    side.Plane          = gp_Pln(side.Settings.Frame);
    const int u = (axis + 1) % 3, v = (axis + 2) % 3;
    side.Settings.UMinimum = low.Coord(u + 1);
    side.Settings.UMaximum = high.Coord(u + 1);
    side.Settings.VMinimum = sign > 0 ? low.Coord(v + 1) : -high.Coord(v + 1);
    side.Settings.VMaximum = sign > 0 ? high.Coord(v + 1) : -low.Coord(v + 1);
    // Reserve tolerance for replacing independent face-endpoint brackets by
    // the canonical box-edge root, never distance-welding unrelated vertices.
    side.Settings.DistanceTolerance = tolerance * .75;
    side.Settings.MaximumCells      = maximumCells;
    sides.push_back(side);
  }

  BRep_Builder builder;
  Boundary     result;
  builder.MakeCompound(result.Shape);
  std::vector<Node> nodes;
  const double      precision = std::max(reserve, std::min(1.e-7, tolerance / 4096));
  const auto&       field     = *static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  const double      k = 2 * std::acos(-1.) / field.Period(), hessian = 4 * k * k;
  require(std::isfinite(hessian) && hessian > 0, "Unrepresentable gyroid derivative bounds");
  const double fieldReserve =
    512 * std::numeric_limits<double>::epsilon() * (1 + std::abs(field.IsoValue()))
    + field.LipschitzBound() * reserve;
  for (int c = 0; c < 8; ++c)
  {
    const gp_Pnt point(c & 1 ? high.X() : low.X(),
                       c & 2 ? high.Y() : low.Y(),
                       c & 4 ? high.Z() : low.Z());
    require(std::abs(field.Value(point)) > fieldReserve, "Box corner on field is unresolved");
    addNode(nodes, point, precision);
  }
  std::vector<Arc>     arcs;
  std::vector<BoxEdge> boxEdges;
  const bool           negative = operation == Operation::Common;
  for (int axis = 0; axis < 3; ++axis)
    for (int corner = 0; corner < 8; ++corner)
    {
      if (corner & (1 << axis))
        continue;
      BoxEdge             edge{corner, corner | (1 << axis), {}, {}};
      std::vector<gp_Pnt> points;
      int                 budget = 65536;
      roots(field,
            nodes[size_t(edge.First)].Point,
            nodes[size_t(edge.Last)].Point,
            hessian,
            tolerance / 256,
            fieldReserve,
            points,
            budget);
      std::vector<int> chain{edge.First};
      for (const auto& point : points)
      {
        const int id = addNode(nodes, point, precision, int(boxEdges.size()));
        edge.Roots.push_back(id);
        chain.push_back(id);
        ++result.BoxEdgeRootCount;
      }
      chain.push_back(edge.Last);
      for (size_t i = 0; i + 1 < chain.size(); ++i)
      {
        const auto& a = nodes[size_t(chain[i])].Point;
        const auto& b = nodes[size_t(chain[i + 1])].Point;
        require(a.Distance(b) > tolerance / 128,
                "Box roots too close to resolve retained intervals");
        const double value = field.Value(gp_Pnt((a.XYZ() + b.XYZ()) * .5));
        require(std::abs(value) > fieldReserve, "Box interval classification unresolved");
        if ((value < 0) == negative)
          edge.Arcs.push_back(addArc(arcs, nodes, chain[i], chain[i + 1], {a, b}, tolerance));
      }
      boxEdges.push_back(edge);
    }
  std::vector<int> rootUses(nodes.size(), 0);
  std::vector<Use> implicitUses;
  for (auto& side : sides)
  {
    auto section = BRep_ImplicitPlaneSection::Trace(fieldFace, side.Settings);
    if (!section.Resolved)
    {
      std::ostringstream message;
      message << "Box face " << side.Source << " unresolved: cells=" << section.UnresolvedCells
              << ", critical=" << section.PotentialCriticalCells
              << ", ambiguous=" << section.AmbiguousCells;
      throw Standard_NotImplemented(message.str().c_str());
    }
    const gp_Vec normal(side.Plane.Axis().Direction());
    gp_Pnt       center((low.XYZ() + high.XYZ()) * .5);
    center.SetCoord(side.Axis + 1, side.Plane.Location().Coord(side.Axis + 1));
    for (const auto& edge : boxEdges)
    {
      const auto& a = nodes[size_t(edge.First)].Point;
      const auto& b = nodes[size_t(edge.Last)].Point;
      if (std::abs(a.Coord(side.Axis + 1) - center.Coord(side.Axis + 1)) > reserve
          || std::abs(b.Coord(side.Axis + 1) - center.Coord(side.Axis + 1)) > reserve)
        continue;
      const gp_Vec tangent(a, b), inward(gp_Pnt((a.XYZ() + b.XYZ()) * .5), center);
      const bool   forward = normal.Crossed(tangent).Dot(inward) > 0;
      for (int arc : edge.Arcs)
        side.Uses.push_back({arc, forward});
    }
    double endpointShift = 0;
    for (auto& curve : section.Curves)
    {
      int first, last;
      if (curve.Closed)
        first = last = addNode(nodes, curve.Points.front(), precision);
      else
      {
        auto endpoint = [&](gp_Pnt& point) {
          int found = -1;
          for (const auto& edge : boxEdges)
          {
            const auto& a = nodes[size_t(edge.First)].Point;
            const auto& b = nodes[size_t(edge.Last)].Point;
            if (std::abs(a.Coord(side.Axis + 1) - center.Coord(side.Axis + 1)) > reserve
                || std::abs(b.Coord(side.Axis + 1) - center.Coord(side.Axis + 1)) > reserve)
              continue;
            for (int root : edge.Roots)
              if (point.Distance(nodes[size_t(root)].Point) <= tolerance / 8)
              {
                require(found < 0, "Ambiguous correspondence between face endpoint and box roots");
                found = root;
              }
          }
          require(found >= 0, "Face endpoint has no certified box-edge root");
          const double shift = point.Distance(nodes[size_t(found)].Point);
          endpointShift      = std::max(endpointShift, shift);
          point              = nodes[size_t(found)].Point;
          ++rootUses[size_t(found)];
          return found;
        };
        first = endpoint(curve.Points.front());
        last  = endpoint(curve.Points.back());
        require(first != last, "Open contour returned to the same box-edge root");
      }
      const int    index = addArc(arcs, nodes, first, last, curve.Points, tolerance);
      const gp_Vec tangent(curve.Points[0], curve.Points[1]);
      const auto   desired = normal.Crossed(field.Gradient(curve.Points[0])) * (negative ? 1 : -1);
      const double alignment = tangent.Dot(desired);
      require(std::isfinite(alignment) && std::abs(alignment) > fieldReserve * precision,
              "Cannot orient gyroid seam");
      const bool forward = alignment > 0;
      side.Uses.push_back({index, forward});
      implicitUses.push_back({index, !forward});
      ++result.SeamCount;
      for (const auto& point : curve.Points)
        section.MaximumResidual = std::max(section.MaximumResidual, std::abs(field.Value(point)));
    }
    section.GeometricBound += endpointShift + reserve;
    require(section.GeometricBound <= tolerance,
            "Canonical endpoint adjustment exhausted seam tolerance");
    result.GeometricBound = std::max(result.GeometricBound, section.GeometricBound);
    result.Sections.push_back(section);
    auto             loops = assemble(side.Uses, arcs, &side.Plane);
    std::vector<int> parent(loops.size(), -1);
    for (size_t i = 0; i < loops.size(); ++i)
    {
      require(std::abs(loops[i].Area) > precision * precision, "Degenerate planar trim region");
      if (loops[i].Area > 0)
        continue;
      double smallest = std::numeric_limits<double>::infinity();
      for (size_t j = 0; j < loops.size(); ++j)
        if (loops[j].Area > 0 && loops[j].Area < smallest
            && contains(loops[j], loops[i].Polygon.front()))
        {
          parent[i] = int(j);
          smallest  = loops[j].Area;
        }
      require(parent[i] >= 0, "Planar trim hole has no enclosing retained region");
    }
    for (size_t i = 0; i < loops.size(); ++i)
    {
      if (loops[i].Area < 0)
        continue;
      TopoDS_Face cap;
      builder.MakeFace(cap, new Geom_Plane(side.Plane), precision);
      auto addLoop = [&](const Loop& loop) {
        for (const auto& use : loop.Uses)
          pcurve(arcs[size_t(use.Index)], cap, side.Plane, tolerance);
        builder.Add(cap, loop.Wire);
      };
      addLoop(loops[i]);
      for (size_t j = 0; j < loops.size(); ++j)
        if (parent[j] == int(i))
          addLoop(loops[j]);
      require(BRepCheck_Analyzer(cap).IsValid(),
              "Constructed planar cap failed conventional validation");
      result.Faces.push_back(cap);
      result.Sources.push_back(side.Source);
      builder.Add(result.Shape, cap);
    }
  }
  for (size_t i = 0; i < rootUses.size(); ++i)
    if (nodes[i].RootEdge >= 0)
      require(rootUses[i] == 2, "Box-edge root must join exactly two adjacent face contours");
  TopoDS_Face implicit;
  builder.MakeImplicitFace(implicit, data->Field(), low, high, TopLoc_Location(), tolerance);
  if (!negative)
    implicit.Reverse();
  for (const auto& loop : assemble(implicitUses, arcs, nullptr))
    builder.Add(implicit, loop.Wire);
  result.Faces.push_back(implicit);
  result.Sources.push_back(-1);
  builder.Add(result.Shape, implicit);
  Map boundaryEdges;
  TopExp::MapShapes(result.Shape, TopAbs_EDGE, boundaryEdges);
  std::vector<int> uses(size_t(boundaryEdges.Extent()), 0), balance(uses.size(), 0);
  for (TopExp_Explorer it(result.Shape, TopAbs_EDGE); it.More(); it.Next())
  {
    const size_t index = size_t(boundaryEdges.FindIndex(it.Current()) - 1);
    ++uses[index];
    require(it.Current().Orientation() == TopAbs_FORWARD
              || it.Current().Orientation() == TopAbs_REVERSED,
            "Unsupported boundary edge orientation");
    balance[index] += it.Current().Orientation() == TopAbs_FORWARD ? 1 : -1;
  }
  for (size_t i = 0; i < uses.size(); ++i)
    require(uses[i] == 2 && balance[i] == 0, "Unpaired or inconsistently oriented hybrid edge");
  result.SeamBoundVerified = true;
  return result;
}
