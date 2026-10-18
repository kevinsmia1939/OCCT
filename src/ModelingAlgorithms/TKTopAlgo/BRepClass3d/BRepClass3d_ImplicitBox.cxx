// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitBox.hxx>
#include <BRepClass3d_ImplicitBoxIO.hxx>
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
#include <set>
#include <sstream>
#include <typeinfo>

namespace
{
using Map       = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
using Boundary  = BRepClass3d_ImplicitBox::Boundary;
using Operation = BRepClass3d_ImplicitBox::Operation;
using CornerPolicy = BRepClass3d_ImplicitBox::CornerPolicy;

bool exactIdentity(const TopLoc_Location& location)
{
  const auto matrix = location.Transformation();
  for (int row = 1; row <= 3; ++row)
    for (int column = 1; column <= 4; ++column)
      if (matrix.Value(row, column) != (row == column ? 1. : 0.))
        return false;
  return true;
}

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

// Endpoint zeros are anchored identities, not interior roots to duplicate on
// each incident edge. A recognized midpoint is inserted once in edge order;
// its two child intervals prove that no additional root was hidden beside it.
void cornerRoots(const Geom_ImplicitGyroid& field, const gp_Pnt& first, const gp_Pnt& last,
                 double hessian, double bracket, double reserve, double gradientReserve, std::vector<gp_Pnt>& out,
                 int& budget, int depth = 0)
{
  require(--budget >= 0 && depth <= 32, "Corner-aware box-edge root budget exhausted");
  const gp_Pnt mid((first.XYZ() + last.XYZ()) * .5);
  const double length = first.Distance(last);
  require(std::isfinite(length) && length > 0, "Corner-aware edge coordinate collapse");
  double va = 0, vb = 0, vm = 0;
  gp_Vec unused;
  const bool za = field.LatticeJet(first, va, unused) && va == 0;
  const bool zb = field.LatticeJet(last, vb, unused) && vb == 0;
  const bool zm = field.LatticeJet(mid, vm, unused) && vm == 0;
  const double fa = za ? 0 : field.Value(first), fb = zb ? 0 : field.Value(last);
  const double fm = zm ? 0 : field.Value(mid);
  require((za || std::abs(fa) > reserve) && (zb || std::abs(fb) > reserve),
          "Unrecognized near-zero corner/edge subdivision node");
  if (std::abs(fm) > field.LipschitzBound() * length * .5 + reserve)
    return;
  const double derivative = field.Gradient(mid).Dot(gp_Vec(first, last) / length);
  const double slope = std::abs(derivative) - hessian * length * .5 - gradientReserve;
  if (slope > 0)
  {
    require(!(za && zb), "Two endpoint zeros contradict edge monotonicity");
    if (za || zb)
    {
      require(std::signbit(za ? fb : fa) == std::signbit(za ? derivative : -derivative),
              "Endpoint zero contradicts the certified inward edge sign");
      return;
    }
    if (std::signbit(fa) == std::signbit(fb))
      return;
    if (zm)
      out.push_back(mid);
    else
    {
      gp_Pnt a = first, b = last;
      double va = fa;
      for (int iteration = 0; iteration < 100 && a.Distance(b) > bracket; ++iteration)
      {
        const gp_Pnt p((a.XYZ() + b.XYZ()) * .5);
        double exact = 0;
        if (field.LatticeJet(p, exact, unused) && exact == 0)
        {
          a = b = p;
          break;
        }
        const double value = field.Value(p);
        if (std::abs(value) <= reserve)
        {
          require(2 * (std::abs(value) + reserve) / slope <= bracket,
                  "Corner-edge root location uncertainty exceeds its bracket budget");
          a = b = p; // A bounded location witness, NOT a newly recognized exact zero.
          break;
        }
        if (std::signbit(va) != std::signbit(value))
          b = p;
        else
        {
          a = p;
          va = value;
        }
      }
      require(a.Distance(b) <= bracket, "Corner-edge root failed its distance bracket");
      out.emplace_back((a.XYZ() + b.XYZ()) * .5);
    }
    return;
  }
  require(length > bracket && mid.Distance(first) > 0 && mid.Distance(last) > 0,
          "Tangent/critical corner-edge interval unresolved");
  cornerRoots(field, first, mid, hessian, bracket, reserve, gradientReserve, out, budget, depth + 1);
  if (zm)
    out.push_back(mid);
  cornerRoots(field, mid, last, hessian, bracket, reserve, gradientReserve, out, budget, depth + 1);
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
  require(boundary.Policy == CornerPolicy::Reject || boundary.Policy == CornerPolicy::Regularized,
          "Unknown mathematical-region corner policy");
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
  CornerModel corners;
  double cornerReserve = coordinateReserve;
  if (boundary.Policy == CornerPolicy::Regularized)
  {
    const double span = std::min({high.X() - low.X(), high.Y() - low.Y(), high.Z() - low.Z()});
    corners = CheckCornerGerms(face, span / 8);
    require(corners.Resolved, "Regularized classification needs checked local corner identities");
    cornerReserve = std::max(coordinateReserve, corners.CoordinateReserve);
    for (size_t index = 0; index < 8; ++index)
    {
      const auto& corner = corners.Corners[index];
      if (!corner.ExactZero || corner.Kind != (negative ? CornerKind::PositiveCone
                                                       : CornerKind::NegativeCone))
        continue;
      double exact = 0;
      gp_Vec jet;
      require(field.LatticeJet(corner.Anchor, exact, jet) && exact == 0,
              "Excluded corner lost its exact field identity");
      double inward = std::numeric_limits<double>::infinity();
      for (int axis = 0; axis < 3; ++axis)
        inward = std::min(inward, jet.Coord(axis + 1) * (index & (1 << axis) ? -1 : 1)
                                  * (negative ? 1 : -1));
      const double error = hessian * 8 * cornerReserve
                           + 1024 * std::numeric_limits<double>::epsilon() * lipschitz;
      // For inward d, G(a+d) >= m*|d| - H*|d|^2/2 since G(a)=0 exactly.
      // Thus the entire clipped ball of radius m/H has G>0 except at a.
      // Its anchor is NOT in closure(interior(material)). This is a region
      // exclusion bound, not residual-as-distance or a nearest-point estimate.
      const double radius = std::max(0., std::min(span / 8, (inward - 4 * error) / hessian)
                                           - 32 * cornerReserve);
      const double exclusion = std::max(0., radius - distance(point, corner.Anchor)
                                            - 16 * cornerReserve);
      result.DistanceLowerBound = std::max(result.DistanceLowerBound, exclusion);
      if (exclusion > tolerance + coordinateReserve)
      {
        result.State = PointState::Outside;
        return result;
      }
    }
  }
  if (result.DistanceLowerBound > tolerance + coordinateReserve)
  {
    result.State = insideBox && retained && signResolved ? PointState::Inside : PointState::Outside;
    return result;
  }
  if (boundary.Policy == CornerPolicy::Regularized)
    for (const auto& corner : corners.Corners)
    {
      // Checked split/retained cones have strict material points arbitrarily
      // close inside the box. Their exact corner is therefore a true boundary
      // witness. Never use an excluded closed-inequality point as a witness.
      if (!corner.ExactZero || corner.Kind == (negative ? CornerKind::PositiveCone
                                                       : CornerKind::NegativeCone))
        continue;
      const double upper = distance(point, corner.Anchor) + 16 * cornerReserve;
      if (upper > tolerance)
        continue;
      result.State = PointState::NearBoundary;
      result.HasBoundaryWitness = true;
      result.WitnessFirst = result.WitnessLast = corner.Anchor;
      result.WitnessDistanceUpperBound = upper;
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

Side canonicalSide(int           axis,
                   int           sign,
                   int           source,
                   const gp_Pnt& low,
                   const gp_Pnt& high,
                   double        tolerance,
                   int           maximumCells)
{
  gp_Pnt origin(0, 0, 0);
  origin.SetCoord(axis + 1, sign > 0 ? high.Coord(axis + 1) : low.Coord(axis + 1));
  Side side;
  side.Axis           = axis;
  side.Sign           = sign;
  side.Source         = source;
  side.Settings.Frame = gp_Ax2(origin, gp_Dir(unit(axis, sign)), gp_Dir(unit((axis + 1) % 3)));
  side.Plane          = gp_Pln(side.Settings.Frame);
  const int u = (axis + 1) % 3, v = (axis + 2) % 3;
  side.Settings.UMinimum          = low.Coord(u + 1);
  side.Settings.UMaximum          = high.Coord(u + 1);
  side.Settings.VMinimum          = sign > 0 ? low.Coord(v + 1) : -high.Coord(v + 1);
  side.Settings.VMaximum          = sign > 0 ? high.Coord(v + 1) : -low.Coord(v + 1);
  side.Settings.DistanceTolerance = tolerance * .75;
  side.Settings.MaximumCells      = maximumCells;
  return side;
}

Boundary construct(const TopoDS_Face&,
                   const gp_Pnt&,
                   const gp_Pnt&,
                   std::vector<Side>,
                   Operation,
                   double,
                   double,
                   CornerPolicy = CornerPolicy::Reject);
} // namespace

BRepClass3d_ImplicitBox::Boundary BRepClass3d_ImplicitBox::Build(const TopoDS_Face&  fieldFace,
                                                                 const TopoDS_Shape& box,
                                                                 Operation           operation,
                                                                 double              tolerance,
                                                                 int                 maximumCells)
{
  return Build(fieldFace, box, operation, tolerance, maximumCells, CornerPolicy::Reject);
}

BRepClass3d_ImplicitBox::Boundary BRepClass3d_ImplicitBox::Build(const TopoDS_Face& fieldFace,
  const TopoDS_Shape& box, Operation operation, double tolerance, int maximumCells, CornerPolicy policy)
{
  if (!std::isfinite(tolerance) || tolerance <= 0)
    throw Standard_DomainError("Mixed Boolean tolerance must be finite and positive");
  if (operation != Operation::Common && operation != Operation::BoxMinusField)
    throw Standard_DomainError("Unknown mixed Boolean operation");
  if (policy != CornerPolicy::Reject && policy != CornerPolicy::Regularized)
    throw Standard_DomainError("Unknown mixed corner policy");
  const auto data = occ::down_cast<BRep_TImplicitFace>(fieldFace.TShape());
  require(!data.IsNull() && typeid(*data) == typeid(BRep_TImplicitFace)
            && typeid(*data->Field()) == typeid(Geom_ImplicitGyroid)
            && data->ExcludedDomains().empty()
            && !TopoDS_Iterator(fieldFace).More()
            && (policy == CornerPolicy::Regularized
                  ? exactIdentity(fieldFace.Location() * data->Location())
                  : fieldFace.Location().IsIdentity() && data->Location().IsIdentity())
            && fieldFace.Orientation() == TopAbs_FORWARD,
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
  if (policy == CornerPolicy::Regularized)
  {
    require(exactIdentity(box.Location()), "Corner-aware input box must have identity placement");
    for (const auto& vertex : vertices)
    {
      require(exactIdentity(vertex.Location()), "Corner-aware box vertices must be local identities");
      const auto p = BRep_Tool::Pnt(TopoDS::Vertex(vertex));
      for (int axis = 1; axis <= 3; ++axis)
        require(p.Coord(axis) == low.Coord(axis) || p.Coord(axis) == high.Coord(axis),
                "Corner-aware box recognition must not move a near-corner vertex");
    }
  }
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
    if (policy == CornerPolicy::Regularized)
    {
      require(exactIdentity(e.Location()), "Corner-aware box edges need identity placement");
      const auto line = curve.Line();
      for (int axis = 1; axis <= 3; ++axis)
      {
        const bool varying = p.Coord(axis) != q.Coord(axis);
        require(varying ? std::abs(line.Direction().Coord(axis)) == 1
                        : line.Direction().Coord(axis) == 0 && line.Location().Coord(axis) == p.Coord(axis),
                "Corner-aware box lines must agree exactly with coordinate-plane intersections");
      }
    }
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
    if (policy == CornerPolicy::Regularized)
    {
      require(exactIdentity(face.Location()) && plane.Location().Coord(axis + 1) == coordinate,
              "Corner-aware box planes need exact intercepts and identity placement");
      for (int j = 0; j < 3; ++j)
        require(normal.Coord(j + 1) == (j == axis ? sign : 0),
                "Corner-aware box planes must be exact coordinate planes");
    }
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
    sides.push_back(canonicalSide(axis, sign, source++, low, high, tolerance, maximumCells));
  }

  return construct(fieldFace, low, high, std::move(sides), operation, tolerance, reserve, policy);
}

namespace
{
Boundary construct(const TopoDS_Face& fieldFace,
                   const gp_Pnt&      low,
                   const gp_Pnt&      high,
                   std::vector<Side>  sides,
                   Operation          operation,
                   double             tolerance,
                   double             reserve,
                   CornerPolicy       policy)
{
  const auto   data = occ::down_cast<BRep_TImplicitFace>(fieldFace.TShape());
  BRep_Builder builder;
  Boundary     result;
  result.Policy = policy;
  builder.MakeCompound(result.Shape);
  std::vector<Node> nodes;
  const double      precision = std::max(reserve, std::min(1.e-7, tolerance / 4096));
  const auto&       field     = *static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  const double      k = 2 * std::acos(-1.) / field.Period(), hessian = 4 * k * k;
  require(std::isfinite(hessian) && hessian > 0, "Unrepresentable gyroid derivative bounds");
  const double fieldReserve =
    512 * std::numeric_limits<double>::epsilon() * (1 + std::abs(field.IsoValue()))
    + field.LipschitzBound() * reserve;
  BRepClass3d_ImplicitBox::CornerModel cornerModel;
  if (policy == CornerPolicy::Regularized)
  {
    TopoDS_Face boxed;
    builder.MakeImplicitFace(boxed, data->Field(), low, high, TopLoc_Location(), tolerance);
    cornerModel = BRepClass3d_ImplicitBox::CheckCornerGerms(boxed, tolerance);
    require(cornerModel.Resolved, "Corner-aware box neighborhoods unresolved");
    result.ExactCornerCount = cornerModel.ExactZeroCount;
    for (const auto& corner : cornerModel.Corners)
      if (corner.ExactZero
          && (operation == Operation::Common
                ? corner.Kind == BRepClass3d_ImplicitBox::CornerKind::PositiveCone
                : corner.Kind == BRepClass3d_ImplicitBox::CornerKind::NegativeCone))
        result.ExcludedCornerContacts.push_back(corner.Anchor);
  }
  for (int c = 0; c < 8; ++c)
  {
    const gp_Pnt point(c & 1 ? high.X() : low.X(),
                       c & 2 ? high.Y() : low.Y(),
                       c & 4 ? high.Z() : low.Z());
    require(std::abs(field.Value(point)) > fieldReserve
              || (policy == CornerPolicy::Regularized && cornerModel.Corners[size_t(c)].ExactZero),
            "Box corner on field is unresolved");
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
      if (policy == CornerPolicy::Regularized)
        cornerRoots(field, nodes[size_t(edge.First)].Point, nodes[size_t(edge.Last)].Point,
                    hessian, tolerance / 256, fieldReserve,
                    hessian * 8 * reserve + 1024 * std::numeric_limits<double>::epsilon() * field.LipschitzBound(),
                    points, budget);
      else
        roots(field,
            nodes[size_t(edge.First)].Point,
            nodes[size_t(edge.Last)].Point,
            hessian,
            tolerance / 256,
            fieldReserve,
            points,
            budget,
            0);
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
    side.Settings.AllowExactLatticeNodes = policy == CornerPolicy::Regularized;
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
    auto cornerAt = [&](const gp_Pnt& p, bool entering) {
      int found = -1;
      for (int index = 0; index < 8; ++index)
      {
        const auto& corner = cornerModel.Corners[size_t(index)];
        if (!corner.ExactZero || p.X() != corner.Anchor.X() || p.Y() != corner.Anchor.Y()
            || p.Z() != corner.Anchor.Z())
          continue;
        const int a = (side.Axis + 1) % 3, b = (side.Axis + 2) % 3;
        const bool positiveA = corner.InwardDerivativeLower[size_t(a)] > 0;
        const bool positiveB = corner.InwardDerivativeLower[size_t(b)] > 0;
        require((positiveA != positiveB) == entering,
                "Corner face-contour incidence disagrees with checked derivative cone");
        require(found < 0, "Duplicate exact corner identity");
        found = index;
      }
      return found;
    };
    for (const auto& contact : section.PointContacts)
    {
      require(policy == CornerPolicy::Regularized && cornerAt(contact, false) >= 0,
              "Plane contact lacks a checked canonical corner");
      ++result.PlanarCornerContacts;
    }
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
          if (policy == CornerPolicy::Regularized)
            found = cornerAt(point, true);
          if (found >= 0)
          {
            ++rootUses[size_t(found)];
            ++result.CornerContourEnds;
            return found; // Reuse the canonical corner vertex, never an edge-specific copy.
          }
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
  if (policy == CornerPolicy::Regularized)
    for (size_t i = 0; i < 8; ++i)
      if (cornerModel.Corners[i].ExactZero)
        require(rootUses[i] == (cornerModel.Corners[i].Kind == BRepClass3d_ImplicitBox::CornerKind::SplitCone ? 2 : 0),
                "Canonical corner must match its checked cross-face contour incidence");
  TopoDS_Face implicit;
  builder.MakeImplicitFace(implicit, data->Field(), low, high, TopLoc_Location(), tolerance);
  if (policy == CornerPolicy::Regularized)
    implicit.TShape(new BRep_TImplicitFace(data->Field(), low, high, tolerance,
      BRepClass3d_ImplicitBox::SurfaceCornerExclusions(implicit)));
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

struct CurveRecord
{
  std::vector<gp_Pnt>      Points;
  std::vector<long double> Knots;
};

CurveRecord curveRecord(std::vector<gp_Pnt> points)
{
  CurveRecord result;
  for (const auto& p : points)
    if (result.Points.empty() || p.Distance(result.Points.back()) > 0)
      result.Points.push_back(p);
  require(result.Points.size() >= 2, "Collapsed correspondence curve");
  result.Knots.push_back(0);
  for (size_t i = 1; i < result.Points.size(); ++i)
  {
    const auto&       a      = result.Points[i - 1];
    const auto&       b      = result.Points[i];
    const long double length = std::hypot(static_cast<long double>(b.X()) - a.X(),
                                          static_cast<long double>(b.Y()) - a.Y(),
                                          static_cast<long double>(b.Z()) - a.Z());
    const long double next   = result.Knots.back() + length;
    require(std::isfinite(next) && next > result.Knots.back(),
            "Correspondence parameter precision exhausted");
    result.Knots.push_back(next);
  }
  const auto total = result.Knots.back();
  for (auto& t : result.Knots)
    t /= total;
  result.Knots.back() = 1;
  for (size_t i = 1; i < result.Knots.size(); ++i)
    require(result.Knots[i] > result.Knots[i - 1], "Correspondence parameter collapsed");
  return result;
}

gp_Pnt at(const CurveRecord& curve, long double t, size_t& segment)
{
  while (segment + 2 < curve.Knots.size() && t > curve.Knots[segment + 1])
    ++segment;
  const long double u =
    std::clamp((t - curve.Knots[segment]) / (curve.Knots[segment + 1] - curve.Knots[segment]),
               0.L,
               1.L);
  gp_Pnt p;
  for (int axis = 1; axis <= 3; ++axis)
    p.SetCoord(axis,
               double((1 - u) * curve.Points[segment].Coord(axis)
                      + u * curve.Points[segment + 1].Coord(axis)));
  return p;
}

double curveBound(const CurveRecord&                       a,
                  const CurveRecord&                       b,
                  double                                   reserve,
                  BRepClass3d_ImplicitBox::Correspondence& report)
{
  if (report.SampleTests + a.Points.size() + b.Points.size() > 8000000)
    throw Standard_NotImplemented("Correspondence sample budget exhausted");
  size_t i = 0, j = 0, sa = 0, sb = 0;
  double bound = 0;
  while (i < a.Knots.size() || j < b.Knots.size())
  {
    const long double t = j == b.Knots.size()        ? a.Knots[i++]
                          : i == a.Knots.size()      ? b.Knots[j++]
                          : a.Knots[i] <= b.Knots[j] ? a.Knots[i++]
                                                     : b.Knots[j++];
    ++report.SampleTests;
    bound = std::max(bound, at(a, t, sa).Distance(at(b, t, sb)) + reserve);
  }
  return bound;
}

std::vector<CurveRecord> curves(const Map& edges)
{
  std::vector<CurveRecord> result;
  for (const auto& shape : edges)
  {
    TopLoc_Location location;
    double          first, last;
    const auto      curve = occ::down_cast<Geom_BSplineCurve>(
      BRep_Tool::Curve(TopoDS::Edge(shape), location, first, last));
    require(!curve.IsNull() && location.IsIdentity(), "Unsupported correspondence curve");
    std::vector<gp_Pnt> points;
    for (int pole = 1; pole <= curve->NbPoles(); ++pole)
      points.push_back(curve->Pole(pole));
    result.push_back(curveRecord(std::move(points)));
  }
  return result;
}

using WireKey = std::vector<int>;
using FaceKey = std::pair<int, std::vector<WireKey>>;

std::vector<WireKey> wires(const TopoDS_Face&      face,
                           const Map&              edges,
                           const std::vector<int>& mapping,
                           const std::vector<int>& directions,
                           bool                    sortResult = true)
{
  std::vector<WireKey> result;
  for (TopoDS_Iterator wire(face); wire.More(); wire.Next())
  {
    WireKey key;
    for (TopoDS_Iterator edge(wire.Value()); edge.More(); edge.Next())
    {
      const int id = edges.FindIndex(edge.Value()) - 1;
      require(id >= 0 && mapping[size_t(id)] >= 0, "Incomplete correspondence edge map");
      key.push_back((mapping[size_t(id)] + 1) * directions[size_t(id)]
                    * (edge.Value().Orientation() == TopAbs_FORWARD ? 1 : -1));
    }
    require(!key.empty(), "Empty correspondence wire");
    std::rotate(key.begin(), std::min_element(key.begin(), key.end()), key.end());
    result.push_back(std::move(key));
  }
  if (sortResult)
    std::sort(result.begin(), result.end());
  return result;
}

std::vector<FaceKey> caps(const Boundary&         boundary,
                          const Map&              edges,
                          const std::vector<int>& mapping,
                          const std::vector<int>& directions)
{
  std::vector<FaceKey> result;
  for (size_t i = 0; i + 1 < boundary.Faces.size(); ++i)
  {
    const auto normal = BRepAdaptor_Surface(boundary.Faces[i]).Plane().Axis().Direction();
    int        axis   = 0;
    for (int j = 1; j <= 3; ++j)
      if (std::abs(normal.Coord(j)) > std::abs(normal.Coord(axis + 1)))
        axis = j - 1;
    const int side = 2 * axis + int(normal.Coord(axis + 1) > 0);
    result.push_back({side, wires(boundary.Faces[i], edges, mapping, directions)});
  }
  std::sort(result.begin(), result.end());
  return result;
}
} // namespace

BRepClass3d_ImplicitBox::Correspondence BRepClass3d_ImplicitBox::CheckCorrespondence(
  const Boundary& source,
  int             maximumCells,
  size_t          maximumComparisons)
{
  if (maximumCells < 4 || maximumCells > 512 || (maximumCells & (maximumCells - 1))
      || !maximumComparisons || maximumComparisons > 10000000)
    throw Standard_DomainError("Invalid trim correspondence settings");
  // Validate/private-copy the supported graph, without touching source caches.
  std::ostringstream record;
  BRepClass3d_ImplicitBoxIO::Write(source, record);
  std::istringstream input(record.str());
  Boundary           copy;
  BRepClass3d_ImplicitBoxIO::Read(input, copy);
  Correspondence report;
  const auto&    original   = copy.Faces.back();
  if (copy.Policy == CornerPolicy::Regularized)
  {
    const auto witness = CheckTopology(original, 32, 16, 4000000, copy.Policy);
    report.RegularityResolved = witness.HomotopyResolved && witness.ReferenceManifold;
    if (!report.RegularityResolved)
    {
      report.Reason = "fixed_corner_homotopy_unresolved: " + witness.Reason;
      return report;
    }
  }
  else
    report.RegularityResolved = CheckRegularity(original).Resolved;
  if (!report.RegularityResolved)
  {
    report.Reason = "field_box_regularity_unresolved";
    return report;
  }
  const auto   data = occ::down_cast<BRep_TImplicitFace>(original.TShape());
  const auto   low = data->Minimum(), high = data->Maximum();
  const double tolerance = data->Tolerance();
  double       magnitude = 1;
  for (int axis = 1; axis <= 3; ++axis)
    magnitude = std::max({magnitude, std::abs(low.Coord(axis)), std::abs(high.Coord(axis))});
  const double reserve     = 8192 * std::numeric_limits<double>::epsilon() * magnitude;
  report.MatchingTolerance = std::max(reserve * 16, tolerance / 1048576);
  try
  {
    TopoDS_Face bare;
    BRep_Builder().MakeImplicitFace(bare, data->Field(), low, high, TopLoc_Location(), tolerance);
    std::vector<Side> sides;
    for (int axis = 0; axis < 3; ++axis)
      for (int sign : {-1, 1})
        sides.push_back(
          canonicalSide(axis, sign, 2 * axis + int(sign > 0), low, high, tolerance, maximumCells));
    const auto reference = construct(
      bare,
      low,
      high,
      std::move(sides),
      original.Orientation() == TopAbs_FORWARD ? Operation::Common : Operation::BoxMinusField,
      tolerance,
      512 * std::numeric_limits<double>::epsilon() * magnitude,
      copy.Policy);
    report.ReferenceSeamBound = reference.GeometricBound;
    Map aEdges, bEdges, aVertices, bVertices, aNative, bNative;
    TopExp::MapShapes(copy.Shape, TopAbs_EDGE, aEdges);
    TopExp::MapShapes(reference.Shape, TopAbs_EDGE, bEdges);
    TopExp::MapShapes(copy.Shape, TopAbs_VERTEX, aVertices);
    TopExp::MapShapes(reference.Shape, TopAbs_VERTEX, bVertices);
    TopExp::MapShapes(original, TopAbs_EDGE, aNative);
    TopExp::MapShapes(reference.Faces.back(), TopAbs_EDGE, bNative);
    if (aEdges.Extent() != bEdges.Extent() || aVertices.Extent() != bVertices.Extent()
        || copy.Faces.size() != reference.Faces.size() || copy.SeamCount != reference.SeamCount
        || copy.BoxEdgeRootCount != reference.BoxEdgeRootCount
        || copy.ExactCornerCount != reference.ExactCornerCount
        || copy.CornerContourEnds != reference.CornerContourEnds
        || copy.PlanarCornerContacts != reference.PlanarCornerContacts)
    {
      report.Reason = "reference_graph_counts_differ";
      return report;
    }
    const auto       ac = curves(aEdges), bc = curves(bEdges);
    std::vector<int> mapping(ac.size(), -1), directions(ac.size(), 0), inverse(bc.size(), -1);
    std::vector<int> vertexMap(size_t(aVertices.Extent()), -1),
      vertexInverse(size_t(bVertices.Extent()), -1);
    for (size_t i = 0; i < ac.size(); ++i)
    {
      for (size_t j = 0; j < bc.size(); ++j)
      {
        if (report.Comparisons >= maximumComparisons)
        {
          report.Reason = "comparison_budget_exhausted";
          return report;
        }
        ++report.Comparisons;
        if (aNative.Contains(aEdges(int(i) + 1)) != bNative.Contains(bEdges(int(j) + 1)))
          continue;
        for (int direction : {1, -1})
        {
          if (ac[i].Points.front().Distance(direction > 0 ? bc[j].Points.front()
                                                          : bc[j].Points.back())
                > report.MatchingTolerance
              || ac[i].Points.back().Distance(direction > 0 ? bc[j].Points.back()
                                                            : bc[j].Points.front())
                   > report.MatchingTolerance)
            continue;
          auto candidate = bc[j];
          if (direction < 0)
          {
            auto p = candidate.Points;
            std::reverse(p.begin(), p.end());
            candidate = curveRecord(std::move(p));
          }
          const double bound = curveBound(ac[i], candidate, reserve, report);
          if (bound > report.MatchingTolerance)
            continue;
          if (mapping[i] >= 0 || inverse[j] >= 0)
          {
            report.Reason = "ambiguous_or_nonbijective_curve_match";
            return report;
          }
          mapping[i]                 = int(j);
          directions[i]              = direction;
          report.CurveDeviationBound = std::max(report.CurveDeviationBound, bound);
        }
      }
      if (mapping[i] < 0)
      {
        report.Reason = "stored_curve_has_no_reference_match";
        return report;
      }
      const size_t j = size_t(mapping[i]);
      inverse[j]     = int(i);
      TopoDS_Vertex av[2], bv[2];
      TopExp::Vertices(TopoDS::Edge(aEdges(int(i) + 1)), av[0], av[1]);
      TopExp::Vertices(TopoDS::Edge(bEdges(int(j) + 1)), bv[0], bv[1]);
      if (directions[i] < 0)
        std::swap(bv[0], bv[1]);
      for (int endpoint = 0; endpoint < 2; ++endpoint)
      {
        const int a = aVertices.FindIndex(av[endpoint]) - 1,
                  b = bVertices.FindIndex(bv[endpoint]) - 1;
        if (a < 0 || b < 0 || (vertexMap[size_t(a)] >= 0 && vertexMap[size_t(a)] != b)
            || (vertexInverse[size_t(b)] >= 0 && vertexInverse[size_t(b)] != a)
            || BRep_Tool::Pnt(av[endpoint]).Distance(BRep_Tool::Pnt(bv[endpoint]))
                 > report.MatchingTolerance)
        {
          report.Reason = "shared_vertex_correspondence_differ";
          return report;
        }
        vertexMap[size_t(a)]     = b;
        vertexInverse[size_t(b)] = a;
      }
      ++report.MatchedEdges;
    }
    for (int vertex : vertexMap)
      if (vertex >= 0)
        ++report.MatchedVertices;
    if (report.MatchedVertices != aVertices.Extent())
    {
      report.Reason = "incomplete_vertex_bijection";
      return report;
    }
    std::vector<int> identity(bc.size()), forward(bc.size(), 1);
    for (size_t i = 0; i < identity.size(); ++i)
      identity[i] = int(i);
    if (caps(copy, aEdges, mapping, directions) != caps(reference, bEdges, identity, forward))
    {
      report.Reason = "cap_grouping_or_directed_wires_differ";
      return report;
    }
    const auto aw = wires(original, aEdges, mapping, directions);
    if (aw != wires(reference.Faces.back(), bEdges, identity, forward))
    {
      report.Reason = "implicit_directed_wires_differ";
      return report;
    }
    report.MatchedCaps  = int(copy.Faces.size()) - 1;
    report.MatchedWires = int(aw.size());
    // Parameter-by-parameter agreement is needed too: matching 3D edges alone
    // cannot establish consistency with a planar cap's p-curve direction.
    for (size_t cap = 0; cap + 1 < copy.Faces.size(); ++cap)
    {
      const auto& face  = copy.Faces[cap];
      const auto  plane = BRepAdaptor_Surface(face).Plane();
      for (TopExp_Explorer use(face, TopAbs_EDGE); use.More(); use.Next())
      {
        const auto      edge = TopoDS::Edge(use.Current());
        TopLoc_Location location;
        double          first, last, pfirst, plast;
        const auto      curve =
          occ::down_cast<Geom_BSplineCurve>(BRep_Tool::Curve(edge, location, first, last));
        const auto uv =
          occ::down_cast<Geom2d_BSplineCurve>(BRep_Tool::CurveOnSurface(edge, face, pfirst, plast));
        if (uv.IsNull() || typeid(*uv) != typeid(Geom2d_BSplineCurve) || uv->Degree() != 1
            || uv->IsRational() || uv->IsPeriodic() || uv->NbPoles() != curve->NbPoles()
            || uv->NbKnots() != curve->NbKnots() || pfirst != first || plast != last)
        {
          report.Reason = "unsupported_cap_parameterization";
          return report;
        }
        if (report.SampleTests + size_t(uv->NbPoles()) > 8000000)
        {
          report.Reason = "cap_parameter_budget_exhausted";
          return report;
        }
        for (int pole = 1; pole <= uv->NbPoles(); ++pole)
        {
          ++report.SampleTests;
          if (uv->Knot(pole) != curve->Knot(pole)
              || uv->Multiplicity(pole) != curve->Multiplicity(pole))
          {
            report.Reason = "cap_parameter_knots_differ";
            return report;
          }
          const auto   p = uv->Pole(pole);
          const double error =
            ElSLib::Value(p.X(), p.Y(), plane).Distance(curve->Pole(pole)) + reserve;
          if (!(error <= report.MatchingTolerance))
          {
            report.Reason = "cap_curve_parameter_mismatch";
            return report;
          }
          report.CapParameterDeviationBound = std::max(report.CapParameterDeviationBound, error);
        }
      }
    }
    if (report.ReferenceSeamBound + report.CurveDeviationBound + report.CapParameterDeviationBound
        > tolerance)
    {
      report.Reason = "combined_seam_budget_exhausted";
      return report;
    }
    const auto&  field = *static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
    const double h     = 4 * std::pow(2 * std::acos(-1.) / field.Period(), 2);
    const double gReserve =
      h * reserve + 1024 * std::numeric_limits<double>::epsilon() * field.LipschitzBound();
    for (size_t cap = 0; cap + 1 < copy.Faces.size(); ++cap)
    {
      const gp_Vec normal(BRepAdaptor_Surface(copy.Faces[cap]).Plane().Axis().Direction());
      for (TopExp_Explorer edge(copy.Faces[cap], TopAbs_EDGE); edge.More(); edge.Next())
      {
        if (!aNative.Contains(edge.Current()))
          continue;
        const int   id     = aEdges.FindIndex(edge.Current()) - 1;
        const auto& points = ac[size_t(id)].Points;
        for (size_t p = 1; p < points.size(); ++p)
        {
          gp_Vec       tangent(points[p - 1], points[p]);
          const double length = tangent.Magnitude();
          if (edge.Current().Orientation() == TopAbs_REVERSED)
            tangent.Reverse();
          const gp_Pnt mid((points[p - 1].XYZ() + points[p].XYZ()) * .5);
          const double lower = (tangent / length).Dot(normal.Crossed(field.Gradient(mid)))
                                 * (original.Orientation() == TopAbs_FORWARD ? 1 : -1)
                               - h * length * .5 - gReserve;
          if (!(lower > 0))
          {
            report.Reason = "seam_tangent_alignment_unresolved";
            return report;
          }
          report.MinimumTangentAlignment = report.MinimumTangentAlignment == 0
                                             ? lower
                                             : std::min(report.MinimumTangentAlignment, lower);
        }
      }
    }
    report.Resolved = true;
    report.Reason   = "matched";
  }
  catch (const Standard_NotImplemented& error)
  {
    report.Reason = error.what();
  }
  return report;
}

namespace
{
using RootKey = std::array<int, 2>;
using ArcKey  = std::array<int, 7>;
using CapPair = std::array<int, 3>; // side, retained cap, complement cap
using Segment = std::pair<int, int>;

int physicalEdge(const gp_Pnt& p, const gp_Pnt& low, const gp_Pnt& high, double reserve)
{
  int fixed[3] = {-1, -1, -1}, count = 0;
  for (int axis = 0; axis < 3; ++axis)
  {
    if (std::abs(p.Coord(axis + 1) - low.Coord(axis + 1)) <= reserve)
      fixed[axis] = 0;
    else if (std::abs(p.Coord(axis + 1) - high.Coord(axis + 1)) <= reserve)
      fixed[axis] = 1;
    if (fixed[axis] >= 0)
      ++count;
  }
  require(count != 3, "Trim attachment has an ambiguous box-corner anchor");
  if (count != 2)
    return -1;
  int axis = 0;
  while (fixed[axis] >= 0)
    ++axis;
  return 4 * axis + fixed[(axis + 1) % 3] + 2 * fixed[(axis + 2) % 3];
}

std::map<int, RootKey> rootLabels(const std::map<int, gp_Pnt>& points,
                                  const gp_Pnt&                low,
                                  const gp_Pnt&                high,
                                  double                       reserve,
                                  std::array<int, 12>&         counts,
                                  const BRepClass3d_ImplicitBox::CornerModel* corners = nullptr)
{
  std::array<std::vector<std::pair<double, int>>, 12> groups;
  std::map<int, RootKey> labels;
  std::set<int> seenCorners;
  for (const auto& entry : points)
  {
    if (corners)
    {
      int corner = 0;
      bool exact = true;
      for (int axis = 0; axis < 3; ++axis)
      {
        if (entry.second.Coord(axis + 1) == high.Coord(axis + 1))
          corner |= 1 << axis;
        else if (entry.second.Coord(axis + 1) != low.Coord(axis + 1))
          exact = false;
      }
      if (exact)
      {
        const auto& germ = corners->Corners[size_t(corner)];
        require(germ.ExactZero && germ.Kind == BRepClass3d_ImplicitBox::CornerKind::SplitCone
                  && seenCorners.insert(corner).second,
                "Trim anchor lacks a unique checked split-corner identity");
        labels.emplace(entry.first, RootKey{12 + corner, 0}); // Edges 0..11; corners 12..19.
        continue; // A corner is NOT one root on each of its three incident edges.
      }
    }
    const int edge = physicalEdge(entry.second, low, high, reserve);
    if (edge >= 0)
      groups[size_t(edge)].emplace_back(entry.second.Coord(edge / 4 + 1), entry.first);
  }
  for (int edge = 0; edge < 12; ++edge)
  {
    auto& group = groups[size_t(edge)];
    std::sort(group.begin(), group.end());
    counts[size_t(edge)] = int(group.size());
    for (size_t i = 0; i < group.size(); ++i)
    {
      require(i == 0 || group[i].first - group[i - 1].first > 4 * reserve,
              "Trim attachment root ordering precision exhausted");
      labels.emplace(group[i].second, RootKey{edge, int(i)});
    }
  }
  return labels;
}

int areaDirection(const std::vector<gp_Pnt>& points, int side, double reserve)
{
  require(points.size() >= 3, "Degenerate closed trim arc");
  const int   u = (side / 2 + 1) % 3 + 1, v = (side / 2 + 2) % 3 + 1;
  long double area     = 0;
  double      diameter = 0;
  const auto& origin   = points.front();
  for (size_t i = 1; i < points.size(); ++i)
  {
    const auto&       a  = points[i - 1];
    const auto&       b  = points[i];
    const long double au = static_cast<long double>(a.Coord(u)) - origin.Coord(u),
                      av = static_cast<long double>(a.Coord(v)) - origin.Coord(v),
                      bu = static_cast<long double>(b.Coord(u)) - origin.Coord(u),
                      bv = static_cast<long double>(b.Coord(v)) - origin.Coord(v);
    area += au * bv - av * bu;
    diameter = std::max(diameter, origin.Distance(b));
  }
  area *= side % 2 ? 1 : -1;
  const double error = points.size()
                       * (4096 * std::numeric_limits<double>::epsilon() * diameter * diameter
                          + 64 * reserve * diameter);
  require(std::abs(area) > error, "Closed trim winding area unresolved");
  return area > 0 ? 1 : -1;
}

ArcKey arcKey(const CapPair& pair, RootKey first, RootKey last, int& direction)
{
  if (last < first)
  {
    std::swap(first, last);
    direction = -1;
  }
  return {pair[0], pair[1], pair[2], first[0], first[1], last[0], last[1]};
}

std::vector<int> seamOwners(const Boundary& boundary, const Map& seams)
{
  std::vector<int> owners(size_t(seams.Extent()), -1);
  for (size_t cap = 0; cap + 1 < boundary.Faces.size(); ++cap)
    for (TopExp_Explorer it(boundary.Faces[cap], TopAbs_EDGE); it.More(); it.Next())
    {
      const int index = seams.FindIndex(it.Current()) - 1;
      if (index < 0)
        continue;
      require(owners[size_t(index)] < 0, "Native seam has multiple cap owners");
      owners[size_t(index)] = int(cap);
    }
  for (int owner : owners)
    require(owner >= 0, "Native seam has no retained cap owner");
  return owners;
}
} // namespace

BRepClass3d_ImplicitBox::TrimAttachment BRepClass3d_ImplicitBox::CheckTrimAttachment(
  const Boundary& source,
  int             resolution,
  int             maximumDepth,
  size_t          maximumCells,
  size_t          maximumTests)
{
  if (resolution < 4 || resolution > 64 || (resolution & (resolution - 1)) || maximumDepth < 0
      || maximumDepth > 16 || !maximumCells || maximumCells > 10000000 || !maximumTests
      || maximumTests > 10000000)
    throw Standard_DomainError("Invalid trim attachment settings");
  TrimAttachment     report;
  std::ostringstream record;
  BRepClass3d_ImplicitBoxIO::Write(source, record);
  std::istringstream input(record.str());
  Boundary           copy;
  BRepClass3d_ImplicitBoxIO::Read(input, copy);
  auto retained = CheckCapAttachment(copy, resolution, maximumDepth, maximumCells, maximumTests);
  report.CapRegionsAttached = retained.Resolved;
  if (!retained.Resolved)
  {
    report.Reason = "retained_caps_unresolved: " + retained.Reason;
    return report;
  }
  const auto& face = copy.Faces.back();
  const auto  data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  const auto  low = data->Minimum(), high = data->Maximum();
  double      magnitude = 1;
  for (int axis = 1; axis <= 3; ++axis)
    magnitude = std::max({magnitude, std::abs(low.Coord(axis)), std::abs(high.Coord(axis))});
  const double reserve = 8192 * std::numeric_limits<double>::epsilon() * magnitude;
  try
  {
    TopoDS_Face bare;
    BRep_Builder()
      .MakeImplicitFace(bare, data->Field(), low, high, TopLoc_Location(), data->Tolerance());
    std::vector<Side> sides;
    for (int axis = 0; axis < 3; ++axis)
      for (int sign : {-1, 1})
        sides.push_back(
          canonicalSide(axis, sign, 2 * axis + int(sign > 0), low, high, data->Tolerance(), 512));
    const auto opposite =
      construct(bare,
                low,
                high,
                std::move(sides),
                face.Orientation() == TopAbs_FORWARD ? Operation::BoxMinusField : Operation::Common,
                data->Tolerance(),
                512 * std::numeric_limits<double>::epsilon() * magnitude,
                copy.Policy);
    const auto complement =
      CheckCapAttachment(opposite, resolution, maximumDepth, maximumCells, maximumTests);
    report.ComplementCapsAttached = complement.Resolved;
    if (!complement.Resolved)
    {
      report.Reason = "complement_caps_unresolved: " + complement.Reason;
      return report;
    }
    const auto&       ref   = retained.Reference;
    const auto&       other = complement.Reference;
    require(ref.Surface.Policy == copy.Policy && other.Surface.Policy == copy.Policy
              && ref.Surface.Resolution == other.Surface.Resolution,
            "Dual material references use different policies or grids");
    const auto&       loops = ref.Surface.BoundaryLoops;
    std::set<Segment> boundaryEdges;
    for (const auto& loop : loops)
      for (size_t i = 1; i < loop.size(); ++i)
        boundaryEdges.insert(std::minmax(loop[i - 1], loop[i]));
    std::map<RootKey, int> originRoots;
    for (size_t i = 0; i < ref.Surface.GridEdges.size(); ++i)
      originRoots.emplace(ref.Surface.GridEdges[i], int(i));
    require(other.Surface.GridEdges.size() == originRoots.size(),
            "Dual reference root counts differ");
    std::vector<int> otherRoots;
    for (const auto& origin : other.Surface.GridEdges)
    {
      const auto found = originRoots.find(origin);
      require(found != originRoots.end(), "Dual reference root provenance differs");
      otherRoots.push_back(found->second);
    }
    auto owners = [&](const CapAttachment& attachment, const std::vector<int>* mapping) {
      std::map<Segment, int> result;
      for (size_t i = 0; i < attachment.Reference.Triangles.size(); ++i)
      {
        const int cap = attachment.ReferenceTriangleCaps[i];
        if (cap < 0)
          continue;
        const auto& t = attachment.Reference.Triangles[i];
        for (int j = 0; j < 3; ++j)
        {
          int a = t[size_t(j)], b = t[size_t((j + 1) % 3)];
          if (size_t(a) >= attachment.Reference.Surface.Nodes.size()
              || size_t(b) >= attachment.Reference.Surface.Nodes.size())
            continue;
          if (mapping)
          {
            a = (*mapping)[size_t(a)];
            b = (*mapping)[size_t(b)];
          }
          const Segment edge = std::minmax(a, b);
          if (!boundaryEdges.count(edge))
            continue;
          require(result.emplace(edge, cap).second, "Ambiguous reference contour cap owner");
        }
      }
      require(result.size() == boundaryEdges.size(), "Reference contour cap owner missing");
      return result;
    };
    const auto                        currentOwners = owners(retained, nullptr),
                                      otherOwners   = owners(complement, &otherRoots);
    std::vector<std::vector<CapPair>> labels;
    std::map<int, gp_Pnt>             refAnchors;
    for (const auto& loop : loops)
    {
      std::vector<CapPair> line;
      for (size_t i = 1; i < loop.size(); ++i)
      {
        const Segment edge = std::minmax(loop[i - 1], loop[i]);
        const int     a = currentOwners.at(edge), b = otherOwners.at(edge);
        const int     side = retained.ReferenceCapSides[size_t(a)];
        require(side == complement.ReferenceCapSides[size_t(b)], "Dual cap physical sides differ");
        line.push_back({side, a, b});
      }
      for (size_t i = 0; i < line.size(); ++i)
        if (line[i] != line[(i + line.size() - 1) % line.size()])
          refAnchors.emplace(loop[i], ref.Nodes[size_t(loop[i])]);
      labels.push_back(std::move(line));
    }
    std::array<int, 12> refCounts{}, sourceCounts{};
    CornerModel corners;
    const CornerModel* cornerModel = nullptr;
    if (copy.Policy == CornerPolicy::Regularized)
    {
      corners = CheckCornerGerms(face, data->Tolerance());
      require(corners.Resolved, "Trim attachment corner identities unresolved");
      cornerModel = &corners;
    }
    const auto referenceRoots = rootLabels(refAnchors, low, high, reserve, refCounts, cornerModel);
    Map                 aEdges, bEdges, aVertices;
    TopExp::MapShapes(face, TopAbs_EDGE, aEdges);
    TopExp::MapShapes(opposite.Faces.back(), TopAbs_EDGE, bEdges);
    TopExp::MapShapes(face, TopAbs_VERTEX, aVertices);
    std::map<int, gp_Pnt> sourcePoints;
    for (int i = 1; i <= aVertices.Extent(); ++i)
      sourcePoints.emplace(i - 1, BRep_Tool::Pnt(TopoDS::Vertex(aVertices(i))));
    const auto sourceRoots = rootLabels(sourcePoints, low, high, reserve, sourceCounts, cornerModel);
    require(refCounts == sourceCounts, "Native/reference ordered box-root counts differ");
    require(int(referenceRoots.size()) == copy.BoxEdgeRootCount + corners.SplitCount
              && referenceRoots.size() == sourceRoots.size(),
            "Native/reference box-root anchors differ");
    std::map<ArcKey, int>                            arcIds;
    std::vector<std::vector<std::pair<ArcKey, int>>> referenceKeys;
    for (size_t wire = 0; wire < loops.size(); ++wire)
    {
      const auto& loop  = loops[wire];
      const auto& line  = labels[wire];
      size_t      start = 0;
      while (start < line.size() && line[start] == line[(start + line.size() - 1) % line.size()])
        ++start;
      const bool closed = start == line.size();
      if (closed)
        start = 0;
      std::vector<std::pair<ArcKey, int>> key;
      size_t                              done = 0;
      while (done < line.size())
      {
        const size_t begin = (start + done) % line.size();
        size_t       span  = 1;
        while (done + span < line.size() && line[(begin + span) % line.size()] == line[begin])
          ++span;
        const size_t end       = (begin + span) % line.size();
        int          direction = 1;
        RootKey      a{-1, -1}, b{-1, -1};
        if (!closed)
        {
          require(referenceRoots.count(loop[begin]) && referenceRoots.count(loop[end]),
                  "Reference arc endpoint lacks root anchor");
          a = referenceRoots.at(loop[begin]);
          b = referenceRoots.at(loop[end]);
        }
        auto descriptor = arcKey(line[begin], a, b, direction);
        if (closed)
        {
          std::vector<gp_Pnt> path;
          for (size_t i = 0; i <= span; ++i)
            path.push_back(ref.Nodes[size_t(loop[(begin + i) % line.size()])]);
          direction = areaDirection(path, line[begin][0], reserve);
        }
        require(arcIds.emplace(descriptor, -1).second, "Nonunique reference contour identity");
        key.emplace_back(descriptor, direction);
        done += span;
      }
      referenceKeys.push_back(std::move(key));
    }
    int next = 0;
    for (auto& entry : arcIds)
      entry.second = next++;
    report.ReferenceArcCount  = next;
    report.ReferenceWireCount = int(loops.size());
    require(aEdges.Extent() == next && aEdges.Extent() == bEdges.Extent(),
            "Native/reference contour arc counts differ");
    const auto       ac = curves(aEdges), bc = curves(bEdges);
    const auto       aOwners = seamOwners(copy, aEdges), bOwners = seamOwners(opposite, bEdges);
    std::vector<int> matching(ac.size(), -1), directions(ac.size(), 1), oppositeUsed(bc.size(), 0),
      used(size_t(next), 0);
    Correspondence samples;
    for (size_t i = 0; i < ac.size(); ++i)
    {
      int otherEdge = -1;
      for (size_t j = 0; j < bc.size(); ++j)
      {
        require(report.Comparisons < maximumTests, "Dual contour comparison budget exhausted");
        ++report.Comparisons;
        for (int dir : {1, -1})
        {
          if (ac[i].Points.front().Distance(dir > 0 ? bc[j].Points.front() : bc[j].Points.back())
                > retained.MatchingTolerance
              || ac[i].Points.back().Distance(dir > 0 ? bc[j].Points.back() : bc[j].Points.front())
                   > retained.MatchingTolerance)
            continue;
          auto candidate = bc[j];
          if (dir < 0)
          {
            auto points = candidate.Points;
            std::reverse(points.begin(), points.end());
            candidate = curveRecord(std::move(points));
          }
          const double deviation = curveBound(ac[i], candidate, reserve, samples);
          report.CurveSamples    = samples.SampleTests;
          if (deviation > retained.MatchingTolerance)
            continue;
          require(otherEdge < 0 && !oppositeUsed[j], "Ambiguous dual contour match");
          otherEdge = int(j);
        }
      }
      require(otherEdge >= 0, "Native contour lacks a complement trace");
      oppositeUsed[size_t(otherEdge)] = 1;
      const int a                     = retained.ReferenceCaps[size_t(aOwners[i])],
                b    = complement.ReferenceCaps[size_t(bOwners[size_t(otherEdge)])],
                side = retained.ReferenceCapSides[size_t(a)];
      require(side == complement.ReferenceCapSides[size_t(b)], "Native dual contour sides differ");
      TopoDS_Vertex first, last;
      TopExp::Vertices(TopoDS::Edge(aEdges(int(i) + 1)), first, last);
      RootKey f{-1, -1}, l{-1, -1};
      if (!first.IsSame(last))
      {
        const int fi = aVertices.FindIndex(first) - 1, li = aVertices.FindIndex(last) - 1;
        require(sourceRoots.count(fi) && sourceRoots.count(li),
                "Native contour endpoint lacks a box-root anchor");
        f = sourceRoots.at(fi);
        l = sourceRoots.at(li);
      }
      auto descriptor = arcKey({side, a, b}, f, l, directions[i]);
      if (first.IsSame(last))
        directions[i] = areaDirection(ac[i].Points, side, reserve);
      const auto found = arcIds.find(descriptor);
      require(found != arcIds.end(), "Native dual cap/root identity has no reference arc");
      require(!used[size_t(found->second)]++, "Nonbijective native contour identity");
      matching[i] = found->second;
    }
    report.CurveSamples      = samples.SampleTests;
    report.SeamReferenceArcs = matching;
    std::map<WireKey, int> referenceWires;
    for (size_t i = 0; i < referenceKeys.size(); ++i)
    {
      WireKey key;
      for (const auto& arc : referenceKeys[i])
        key.push_back((arcIds.at(arc.first) + 1) * arc.second);
      std::rotate(key.begin(), std::min_element(key.begin(), key.end()), key.end());
      require(referenceWires.emplace(key, int(i)).second, "Ambiguous reference directed wire");
    }
    std::vector<int> parent(ref.Surface.Nodes.size());
    for (size_t i = 0; i < parent.size(); ++i)
      parent[i] = int(i);
    auto root = [&](int i) {
      while (parent[size_t(i)] != i)
      {
        parent[size_t(i)] = parent[size_t(parent[size_t(i)])];
        i                 = parent[size_t(i)];
      }
      return i;
    };
    for (const auto& t : ref.Surface.Triangles)
    {
      parent[size_t(root(t[0]))] = root(t[1]);
      parent[size_t(root(t[1]))] = root(t[2]);
    }
    std::map<int, int> components;
    for (size_t i = 0; i < parent.size(); ++i)
      components.emplace(root(int(i)), 0);
    for (auto& entry : components)
      entry.second = report.ReferenceSurfaceComponentCount++;
    report.SurfaceNodeComponents.resize(parent.size());
    report.SurfaceComponentShells.assign(size_t(report.ReferenceSurfaceComponentCount), -1);
    report.SurfaceComponentRegions.assign(size_t(report.ReferenceSurfaceComponentCount), -1);
    report.SurfaceComponentEuler.assign(size_t(report.ReferenceSurfaceComponentCount), 0);
    for (size_t i = 0; i < parent.size(); ++i)
    {
      const int component             = components.at(root(int(i)));
      report.SurfaceNodeComponents[i] = component;
      ++report.SurfaceComponentEuler[size_t(component)];
    }
    std::set<Segment> surfaceEdges;
    for (size_t i = 0; i < ref.Surface.Triangles.size(); ++i)
    {
      const auto&  t     = ref.Surface.Triangles[i];
      const size_t c     = size_t(report.SurfaceNodeComponents[size_t(t[0])]);
      const int    shell = ref.TriangleComponents[i], region = ref.ComponentRegions[size_t(shell)];
      require(report.SurfaceComponentShells[c] < 0 || report.SurfaceComponentShells[c] == shell,
              "Reference surface component spans multiple boundary shells");
      report.SurfaceComponentShells[c]  = shell;
      report.SurfaceComponentRegions[c] = region;
      ++report.SurfaceComponentEuler[c];
      for (int j = 0; j < 3; ++j)
        surfaceEdges.insert(std::minmax(t[size_t(j)], t[size_t((j + 1) % 3)]));
    }
    for (const auto& edge : surfaceEdges)
      --report.SurfaceComponentEuler[size_t(report.SurfaceNodeComponents[size_t(edge.first)])];
    std::set<int> usedWires;
    for (const auto& key : wires(face, aEdges, matching, directions, false))
    {
      const auto found = referenceWires.find(key);
      require(found != referenceWires.end() && usedWires.insert(found->second).second,
              "Native/reference directed wire cycles differ");
      report.WireReferenceLoops.push_back(found->second);
      report.WireSurfaceComponents.push_back(
        components.at(root(loops[size_t(found->second)].front())));
    }
    require(usedWires.size() == loops.size(), "Native/reference wire bijection incomplete");
    report.Resolved     = true;
    report.Reason       = "trim_wires_attached";
    report.RetainedCaps = std::move(retained);
  }
  catch (const Standard_NotImplemented& error)
  {
    report.Reason = error.what();
  }
  return report;
}
