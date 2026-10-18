// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRep_ImplicitSphereOps.hxx>
#include <BRep_TImplicitSpherePatch.hxx>
#include <BRep_Builder.hxx>
#include <Geom_Circle.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <Poly_Triangle.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <typeinfo>

IMPLEMENT_STANDARD_RTTIEXT(BRep_TImplicitSpherePatch, BRep_TImplicitFace)

namespace
{
constexpr double Pi = 3.1415926535897932384626433832795;

void require(bool valid, const char* message)
{
  if (!valid)
    throw Standard_DomainError(message);
}

gp_Pnt corner(const occ::handle<Geom_ImplicitSphere>& field, double sign, double tol)
{
  require(!field.IsNull(), "Null spherical patch field");
  const auto&  s = *field;
  const double d = sign * (s.Radius() + tol);
  return s.Center().Translated(gp_Vec(d, d, d));
}

struct Sphere
{
  gp_Pnt Center;
  double Radius, Tolerance;
};

Sphere sphere(const TopoDS_Face& face)
{
  auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  if (data.IsNull() || !data->ExcludedDomains().empty() || typeid(*data) != typeid(BRep_TImplicitFace)
      || face.Orientation() != TopAbs_FORWARD || TopoDS_Iterator(face).More()
      || typeid(*data->Field()) != typeid(Geom_ImplicitSphere))
    throw Standard_NotImplemented("Native sphere operation requires an untrimmed forward sphere");
  const auto* s = static_cast<const Geom_ImplicitSphere*>(data->Field().get());
  for (int axis = 1; axis <= 3; ++axis)
    require(data->Minimum().Coord(axis) <= s->Center().Coord(axis) - s->Radius()
              && data->Maximum().Coord(axis) >= s->Center().Coord(axis) + s->Radius(),
            "Sphere domain does not enclose the complete surface");
  const gp_Trsf t = (face.Location() * data->Location()).Transformation();
  const double  r = s->Radius() * std::abs(t.ScaleFactor());
  const gp_Pnt  c = s->Center().Transformed(t);
  require(std::isfinite(r) && r > 0 && std::isfinite(c.X()) && std::isfinite(c.Y())
            && std::isfinite(c.Z()) && std::isfinite(data->Tolerance()) && data->Tolerance() > 0,
          "Non-finite placed sphere");
  return {c, r, data->Tolerance()};
}

TopoDS_Face full(const Sphere& s)
{
  auto        field = occ::handle<Geom_ImplicitSphere>(new Geom_ImplicitSphere(s.Center, s.Radius));
  TopoDS_Face face;
  BRep_Builder().MakeImplicitFace(face,
                                  field,
                                  corner(field, -1, s.Tolerance),
                                  corner(field, 1, s.Tolerance),
                                  TopLoc_Location(),
                                  s.Tolerance);
  return face;
}

double signedDistance(const gp_Pln& plane, const gp_Pnt& p)
{
  return gp_Vec(plane.Location(), p).Dot(gp_Vec(plane.Axis().Direction()));
}

TopoDS_Edge circleEdge(const gp_Pnt& c, const gp_Dir& n, double radius, double tolerance)
{
  BRep_Builder builder;
  auto         curve = occ::handle<Geom_Circle>(new Geom_Circle(gp_Ax2(c, n), radius));
  TopoDS_Edge  edge;
  builder.MakeEdge(edge, curve, tolerance);
  builder.Range(edge, 0, 2 * Pi);
  TopoDS_Vertex vertex;
  builder.MakeVertex(vertex, curve->Value(0), tolerance);
  builder.Add(edge, vertex.Oriented(TopAbs_FORWARD));
  builder.Add(edge, vertex.Oriented(TopAbs_REVERSED));
  builder.UpdateVertex(vertex, 0, edge, tolerance);
  edge.Closed(true);
  return edge;
}

TopoDS_Face patch(const Sphere& s, const gp_Pln& plane, const TopoDS_Edge& edge, bool reverseEdge)
{
  auto        field = occ::handle<Geom_ImplicitSphere>(new Geom_ImplicitSphere(s.Center, s.Radius));
  TopoDS_Face face;
  face.TShape(new BRep_TImplicitSpherePatch(field, plane, s.Tolerance));
  face.Orientation(TopAbs_FORWARD);
  BRep_Builder builder;
  TopoDS_Wire  wire;
  builder.MakeWire(wire);
  builder.Add(wire, edge.Oriented(reverseEdge ? TopAbs_REVERSED : TopAbs_FORWARD));
  wire.Closed(true);
  builder.Add(face, wire);
  return face;
}
} // namespace

BRep_TImplicitSpherePatch::BRep_TImplicitSpherePatch(const occ::handle<Geom_ImplicitSphere>& field,
                                                     const gp_Pln&                           plane,
                                                     double                                  tol)
    : BRep_TImplicitFace(field, corner(field, -1, tol), corner(field, 1, tol), tol),
      myPlane(plane)
{
  if (typeid(*field) != typeid(Geom_ImplicitSphere))
    throw Standard_NotImplemented("Patch requires exact sphere field");
  require(std::isfinite(plane.Location().X()) && std::isfinite(plane.Location().Y())
            && std::isfinite(plane.Location().Z()),
          "Non-finite trim plane");
  require(std::abs(signedDistance(plane, field->Center())) < field->Radius() - tol,
          "Native sphere trim requires a transverse circle away from tangency");
}

occ::handle<TopoDS_TShape> BRep_TImplicitSpherePatch::EmptyCopy() const
{
  auto                                   field = occ::down_cast<Geom_ImplicitSphere>(Field());
  occ::handle<BRep_TImplicitSpherePatch> copy =
    new BRep_TImplicitSpherePatch(field, myPlane, Tolerance());
  copy->Location(Location());
  return copy;
}

void BRep_TImplicitSpherePatch::DumpJson(Standard_OStream& out, int depth) const {
  OCCT_DUMP_TRANSIENT_CLASS_BEGIN(out) OCCT_DUMP_BASE_CLASS(out, depth, BRep_TImplicitFace)
    OCCT_DUMP_FIELD_VALUES_DUMPED(out, depth, &myPlane)}

TopoDS_Face BRep_ImplicitSphereOps::Trim(const TopoDS_Face& input, const gp_Pln& plane)
{
  const Sphere s = sphere(input);
  const double h = -signedDistance(plane, s.Center);
  require(std::isfinite(h) && std::abs(h) < s.Radius - s.Tolerance,
          "Trim is non-transverse or within tolerance of tangency");
  const gp_Dir n      = plane.Axis().Direction();
  const gp_Pnt c      = s.Center.Translated(gp_Vec(n) * h);
  const double radius = std::sqrt((s.Radius - h) * (s.Radius + h));
  require(std::isfinite(radius) && radius > 0, "Unrepresentable intersection circle");
  return patch(s, gp_Pln(c, n), circleEdge(c, n, radius, s.Tolerance), true);
}

BRep_ImplicitSphereOps::Result BRep_ImplicitSphereOps::Boolean(const TopoDS_Face& left,
                                                               const TopoDS_Face& right,
                                                               Operation          op)
{
  require(op == Operation::Common || op == Operation::Fuse || op == Operation::Cut,
          "Unknown sphere Boolean operation");
  Sphere       a = sphere(left), b = sphere(right);
  const double tol = std::max(a.Tolerance, b.Tolerance);
  a.Tolerance = b.Tolerance = tol;
  const double d            = a.Center.Distance(b.Center);
  require(std::isfinite(d), "Non-finite sphere separation");
  require(std::abs(d - (a.Radius + b.Radius)) > tol
            && std::abs(d - std::abs(a.Radius - b.Radius)) > tol,
          "Tangency/coincidence or tolerance-ambiguous sphere Boolean");
  Result       result;
  BRep_Builder builder;
  builder.MakeCompound(result.Boundary);
  auto add = [&](TopoDS_Face face, int source, bool reversed = false) {
    if (reversed)
      face.Reverse();
    builder.Add(result.Boundary, face);
    result.Faces.push_back(face);
    result.Sources.push_back(source);
  };
  if (d > a.Radius + b.Radius)
  {
    if (op != Operation::Common)
      add(full(a), 0);
    if (op == Operation::Fuse)
      add(full(b), 1);
    return result;
  }
  if (d < std::abs(a.Radius - b.Radius))
  {
    const bool aInside = a.Radius < b.Radius;
    if (op == Operation::Common)
      add(full(aInside ? a : b), aInside ? 0 : 1);
    else if (op == Operation::Fuse)
      add(full(aInside ? b : a), aInside ? 1 : 0);
    else if (!aInside)
    {
      add(full(a), 0);
      add(full(b), 1, true);
    }
    return result;
  }
  const gp_Dir n(gp_Vec(a.Center, b.Center));
  // Scale the computation to avoid squaring large model lengths.
  const double x            = 0.5 * (d + (a.Radius - b.Radius) * ((a.Radius + b.Radius) / d));
  const gp_Pnt c            = a.Center.Translated(gp_Vec(n) * x);
  const double circleRadius = std::sqrt((a.Radius - x) * (a.Radius + x));
  require(std::isfinite(circleRadius) && circleRadius > tol,
          "Intersection circle is unresolved at the requested tolerance");
  const auto   edge = circleEdge(c, n, circleRadius, tol);
  const gp_Pln positive(c, n), negative(c, n.Reversed());
  if (op == Operation::Common)
  {
    add(patch(a, negative, edge, false), 0);
    add(patch(b, positive, edge, true), 1);
  }
  else if (op == Operation::Fuse)
  {
    add(patch(a, positive, edge, true), 0);
    add(patch(b, negative, edge, false), 1);
  }
  else
  {
    add(patch(a, positive, edge, true), 0);
    add(patch(b, positive, edge, true), 1, true);
  }
  return result;
}

BRep_ImplicitSphereOps::Mesh BRep_ImplicitSphereOps::Tessellate(const TopoDS_Face& face,
                                                                double             deflection)
{
  auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  if (data.IsNull() || !data->ExcludedDomains().empty() || typeid(*data->Field()) != typeid(Geom_ImplicitSphere)
      || (typeid(*data) != typeid(BRep_TImplicitFace)
          && typeid(*data) != typeid(BRep_TImplicitSpherePatch)))
    throw Standard_NotImplemented("Analytic mesher supports only native sphere/plane patches");
  require(face.Orientation() == TopAbs_FORWARD || face.Orientation() == TopAbs_REVERSED,
          "Unsupported mesh orientation");
  const auto*  field = static_cast<const Geom_ImplicitSphere*>(data->Field().get());
  const double r     = field->Radius();
  const double scale =
    std::abs((face.Location() * data->Location()).Transformation().ScaleFactor());
  require(std::isfinite(deflection) && deflection > 0 && std::isfinite(scale) && scale > 0,
          "Invalid mesh geometric deflection/placement");
  auto trimmed = occ::down_cast<BRep_TImplicitSpherePatch>(data);
  if (!trimmed.IsNull())
    require(data->Location().IsIdentity(),
            "Patch geometry location must remain identity; move its face instance");
  if (trimmed.IsNull())
  {
    TopoDS_Face forward = face;
    forward.Orientation(TopAbs_FORWARD);
    (void)sphere(forward);
  }
  gp_Dir axis(0, 0, 1);
  double thetaMax = Pi;
  if (!trimmed.IsNull())
  {
    axis           = trimmed->Plane().Axis().Direction().Reversed();
    const double h = -signedDistance(trimmed->Plane(), field->Center());
    thetaMax       = std::acos(std::clamp(-h / r, -1.0, 1.0));
  }
  const auto transform = (face.Location() * data->Location()).Transformation();
  for (int row = 1; row <= 3; ++row)
    for (int col = 1; col <= 4; ++col)
      require(std::isfinite(transform.Value(row, col)), "Non-finite mesh placement");
  const double magnitude = std::max({std::abs(field->Center().X()) * scale,
                                     std::abs(field->Center().Y()) * scale,
                                     std::abs(field->Center().Z()) * scale,
                                     r * scale,
                                     std::abs(transform.Value(1, 4)),
                                     std::abs(transform.Value(2, 4)),
                                     std::abs(transform.Value(3, 4))});
  // Refuse numerically unresolved coordinates rather than advertise a tiny
  // analytic bound at enormous model offsets. This is a floating-point safety
  // reserve, not an interval-arithmetic certificate for rounded computations.
  const double reserve = 64 * std::numeric_limits<double>::epsilon() * magnitude;
  require(std::isfinite(reserve) && deflection > reserve,
          "Deflection below coordinate precision reserve");
  const double step = std::min(0.5, std::sqrt(((deflection - reserve) / scale) / (2 * r)));
  require(std::isfinite(step) && step > 0, "Unrepresentable mesh deflection");
  const double nPhi   = std::max(3.0, std::ceil(2 * Pi / step));
  const double nTheta = std::max(2.0, std::ceil(thetaMax / step));
  require(nPhi * nTheta < 119998, "Native sphere mesh exceeds vertex budget");
  const int  columns = int(nPhi), rows = int(nTheta);
  const bool whole     = trimmed.IsNull();
  const int  rings     = whole ? rows - 1 : rows;
  const int  vertices  = 1 + rings * columns + (whole ? 1 : 0);
  const int  triangles = columns * (1 + 2 * (rings - 1) + (whole ? 1 : 0));
  require(vertices <= 120000 && triangles <= 240000, "Native sphere mesh budget exceeded");
  auto mesh = occ::handle<Poly_Triangulation>(new Poly_Triangulation(vertices, triangles, false));
  const gp_Ax2 frame(field->Center(), axis);
  const gp_Vec z(axis), u(frame.XDirection()), v(frame.YDirection());
  mesh->SetNode(1, field->Center().Translated(z * r));
  auto node = [columns](int ring, int col) { return 2 + (ring - 1) * columns + col % columns; };
  for (int ring = 1; ring <= rings; ++ring)
    for (int col = 0; col < columns; ++col)
    {
      const double theta = thetaMax * ring / rows, phi = 2 * Pi * col / columns;
      mesh->SetNode(
        node(ring, col),
        field->Center().Translated(
          r * (z * std::cos(theta) + (u * std::cos(phi) + v * std::sin(phi)) * std::sin(theta))));
    }
  if (whole)
    mesh->SetNode(vertices, field->Center().Translated(z * -r));
  int  index    = 1;
  auto triangle = [&](int a, int b, int c) {
    if (face.Orientation() == TopAbs_REVERSED)
      std::swap(b, c);
    mesh->SetTriangle(index++, Poly_Triangle(a, b, c));
  };
  for (int col = 0; col < columns; ++col)
  {
    triangle(1, node(1, col), node(1, col + 1));
    for (int ring = 1; ring < rings; ++ring)
    {
      triangle(node(ring, col), node(ring + 1, col), node(ring + 1, col + 1));
      triangle(node(ring, col), node(ring + 1, col + 1), node(ring, col + 1));
    }
    if (whole)
      triangle(node(rings, col), vertices, node(rings, col + 1));
  }
  // Spherical-map second derivatives have norm <= r. Barycentric Taylor
  // remainder <= r/2 * (deltaTheta + deltaPhi)^2, in both directions between
  // the parameter triangle and its affine mesh triangle. Includes trim edge.
  const double delta = thetaMax / rows + 2 * Pi / columns;
  const double bound = scale * (r * 0.5 * delta * delta) + reserve;
  require(std::isfinite(bound) && bound <= deflection * (1 + 1.e-12), "Mesh bound overflow");
  mesh->Deflection(bound);
  return {mesh, bound};
}
