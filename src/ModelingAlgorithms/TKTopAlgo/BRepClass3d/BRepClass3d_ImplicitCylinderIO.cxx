// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitCylinderIO.hxx>
#include <BRep_ImplicitIO.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_TImplicitCylinderPatch.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_Plane.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Line.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom2d_Circle.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <gp_Lin.hxx>
#include <gp_Lin2d.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <Standard_DomainError.hxx>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <typeinfo>

namespace
{
using IO = BRepClass3d_ImplicitCylinderIO;
using Map = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
constexpr size_t MaxBytes = 16 * 1024 * 1024, MaxMath = 65536, MaxPoles = 65536;
void check(bool ok, const char* reason) { if (!ok) throw Standard_DomainError(reason); }
double real(std::istream& in) { double v; check(bool(in >> v) && std::isfinite(v), "cylinder_io_real"); return v; }
int integer(std::istream& in, int low, int high)
{ int v; check(bool(in >> v) && v >= low && v <= high, "cylinder_io_integer"); return v; }
void word(std::istream& in, const char* expected)
{ std::string token; check(bool(in >> token) && token == expected, "cylinder_io_token"); }
gp_Pnt point(std::istream& in) { const double x = real(in), y = real(in), z = real(in); return gp_Pnt(x, y, z); }
gp_Pnt2d point2(std::istream& in) { const double x = real(in), y = real(in); return gp_Pnt2d(x, y); }
gp_Dir direction(std::istream& in)
{ const auto p = point(in); check(std::abs(p.Distance(gp_Pnt()) - 1) < 1.e-12, "cylinder_io_direction"); return gp_Dir(p.XYZ()); }
gp_Dir2d direction2(std::istream& in)
{ const auto p = point2(in); check(std::abs(p.Distance(gp_Pnt2d()) - 1) < 1.e-12, "cylinder_io_direction2"); return gp_Dir2d(p.XY()); }
gp_Ax2 frame(std::istream& in)
{ const auto p = point(in); const auto z = direction(in), x = direction(in);
  check(std::abs(z.Dot(x)) < 1.e-12, "cylinder_io_frame"); return gp_Ax2(p, z, x); }
void put(std::ostream& out, const gp_Pnt& p) { out << p.X() << ' ' << p.Y() << ' ' << p.Z() << ' '; }
void put(std::ostream& out, const gp_Pnt2d& p) { out << p.X() << ' ' << p.Y() << ' '; }
void put(std::ostream& out, const gp_Dir& d) { put(out, gp_Pnt(d.XYZ())); }
void put(std::ostream& out, const gp_Dir2d& d) { put(out, gp_Pnt2d(d.XY())); }
void put(std::ostream& out, const gp_Ax2& a) { put(out, a.Location()); put(out, a.Direction()); put(out, a.XDirection()); }
void put(std::ostream& out, const gp_Ax3& a)
{ check(a.Direct(), "cylinder_io_indirect_surface"); put(out, gp_Ax2(a.Location(), a.Direction(), a.XDirection())); }
int orientation(const TopoDS_Shape& s)
{ check(s.Orientation() == TopAbs_FORWARD || s.Orientation() == TopAbs_REVERSED, "cylinder_io_orientation");
  return s.Orientation() == TopAbs_REVERSED ? 1 : 0; }
TopoDS_Wire effectiveWire(const TopoDS_Face& face)
{ TopoDS_Iterator it(face); check(it.More() && it.Value().ShapeType() == TopAbs_WIRE, "cylinder_io_wire_count");
  auto wire = TopoDS::Wire(it.Value()); it.Next(); check(!it.More(), "cylinder_io_wire_count"); return wire; }
std::vector<TopoDS_Edge> ordered(const TopoDS_Wire& wire)
{ std::vector<TopoDS_Edge> edges; for (TopoDS_Iterator it(wire); it.More(); it.Next())
  { check(it.Value().ShapeType() == TopAbs_EDGE, "cylinder_io_wire_child"); edges.push_back(TopoDS::Edge(it.Value())); }
  if (wire.Orientation() == TopAbs_REVERSED) std::reverse(edges.begin(), edges.end());
  return edges; }
void pcWrite(std::ostream& out, const occ::handle<Geom2d_Curve>& pc)
{
  if (auto c = occ::down_cast<Geom2d_Circle>(pc); !c.IsNull())
  { out << "C "; put(out, c->Location()); put(out, c->XAxis().Direction()); put(out, c->YAxis().Direction()); out << c->Radius() << ' '; }
  else if (auto l = occ::down_cast<Geom2d_Line>(pc); !l.IsNull())
  { out << "L "; put(out, l->Lin2d().Location()); put(out, l->Lin2d().Direction()); }
  else
  { auto b = occ::down_cast<Geom2d_BSplineCurve>(pc);
    check(!b.IsNull() && b->Degree() == 1 && !b->IsPeriodic() && b->NbPoles() <= int(MaxPoles)
      && b->NbPoles() == b->NbKnots(), "cylinder_io_pcurve_type"); out << "B " << b->NbPoles() << ' ';
    for (int i = 1; i <= b->NbPoles(); ++i)
    { check(b->Knot(i) == i - 1 && b->Multiplicity(i) == (i == 1 || i == b->NbPoles() ? 2 : 1), "cylinder_io_pcurve_knots"); put(out, b->Pole(i)); }
  }
}
occ::handle<Geom2d_Curve> pcRead(std::istream& in)
{
  std::string kind; check(bool(in >> kind), "cylinder_io_pcurve");
  if (kind == "C")
  { const auto p = point2(in); const auto x = direction2(in), y = direction2(in); const double r = real(in);
    check(r > 0 && std::abs(x.Dot(y)) < 1.e-12, "cylinder_io_circle2"); return new Geom2d_Circle(gp_Ax22d(p, x, y), r); }
  if (kind == "L") { const auto p = point2(in); const auto d = direction2(in); return new Geom2d_Line(p, d); }
  check(kind == "B", "cylinder_io_pcurve_kind"); const int n = integer(in, 2, int(MaxPoles));
  NCollection_Array1<gp_Pnt2d> p(1, n); NCollection_Array1<double> k(1, n); NCollection_Array1<int> m(1, n);
  for (int i = 1; i <= n; ++i) { p.SetValue(i, point2(in)); k.SetValue(i, i - 1); m.SetValue(i, i == 1 || i == n ? 2 : 1); }
  return new Geom2d_BSplineCurve(p, k, m, 1);
}
void validate(const IO::StoredSolid& value)
{
  const auto& s = value.Parameters;
  for (double number : {s.Radius, s.HeightMinimum, s.HeightMaximum, s.DistanceTolerance,
      s.AngleMinimum, s.AngleMaximum, s.Frame.Location().X(), s.Frame.Location().Y(), s.Frame.Location().Z(),
      value.StoredBoundaryBound, value.StoredVolumeLower, value.StoredVolumeUpper, value.StoredVolumeTolerance})
    check(std::isfinite(number), "cylinder_io_nonfinite_metadata");
  check(!value.Shape.IsNull() && value.Shape.Location().IsIdentity() && orientation(value.Shape) == 0
    && value.Faces.size() == 3 && (value.StoredAxialDirection == 1 || value.StoredAxialDirection == -1)
    && s.Radius > 0 && s.HeightMaximum > s.HeightMinimum && s.DistanceTolerance > 0
    && s.AngleMaximum == s.AngleMinimum + 2 * std::acos(-1.) && std::abs(s.AngleMinimum) < 1.e6,
    "cylinder_io_recipe");
  const auto source = occ::down_cast<BRep_TImplicitFace>(value.Primitive.TShape());
  const auto data = occ::down_cast<BRep_TImplicitCylinderPatch>(value.Faces[2].TShape());
  check(!source.IsNull() && typeid(*source) == typeid(BRep_TImplicitFace) && !data.IsNull()
    && !source->Field().IsNull() && !data->Field().IsNull()
    && typeid(*source->Field()) == typeid(Geom_ImplicitGyroid) && typeid(*data->Field()) == typeid(Geom_ImplicitGyroid)
    && !TopoDS_Iterator(value.Primitive).More() && source->ExcludedDomains().empty() && data->ExcludedDomains().empty(), "cylinder_io_field");
  const auto* a = static_cast<const Geom_ImplicitGyroid*>(source->Field().get());
  const auto* b = static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  check(a->Period() == b->Period() && a->IsoValue() == b->IsoValue()
    && value.Negative == (value.Faces[2].Orientation() == TopAbs_FORWARD), "cylinder_io_field_identity");
  for (int axis = 1; axis <= 3; ++axis)
    check(source->Minimum().Coord(axis) == data->Minimum().Coord(axis)
      && source->Maximum().Coord(axis) == data->Maximum().Coord(axis), "cylinder_io_field_box");
  for (const auto& pair : {std::make_pair(value.Primitive.Location(), value.Faces[2].Location()),
                           std::make_pair(source->Location(), data->Location())})
    for (int r = 1; r <= 3; ++r) for (int c = 1; c <= 4; ++c)
      check(pair.first.Transformation().Value(r,c) == pair.second.Transformation().Value(r,c), "cylinder_io_field_placement");
  const auto transform = (value.Primitive.Location() * source->Location()).Transformation();
  const double scale = std::abs(transform.ScaleFactor());
  const double reserve = 32768 * std::numeric_limits<double>::epsilon()
    * (1 + s.Radius + std::abs(s.HeightMinimum) + std::abs(s.HeightMaximum) + s.Frame.Location().Distance(gp_Pnt()));
  check(data->Frame().Location().Transformed(transform).Distance(s.Frame.Location()) <= reserve
    && std::abs(data->Radius() * scale - s.Radius) <= reserve
    && std::abs(data->HeightMinimum() * scale - s.HeightMinimum) <= reserve
    && std::abs(data->HeightMaximum() * scale - s.HeightMaximum) <= reserve
    && gp_Dir(gp_Vec(data->Frame().Direction()).Transformed(transform)).Dot(s.Frame.Direction()) > 1 - 1.e-12,
    "cylinder_io_mask");
  check(std::isfinite(value.StoredBoundaryBound) && value.StoredBoundaryBound > 0
    && value.StoredBoundaryBound <= s.DistanceTolerance && value.StoredVolumeLower > 0
    && std::isfinite(value.StoredVolumeUpper) && value.StoredVolumeUpper >= value.StoredVolumeLower
    && value.StoredVolumeTolerance > 0 && value.StoredVolumeUpper - value.StoredVolumeLower <= value.StoredVolumeTolerance,
    "cylinder_io_stored_metadata");
  check(!occ::down_cast<Geom_Plane>(BRep_Tool::Surface(value.Faces[0])).IsNull()
    && !occ::down_cast<Geom_CylindricalSurface>(BRep_Tool::Surface(value.Faces[1])).IsNull()
    && value.Faces[0].Location().IsIdentity() && value.Faces[1].Location().IsIdentity(), "cylinder_io_surfaces");
  check(BRepCheck_Analyzer(value.Faces[0]).IsValid() && BRepCheck_Analyzer(value.Faces[1]).IsValid(), "cylinder_io_conventional_faces");
  TopoDS_Iterator shell(value.Shape); check(shell.More() && shell.Value().ShapeType() == TopAbs_SHELL, "cylinder_io_shell");
  const auto shape = TopoDS::Shell(shell.Value()); shell.Next(); check(!shell.More(), "cylinder_io_shell_count");
  BRepClass3d_ImplicitCylinder::ValidateStoredShellTopology(shape, value.Faces);
}
}

void BRepClass3d_ImplicitCylinderIO::Write(const BRepClass3d_ImplicitCylinder::NativeSolid& solid, Standard_OStream& output)
{
  check(solid.SolidValidated && solid.Volume.Verified && solid.Source.BoundaryVerified, "cylinder_io_requires_qualified_solid");
  const auto& shell = solid.Source.Source.Source;
  const auto& recipe = shell.Source.Reference.Source.Boundary.Source.Boundary.Source.Regions;
  StoredSolid value; value.Shape = solid.Solid; value.Faces = shell.Faces; value.Primitive = recipe.Field;
  value.Parameters = recipe.Parameters; value.Negative = recipe.Negative;
  value.StoredAxialDirection = shell.Source.Reference.Source.Direction;
  value.StoredBoundaryBound = solid.Source.BoundaryBound; value.StoredVolumeLower = solid.Volume.Lower;
  value.StoredVolumeUpper = solid.Volume.Upper; value.StoredVolumeTolerance = solid.Volume.AbsoluteTolerance;
  validate(value);
  std::ostringstream primitive, math; BRep_ImplicitIO::Write(value.Primitive, primitive);
  auto bare = value.Faces[2]; bare.TShape(bare.TShape()->EmptyCopy()); BRep_ImplicitIO::Write(bare, math);
  const auto p = primitive.str(), m = math.str(); check(p.size() <= MaxMath && m.size() <= MaxMath, "cylinder_io_math_size");
  std::ostringstream out; out.imbue(std::locale::classic()); out << std::setprecision(std::numeric_limits<double>::max_digits10);
  out << "OCCT_IMPLICIT_CYLINDER_SOLID 1 " << p.size() << ' ' << m.size() << '\n' << p << m;
  out << "META "; put(out, value.Parameters.Frame);
  const auto& s = value.Parameters;
  out << s.Radius << ' ' << s.AngleMinimum << ' ' << s.AngleMaximum << ' ' << s.HeightMinimum << ' ' << s.HeightMaximum << ' '
      << s.DistanceTolerance << ' ' << s.MinimumCells << ' ' << s.MaximumCells << ' ' << value.Negative << ' ' << value.StoredAxialDirection << ' '
      << value.StoredBoundaryBound << ' ' << value.StoredVolumeLower << ' ' << value.StoredVolumeUpper << ' ' << value.StoredVolumeTolerance << '\n';
  Map vertices, edges; TopExp::MapShapes(value.Shape, TopAbs_VERTEX, vertices); TopExp::MapShapes(value.Shape, TopAbs_EDGE, edges);
  check(vertices.Size() == 3 && edges.Size() == 4, "cylinder_io_graph_counts");
  out << "VERTICES 3\n";
  for (const auto& v : vertices) { check(v.Location().IsIdentity(), "cylinder_io_vertex_placement");
    put(out, BRep_Tool::Pnt(TopoDS::Vertex(v))); out << BRep_Tool::Tolerance(TopoDS::Vertex(v)) << '\n'; }
  out << "EDGES 4\n";
  for (const auto& shape : edges)
  {
    const auto edge = TopoDS::Edge(shape.Oriented(TopAbs_FORWARD)); TopoDS_Vertex a, b; TopExp::Vertices(edge, a, b);
    TopLoc_Location location; double first, last; const auto curve = BRep_Tool::Curve(edge, location, first, last);
    check(location.IsIdentity(), "cylinder_io_curve_placement");
    out << "EDGE " << first << ' ' << last << ' ' << BRep_Tool::Tolerance(edge) << ' ' << vertices.FindIndex(a) << ' ' << vertices.FindIndex(b) << ' ';
    if (auto c = occ::down_cast<Geom_Circle>(curve); !c.IsNull()) { out << "C "; put(out, c->Position()); out << c->Radius() << ' '; }
    else if (auto l = occ::down_cast<Geom_Line>(curve); !l.IsNull()) { out << "L "; put(out, l->Lin().Location()); put(out, l->Lin().Direction()); }
    else { const auto b = occ::down_cast<Geom_BSplineCurve>(curve);
      check(!b.IsNull() && b->Degree() == 1 && !b->IsPeriodic() && b->NbPoles() == b->NbKnots() && b->NbPoles() <= int(MaxPoles), "cylinder_io_curve_type");
      out << "B " << b->NbPoles() << ' ';
      for (int i = 1; i <= b->NbPoles(); ++i) { check(b->Knot(i) == i - 1 && b->Multiplicity(i) == (i == 1 || i == b->NbPoles() ? 2 : 1), "cylinder_io_curve_knots"); put(out, b->Pole(i)); }
    }
    out << '\n';
  }
  for (int f = 0; f < 3; ++f)
  {
    const auto& face = value.Faces[size_t(f)]; out << "FACE " << f << ' ' << orientation(face) << ' ' << BRep_Tool::Tolerance(face) << ' ';
    if (f == 0) put(out, occ::down_cast<Geom_Plane>(BRep_Tool::Surface(face))->Position());
    if (f == 1) { auto c = occ::down_cast<Geom_CylindricalSurface>(BRep_Tool::Surface(face)); put(out, c->Position()); out << c->Radius() << ' '; }
    const auto uses = ordered(effectiveWire(face)); check(uses.size() == size_t(f == 0 ? 2 : f == 1 ? 5 : 1), "cylinder_io_face_uses");
    out << uses.size() << '\n';
    for (const auto& edge : uses)
    {
      out << edges.FindIndex(edge) << ' ' << orientation(edge) << ' ';
      if (f != 2) { double first, last; const auto pc = BRep_Tool::CurveOnSurface(edge, face, first, last);
        double a, b; BRep_Tool::Range(edge, a, b); check(first == a && last == b, "cylinder_io_pcurve_range"); pcWrite(out, pc); }
      out << '\n';
    }
  }
  out << "END\n"; const auto text = out.str(); check(bool(out) && text.size() <= MaxBytes, "cylinder_io_record_size");
  output.write(text.data(), std::streamsize(text.size())); check(bool(output), "cylinder_io_output");
}

void BRepClass3d_ImplicitCylinderIO::Read(Standard_IStream& input, StoredSolid& destination)
{
  check(destination.Shape.IsNull() || !destination.Shape.Locked(), "cylinder_io_locked_destination");
  std::string text; char block[4096];
  while (input.read(block, sizeof(block)) || input.gcount())
  { check(text.size() + size_t(input.gcount()) <= MaxBytes, "cylinder_io_input_size"); text.append(block, size_t(input.gcount())); }
  check(input.eof() && !input.bad(), "cylinder_io_input");
  std::istringstream in(text); in.imbue(std::locale::classic()); word(in, "OCCT_IMPLICIT_CYLINDER_SOLID");
  check(integer(in, 1, 1) == 1, "cylinder_io_version"); const int ps = integer(in, 1, int(MaxMath)), ms = integer(in, 1, int(MaxMath));
  check(in.get() == '\n', "cylinder_io_header"); std::string p(size_t(ps), '\0'), m(size_t(ms), '\0');
  check(bool(in.read(p.data(), ps)) && bool(in.read(m.data(), ms)), "cylinder_io_math_record");
  StoredSolid value; std::istringstream pi(p), mi(m); BRep_ImplicitIO::Read(pi, value.Primitive); TopoDS_Face bare; BRep_ImplicitIO::Read(mi, bare);
  word(in, "META"); auto& s = value.Parameters; s.Frame = frame(in);
  s.Radius = real(in); s.AngleMinimum = real(in); s.AngleMaximum = real(in); s.HeightMinimum = real(in); s.HeightMaximum = real(in);
  s.DistanceTolerance = real(in); s.MinimumCells = integer(in, 4, 512); s.MaximumCells = integer(in, s.MinimumCells, 512);
  check((s.MinimumCells & (s.MinimumCells - 1)) == 0 && (s.MaximumCells & (s.MaximumCells - 1)) == 0, "cylinder_io_resolution");
  value.Negative = integer(in, 0, 1) != 0; value.StoredAxialDirection = integer(in, -1, 1);
  value.StoredBoundaryBound = real(in); value.StoredVolumeLower = real(in); value.StoredVolumeUpper = real(in); value.StoredVolumeTolerance = real(in);
  word(in, "VERTICES"); integer(in, 3, 3); BRep_Builder builder; std::vector<TopoDS_Vertex> vertices(3);
  for (auto& vertex : vertices) { const auto pt = point(in); const double tol = real(in);
    check(tol > 0 && tol <= s.DistanceTolerance, "cylinder_io_vertex_tolerance"); builder.MakeVertex(vertex, pt, tol); }
  word(in, "EDGES"); integer(in, 4, 4); std::vector<TopoDS_Edge> edges(4);
  for (auto& edge : edges)
  {
    word(in, "EDGE"); const double first = real(in), last = real(in), tol = real(in);
    check(last > first && tol > 0 && tol <= s.DistanceTolerance, "cylinder_io_edge_range");
    auto a = vertices[size_t(integer(in, 1, 3) - 1)], b = vertices[size_t(integer(in, 1, 3) - 1)];
    std::string kind; check(bool(in >> kind), "cylinder_io_curve_kind"); occ::handle<Geom_Curve> curve;
    if (kind == "C") { const auto position = frame(in); const double r = real(in); check(r > 0, "cylinder_io_radius"); curve = new Geom_Circle(position, r); }
    else if (kind == "L") { const auto position = point(in); const auto dir = direction(in); curve = new Geom_Line(position, dir); }
    else { check(kind == "B", "cylinder_io_curve_kind"); const int n = integer(in, 2, int(MaxPoles));
      check(first == 0 && last == n - 1, "cylinder_io_spline_range"); NCollection_Array1<gp_Pnt> poles(1,n);
      NCollection_Array1<double> k(1,n); NCollection_Array1<int> m(1,n);
      for (int i=1;i<=n;++i) { poles.SetValue(i,point(in)); k.SetValue(i,i-1); m.SetValue(i,i==1 || i==n ? 2:1); }
      curve = new Geom_BSplineCurve(poles,k,m,1); }
    check(curve->Value(first).Distance(BRep_Tool::Pnt(a)) <= tol && curve->Value(last).Distance(BRep_Tool::Pnt(b)) <= tol, "cylinder_io_edge_endpoints");
    builder.MakeEdge(edge, curve, tol); a.Orientation(TopAbs_FORWARD); b.Orientation(TopAbs_REVERSED);
    builder.Add(edge,a); builder.Add(edge,b); builder.UpdateVertex(a,first,edge,BRep_Tool::Tolerance(a)); builder.UpdateVertex(b,last,edge,BRep_Tool::Tolerance(b)); builder.Range(edge,first,last);
  }
  for (int f=0;f<3;++f)
  {
    word(in,"FACE"); integer(in,f,f); const int orient=integer(in,0,1); const double tol=real(in);
    check(tol>0 && tol<=s.DistanceTolerance,"cylinder_io_face_tolerance"); TopoDS_Face face; occ::handle<Geom_Surface> surface;
    if (f==0) surface=new Geom_Plane(frame(in));
    if (f==1) { const auto position=frame(in); const double r=real(in); check(r>0,"cylinder_io_radius"); surface=new Geom_CylindricalSurface(position,r); }
    if (f<2) builder.MakeFace(face,surface,tol);
    else { face=bare; check(orientation(face)==orient && BRep_Tool::Tolerance(face)==tol,"cylinder_io_bare_face"); }
    const int n=integer(in,f==0?2:f==1?5:1,f==0?2:f==1?5:1); TopoDS_Wire wire; builder.MakeWire(wire);
    struct PC { occ::handle<Geom2d_Curve> Forward, Reverse; int Count=0; }; std::array<PC,4> curves;
    for(int i=0;i<n;++i) { const int id=integer(in,1,4)-1, direction=integer(in,0,1); auto edge=edges[size_t(id)];
      edge.Orientation(direction?TopAbs_REVERSED:TopAbs_FORWARD); builder.Add(wire,edge);
      if(f<2) { auto& entry=curves[size_t(id)]; ++entry.Count; auto pc=pcRead(in);
        if(direction) { check(entry.Reverse.IsNull(),"cylinder_io_duplicate_use"); entry.Reverse=pc; }
        else { check(entry.Forward.IsNull(),"cylinder_io_duplicate_use"); entry.Forward=pc; } }
    }
    if(f<2) for(size_t i=0;i<curves.size();++i) { const auto& pc=curves[i]; if(!pc.Count)continue;
      const auto& edge=edges[i]; const double et=BRep_Tool::Tolerance(edge); double first,last; BRep_Tool::Range(edge,first,last);
      if(pc.Count==1) builder.UpdateEdge(edge,pc.Forward.IsNull()?pc.Reverse:pc.Forward,surface,TopLoc_Location(),et);
      else { check(f==1 && pc.Count==2 && !pc.Forward.IsNull() && !pc.Reverse.IsNull(),"cylinder_io_seam_uses");
        builder.UpdateEdge(edge,pc.Forward,pc.Reverse,surface,TopLoc_Location(),et); }
      builder.Range(edge,surface,TopLoc_Location(),first,last); builder.SameRange(edge,true); builder.SameParameter(edge,true); }
    wire.Closed(true); face.Orientation(TopAbs_FORWARD);
    if(f==2 && orient) wire.Reverse();
    builder.Add(face,wire); face.Orientation(orient?TopAbs_REVERSED:TopAbs_FORWARD);
    value.Faces.push_back(face);
  }
  word(in,"END"); in>>std::ws; check(in.eof(),"cylinder_io_trailing_data");
  TopoDS_Shell shell; builder.MakeShell(shell); for(const auto& face:value.Faces)builder.Add(shell,face);
  shell.Closed(true); builder.MakeSolid(value.Shape); builder.Add(value.Shape,shell); value.Shape.Closed(true);
  validate(value); value.TopologyValidated=true; destination=std::move(value);
}

BRepClass3d_ImplicitCylinderIO::CurveCoherence
BRepClass3d_ImplicitCylinderIO::CheckCurveCoherence(const StoredSolid& value,
                                                 size_t maximumIntervals)
{
  check(maximumIntervals > 0 && maximumIntervals <= 2000000,
        "cylinder_io_invalid_coherence_budget");
  CurveCoherence result;
  try
  {
    validate(value);
    for (size_t f = 0; f < 2; ++f)
    {
      const auto& face = value.Faces[f];
      const auto surface = BRep_Tool::Surface(face);
      const auto cylinder = occ::down_cast<Geom_CylindricalSurface>(surface);
      const double surfaceLipschitz = cylinder.IsNull() ? 1. : std::max(1., cylinder->Radius());
      for (const auto& edge : ordered(effectiveWire(face)))
      {
        TopLoc_Location location;
        double first, last, pcFirst, pcLast;
        const auto curve = BRep_Tool::Curve(edge, location, first, last);
        const auto pc = BRep_Tool::CurveOnSurface(edge, face, pcFirst, pcLast);
        check(!curve.IsNull() && !pc.IsNull() && location.IsIdentity()
                && first == pcFirst && last == pcLast, "cylinder_io_coherence_range");
        const auto spline = occ::down_cast<Geom_BSplineCurve>(curve);
        const auto pcSpline = occ::down_cast<Geom2d_BSplineCurve>(pc);
        const auto circle = occ::down_cast<Geom_Circle>(curve);
        const auto pcCircle = occ::down_cast<Geom2d_Circle>(pc);
        check((!spline.IsNull() && spline->Degree() == 1 && !spline->IsRational())
                || !circle.IsNull() || !occ::down_cast<Geom_Line>(curve).IsNull(),
              "cylinder_io_coherence_curve");
        check((!pcSpline.IsNull() && pcSpline->Degree() == 1 && !pcSpline->IsRational())
                || !pcCircle.IsNull() || !occ::down_cast<Geom2d_Line>(pc).IsNull(),
              "cylinder_io_coherence_pcurve");
        std::vector<double> breaks {first, last};
        if (!spline.IsNull()) for (int i = 1; i <= spline->NbKnots(); ++i)
          if (spline->Knot(i) > first && spline->Knot(i) < last) breaks.push_back(spline->Knot(i));
        if (!pcSpline.IsNull()) for (int i = 1; i <= pcSpline->NbKnots(); ++i)
          if (pcSpline->Knot(i) > first && pcSpline->Knot(i) < last) breaks.push_back(pcSpline->Knot(i));
        std::sort(breaks.begin(), breaks.end());
        breaks.erase(std::unique(breaks.begin(), breaks.end()), breaks.end());
        // Each initial interval lies in one affine spline span. Circle and line
        // speeds are constant. Their sum bounds the difference everywhere.
        for (size_t i = 1; i < breaks.size(); ++i)
        {
          const double a = breaks[i-1], b = breaks[i];
          const double speed = spline.IsNull() ? (circle.IsNull() ? 1. : circle->Radius())
            : curve->Value(a).Distance(curve->Value(b)) / (b-a);
          const double pcSpeed = pcSpline.IsNull() ? (pcCircle.IsNull() ? 1. : pcCircle->Radius())
            : pc->Value(a).Distance(pc->Value(b)) / (b-a);
          const double lipschitz = speed + surfaceLipschitz * pcSpeed;
          const double acceleration = (circle.IsNull() ? 0. : circle->Radius())
            + (cylinder.IsNull() ? 0. : cylinder->Radius() * pcSpeed * pcSpeed)
            + surfaceLipschitz * (pcCircle.IsNull() ? 0. : pcCircle->Radius());
          check(std::isfinite(lipschitz), "cylinder_io_coherence_speed");
          std::vector<std::pair<double,double>> pending {{a,b}};
          while (!pending.empty())
          {
            check(result.Intervals < maximumIntervals, "cylinder_io_coherence_budget");
            ++result.Intervals;
            const auto interval = pending.back(); pending.pop_back();
            const double middle = interval.first + (interval.second-interval.first)*.5;
            const auto p = curve->Value(middle);
            const auto uv = pc->Value(middle);
            const auto q = surface->Value(uv.X(),uv.Y());
            const double reserve = 32768 * std::numeric_limits<double>::epsilon()
              * (1 + p.Distance(gp_Pnt()) + q.Distance(gp_Pnt()) + lipschitz
                 * (1 + std::abs(interval.first) + std::abs(interval.second)));
            const double distance = p.Distance(q);
            const auto residual = [&](double parameter) {
              const auto uvEnd = pc->Value(parameter);
              return curve->Value(parameter).Distance(surface->Value(uvEnd.X(),uvEnd.Y()));
            };
            const double width = interval.second-interval.first;
            // Linear interpolation of the residual vector has a whole-interval
            // remainder <= sup|residual''| * width^2 / 8. No sample-only claim.
            const double interpolation = std::max(residual(interval.first), residual(interval.second))
              + acceleration * width * width / 8;
            const double bound = std::min(distance + lipschitz*width*.5, interpolation) + reserve;
            const double tolerance = BRep_Tool::Tolerance(edge);
            check(std::isfinite(bound) && distance <= tolerance,
                  "cylinder_io_coherence_mismatch");
            if (bound <= tolerance) { result.Bound = std::max(result.Bound,bound); continue; }
            check(middle > interval.first && middle < interval.second,
                  "cylinder_io_coherence_precision");
            pending.emplace_back(middle,interval.second);
            pending.emplace_back(interval.first,middle);
          }
        }
      }
    }
    result.Verified = true; result.Reason = "cylinder_stored_curve_coherence_verified";
  }
  catch (const Standard_Failure& error)
  { result.Reason = error.what() ? error.what() : "cylinder_io_coherence_unresolved"; }
  return result;
}

namespace
{
IO::CurveCoherence checkStoredSeam(const IO::StoredSolid& value,
                                 size_t maximumIntervals, bool masked)
{
  check(maximumIntervals > 0 && maximumIntervals <= 2000000,
        "cylinder_io_invalid_seam_budget");
  IO::CurveCoherence result;
  try
  {
    validate(value);
    const auto data = occ::down_cast<BRep_TImplicitFace>(value.Primitive.TShape());
    const auto* field = static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
    const double scale = std::abs((value.Primitive.Location() * data->Location()).Transformation().ScaleFactor());
    const double k = 2 * std::acos(-1.) / field->Period() / scale;
    const double hessian = 4 * k * k;
    check(std::isfinite(hessian) && hessian > 0, "cylinder_io_seam_hessian");
    for (const auto& edge : ordered(effectiveWire(value.Faces[2])))
    {
      TopLoc_Location location;
      double first,last;
      const auto curve = occ::down_cast<Geom_BSplineCurve>(BRep_Tool::Curve(edge,location,first,last));
      check(!curve.IsNull() && location.IsIdentity() && curve->Degree() == 1
        && !curve->IsRational() && !curve->IsPeriodic()
        && first == curve->FirstParameter() && last == curve->LastParameter(),
        "cylinder_io_seam_curve");
      const double tolerance = std::min(value.Parameters.DistanceTolerance, BRep_Tool::Tolerance(value.Faces[2]));
      for (int i = 2; i <= curve->NbPoles(); ++i)
      {
        std::vector<std::pair<gp_Pnt,gp_Pnt>> pending {{curve->Pole(i-1),curve->Pole(i)}};
        while (!pending.empty())
        {
          check(result.Intervals < maximumIntervals, "cylinder_io_seam_budget");
          ++result.Intervals;
          const auto segment = pending.back(); pending.pop_back();
          const gp_Pnt originalMiddle((segment.first.XYZ()+segment.second.XYZ())*.5);
          const double originalRadius = segment.first.Distance(segment.second)*.5;
          const double reserve = 32768 * std::numeric_limits<double>::epsilon()
            * (1 + originalMiddle.Distance(gp_Pnt()) + originalRadius + tolerance)
            * (1 + k + hessian + std::abs(field->IsoValue()));
          auto start = segment.first, finish = segment.second;
          double displacement = 0;
          const gp_Vec axis(value.Parameters.Frame.Direction());
          if (masked)
          {
            const auto& frame = value.Parameters.Frame;
            const auto radial = [&](const gp_Pnt& p) {
              const gp_Vec delta(frame.Location(),p);
              return delta-axis*delta.Dot(axis);
            };
            const double maximumRadius = std::max(radial(start).Magnitude(),radial(finish).Magnitude());
            check(value.Parameters.Radius > reserve, "cylinder_io_masked_seam_precision");
            const double factor = std::min(1., (value.Parameters.Radius-reserve)/(maximumRadius+reserve));
            displacement = (1-factor)*maximumRadius + reserve;
            start.Translate(radial(start)*(factor-1));
            finish.Translate(radial(finish)*(factor-1));
          }
          const gp_Pnt middle((start.XYZ()+finish.XYZ())*.5);
          const double radius = start.Distance(finish)*.5;
          const double reach = tolerance-displacement;
          check(reach > reserve, "cylinder_io_masked_seam_displacement");
          if (masked)
          {
            // The source box and finite cylinder are convex. Four endpoints
            // enclose every axial witness segment for the whole affine span.
            for (const auto& p : {start,finish}) for (double sign : {-1.,1.})
              check(BRep_ImplicitTool::IsInDomain(value.Faces[2],p.Translated(axis*(sign*reach))),
                    "cylinder_io_masked_seam_domain");
          }
          const double a = BRep_ImplicitTool::Value(value.Primitive,start);
          const double b = BRep_ImplicitTool::Value(value.Primitive,finish);
          const auto gradient = BRep_ImplicitTool::Gradient(value.Primitive,middle);
          // The fixed direction is the cylinder axis for masked witnesses,
          // otherwise the midpoint gradient. Its signed derivative is >= lower
          // in every segment point's +/- reach neighbourhood.
          // Taylor's endpoint interpolation bounds F over the complete segment.
          const double derivative = masked ? std::abs(gradient.Dot(axis)) : gradient.Magnitude();
          const double lower = derivative - hessian*(radius+reach) - reserve;
          const double residual = std::max(std::abs(a),std::abs(b)) + hessian*radius*radius*.5 + reserve;
          const double bound = lower > 0 ? displacement+residual/lower : std::numeric_limits<double>::infinity();
          check(std::isfinite(a) && std::isfinite(b) && std::isfinite(residual)
            && std::isfinite(gradient.Magnitude()), "cylinder_io_seam_nonfinite");
          check(std::max(std::abs(a),std::abs(b)) <= field->LipschitzBound()/scale*tolerance + reserve,
                "cylinder_io_seam_outside_tolerance");
          if (bound <= tolerance)
          { result.Bound = std::max(result.Bound,bound); continue; }
          check(originalRadius > reserve && originalMiddle.Distance(segment.first) > 0
            && originalMiddle.Distance(segment.second) > 0, "cylinder_io_seam_unresolved");
          pending.emplace_back(originalMiddle,segment.second);
          pending.emplace_back(segment.first,originalMiddle);
        }
      }
    }
    result.Verified = true;
    result.Reason = masked ? "cylinder_stored_masked_seam_verified" : "cylinder_stored_primitive_seam_verified";
  }
  catch (const Standard_Failure& error)
  { result.Reason = error.what() ? error.what() : "cylinder_io_seam_unresolved"; }
  return result;
}
}

IO::CurveCoherence IO::CheckPrimitiveSeam(const StoredSolid& value, size_t maximumIntervals)
{ return checkStoredSeam(value,maximumIntervals,false); }

IO::CurveCoherence IO::CheckMaskedSeam(const StoredSolid& value, size_t maximumIntervals)
{ return checkStoredSeam(value,maximumIntervals,true); }

IO::CurveCoherence IO::CheckSeamCoverage(const StoredSolid& value, size_t maximumIntervals)
{
  check(maximumIntervals > 0 && maximumIntervals <= 2000000,
        "cylinder_io_invalid_coverage_budget");
  CurveCoherence result;
  try
  {
    validate(value);
    const auto& s = value.Parameters;
    const auto axial = BRepClass3d_ImplicitCylinder::CheckAxialField(value.Primitive,s,24,maximumIntervals);
    check(axial.MonotonicityVerified, axial.Reason.c_str());
    result.Intervals = axial.Cells;
    const auto data = occ::down_cast<BRep_TImplicitFace>(value.Primitive.TShape());
    const auto* field = static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
    const double scale = std::abs((value.Primitive.Location()*data->Location()).Transformation().ScaleFactor());
    const double k = 2*std::acos(-1.)/field->Period()/scale;
    const double hessian = 4*k*k, lipschitz = field->LipschitzBound()/scale;
    const gp_Vec x(s.Frame.XDirection()), y(s.Frame.YDirection()), axis(s.Frame.Direction());
    const double reserve = 32768*std::numeric_limits<double>::epsilon()
      * (1+s.Frame.Location().Distance(gp_Pnt())+s.Radius+std::abs(s.HeightMinimum)+std::abs(s.HeightMaximum))
      * (1+hessian+lipschitz+std::abs(field->IsoValue()));
    double winding = 0;
    int sense = 0;
    const auto uses = ordered(effectiveWire(value.Faces[2]));
    check(uses.size() == 1, "cylinder_io_coverage_wire");
    double first,last;
    const auto curve = occ::down_cast<Geom_BSplineCurve>(BRep_Tool::Curve(uses[0],first,last));
    check(!curve.IsNull() && curve->Degree() == 1 && !curve->IsRational() && !curve->IsPeriodic()
      && first == curve->FirstParameter() && last == curve->LastParameter()
      && curve->Pole(1).IsEqual(curve->Pole(curve->NbPoles()),0.), "cylinder_io_coverage_closed_curve");
    for (int i = 2; i <= curve->NbPoles(); ++i)
    {
      check(result.Intervals < maximumIntervals, "cylinder_io_coverage_budget"); ++result.Intervals;
      const gp_Vec a(s.Frame.Location(),curve->Pole(i-1)), b(s.Frame.Location(),curve->Pole(i));
      const double ax=a.Dot(x), ay=a.Dot(y), bx=b.Dot(x), by=b.Dot(y);
      const double ra=std::hypot(ax,ay), rb=std::hypot(bx,by);
      check(ra > reserve && rb > reserve, "cylinder_io_coverage_axis");
      const double angle=std::atan2(ax*by-ay*bx,ax*bx+ay*by);
      check(std::abs(angle)*s.Radius > reserve && std::abs(angle) < std::acos(-1.)*.5,
            "cylinder_io_coverage_angular_span");
      const int direction=angle > 0 ? 1 : -1;
      if (!sense) sense=direction;
      check(sense == direction,"cylinder_io_coverage_angular_reversal");
      winding+=angle;
      const double za=a.Dot(axis), zb=b.Dot(axis), theta=std::atan2(ay,ax);
      const auto point = [&](double fraction) {
        const double t=theta+angle*fraction, z=za+(zb-za)*fraction;
        return s.Frame.Location().Translated(x*(s.Radius*std::cos(t))
          +y*(s.Radius*std::sin(t))+axis*z);
      };
      const double deviation=std::max(std::abs(ra-s.Radius),std::abs(rb-s.Radius))
        +s.Radius*angle*angle/8+reserve;
      const double acceleration=hessian*(s.Radius*s.Radius*angle*angle+(zb-za)*(zb-za))
        +lipschitz*s.Radius*angle*angle;
      const double residual=std::max(std::abs(BRep_ImplicitTool::Value(value.Primitive,point(0))),
        std::abs(BRep_ImplicitTool::Value(value.Primitive,point(1))))+acceleration/8+reserve;
      const double halfLength=.5*std::hypot(s.Radius*angle,zb-za);
      const double localLower=axial.Direction*BRep_ImplicitTool::Gradient(value.Primitive,point(.5)).Dot(axis)
        -hessian*(halfLength+s.DistanceTolerance)-reserve;
      const double root=residual/std::max(axial.DerivativeLowerBound,localLower);
      check(std::isfinite(root) && std::min(za,zb)-root >= s.HeightMinimum+reserve
        && std::max(za,zb)+root <= s.HeightMaximum-reserve, "cylinder_io_coverage_height_bracket");
      const double bound=deviation+root+reserve;
      check(bound <= s.DistanceTolerance,"cylinder_io_coverage_distance");
      result.Bound=std::max(result.Bound,bound);
    }
    check(std::abs(std::abs(winding)-2*std::acos(-1.))*s.Radius <= reserve*curve->NbPoles(),
          "cylinder_io_coverage_winding");
    result.Verified=true; result.Reason="cylinder_stored_seam_reverse_coverage_verified";
  }
  catch (const Standard_Failure& error)
  { result.Reason=error.what() ? error.what() : "cylinder_io_coverage_unresolved"; }
  return result;
}

IO::CurveCoherence IO::CheckRetainedCap(const StoredSolid& value, size_t maximumIntervals)
{
  check(maximumIntervals > 0 && maximumIntervals <= 2000000,
        "cylinder_io_invalid_cap_budget");
  CurveCoherence result;
  try
  {
    validate(value);
    const auto& s = value.Parameters;
    const auto axial = BRepClass3d_ImplicitCylinder::CheckAxialField(value.Primitive,s,24,maximumIntervals);
    check(axial.MonotonicityVerified,axial.Reason.c_str());
    result.Intervals=axial.Cells;
    const auto data=occ::down_cast<BRep_TImplicitFace>(value.Primitive.TShape());
    const auto* field=static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
    const double scale=std::abs((value.Primitive.Location()*data->Location()).Transformation().ScaleFactor());
    const double lipschitz=field->LipschitzBound()/scale;
    const double precision=32768*std::numeric_limits<double>::epsilon()
      *(1+s.Frame.Location().Distance(gp_Pnt())+s.Radius+std::abs(s.HeightMinimum)+std::abs(s.HeightMaximum));
    const double reserve=precision*(1+lipschitz+std::abs(field->IsoValue()));
    check(std::isfinite(reserve) && precision < s.DistanceTolerance,"cylinder_io_cap_precision");
    const gp_Vec x(s.Frame.XDirection()),y(s.Frame.YDirection()),axis(s.Frame.Direction());
    struct Cell { double X0,X1,Y0,Y1; int Depth; };
    std::vector<Cell> pending {{-s.Radius,s.Radius,-s.Radius,s.Radius,0}};
    while (!pending.empty())
    {
      check(result.Intervals < maximumIntervals,"cylinder_io_cap_budget"); ++result.Intervals;
      const auto cell=pending.back(); pending.pop_back();
      const auto nearest=[](double a,double b) { return a > 0 ? a : b < 0 ? -b : 0.; };
      if (std::hypot(nearest(cell.X0,cell.X1),nearest(cell.Y0,cell.Y1)) > s.Radius+precision) continue;
      const double dx=cell.X1-cell.X0,dy=cell.Y1-cell.Y0;
      const auto p=s.Frame.Location().Translated(x*(cell.X0+dx*.5)+y*(cell.Y0+dy*.5));
      const double low=axial.Direction*BRep_ImplicitTool::Value(value.Primitive,p.Translated(axis*s.HeightMinimum));
      const double high=axial.Direction*BRep_ImplicitTool::Value(value.Primitive,p.Translated(axis*s.HeightMaximum));
      const double uncertainty=lipschitz*(std::hypot(dx,dy)*.5+precision)+reserve;
      check(std::isfinite(low) && std::isfinite(high) && std::isfinite(uncertainty),"cylinder_io_cap_nonfinite");
      if (low < -uncertainty && high > uncertainty) continue;
      check(low-uncertainty < 0 && high+uncertainty > 0,"cylinder_io_cap_endpoint_signs");
      check(cell.Depth < 24,"cylinder_io_cap_depth");
      auto a=cell,b=cell; a.Depth=b.Depth=cell.Depth+1;
      if (dx >= dy) { const double cut=cell.X0+dx*.5;
        check(cut > cell.X0 && cut < cell.X1,"cylinder_io_cap_split"); a.X1=b.X0=cut; }
      else { const double cut=cell.Y0+dy*.5;
        check(cut > cell.Y0 && cut < cell.Y1,"cylinder_io_cap_split"); a.Y1=b.Y0=cut; }
      pending.push_back(b); pending.push_back(a);
    }
    const bool bottom=value.Negative == (axial.Direction > 0);
    const double height=bottom ? s.HeightMinimum : s.HeightMaximum;
    const auto center=s.Frame.Location().Translated(axis*height);
    const auto& cap=value.Faces[0];
    const auto plane=occ::down_cast<Geom_Plane>(BRep_Tool::Surface(cap));
    const double outward=bottom ? -1. : 1.;
    check(!plane.IsNull() && plane->Location().Distance(center) <= precision
      && gp_Vec(plane->Axis().Direction()).Crossed(axis).Magnitude()*s.Radius <= precision
      && plane->Axis().Direction().Dot(s.Frame.Direction())*(orientation(cap) ? -1. : 1.)*outward > 1-1.e-12,
      "cylinder_io_cap_surface_orientation");
    const auto uses=ordered(effectiveWire(cap));
    check(uses.size() == 2,"cylinder_io_cap_arc_count");
    double winding=0;
    for (const auto& edge:uses)
    {
      double first,last; TopLoc_Location location;
      const auto circle=occ::down_cast<Geom_Circle>(BRep_Tool::Curve(edge,location,first,last));
      check(!circle.IsNull() && location.IsIdentity()
        && circle->Location().Distance(center) <= precision
        && std::abs(circle->Radius()-s.Radius) <= precision
        && gp_Vec(circle->Axis().Direction()).Crossed(axis).Magnitude()*s.Radius <= precision
        && std::abs(circle->Axis().Direction().Dot(s.Frame.Direction())) > 1-1.e-12
        && last > first && last-first <= 2*std::acos(-1.)+precision/s.Radius,
        "cylinder_io_cap_circle");
      const double direction=circle->Axis().Direction().Dot(s.Frame.Direction()) > 0 ? 1. : -1.;
      const double span=(last-first)*direction*(orientation(edge) ? -1. : 1.);
      check(span*outward > 0,"cylinder_io_cap_arc_orientation");
      winding+=span;
      TopoDS_Vertex a,b; TopExp::Vertices(TopoDS::Edge(edge.Oriented(TopAbs_FORWARD)),a,b);
      check(circle->Value(first).Distance(BRep_Tool::Pnt(a)) <= precision
        && circle->Value(last).Distance(BRep_Tool::Pnt(b)) <= precision,"cylinder_io_cap_arc_endpoints");
      double pcFirst,pcLast;
      const auto pc=occ::down_cast<Geom2d_Circle>(BRep_Tool::CurveOnSurface(edge,cap,pcFirst,pcLast));
      check(!pc.IsNull() && pcFirst == first && pcLast == last,"cylinder_io_cap_pcurve");
      const auto mapped=[&](const gp_Dir2d& d) {
        return gp_Vec(plane->Position().XDirection())*d.X()+gp_Vec(plane->Position().YDirection())*d.Y();
      };
      const double coefficients=plane->Value(pc->Location().X(),pc->Location().Y()).Distance(circle->Location())
        +(mapped(pc->XAxis().Direction())*pc->Radius()-gp_Vec(circle->Position().XDirection())*circle->Radius()).Magnitude()
        +(mapped(pc->YAxis().Direction())*pc->Radius()-gp_Vec(circle->Position().YDirection())*circle->Radius()).Magnitude();
      check(coefficients <= precision,"cylinder_io_cap_pcurve_coefficients");
    }
    check(std::abs(winding*outward-2*std::acos(-1.))*s.Radius <= 8*precision,
      "cylinder_io_cap_circle_coverage");
    check(BRepClass_FaceClassifier(cap,center,precision).State() == TopAbs_IN,
      "cylinder_io_cap_complement");
    result.Bound=8*precision;
    check(result.Bound <= s.DistanceTolerance && result.Bound <= BRep_Tool::Tolerance(cap),
      "cylinder_io_cap_tolerance");
    result.Verified=true; result.Reason="cylinder_stored_retained_cap_verified";
  }
  catch (const Standard_Failure& error)
  { result.Reason=error.what() ? error.what() : "cylinder_io_cap_unresolved"; }
  return result;
}
