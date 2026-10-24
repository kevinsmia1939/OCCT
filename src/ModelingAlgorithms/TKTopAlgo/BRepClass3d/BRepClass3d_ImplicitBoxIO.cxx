// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitBoxIO.hxx>
#include <BRep_ImplicitIO.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <BRepTools_ShapeSet.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_BSplineCurve.hxx>
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
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <typeinfo>

namespace
{
constexpr size_t MaxBytes = 16 * 1024 * 1024;
using Boundary            = BRepClass3d_ImplicitBox::Boundary;
using Map                 = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

void requireRecord(bool condition, const char* invariant)
{
  if (!condition)
  {
    const std::string message = std::string("Invalid native mixed-box topology record or stream: ")
                                + invariant;
    throw Standard_DomainError(message.c_str());
  }
}

// Keep failed graph/stream invariants diagnosable without weakening validation.
#define require(condition) requireRecord((condition), #condition)

// TopTools_ShapeSet selects 15 significant digits for per-shape geometry.
// Exact corner/plane identities need round-trip binary64 vertex coordinates.
// Keep this policy local to native boundary records, without changing general
// B-Rep persistence or introducing a new payload format.
class ExactShapeSet : public BRepTools_ShapeSet
{
public:
  ExactShapeSet() : BRepTools_ShapeSet(false, false) {}

  void WriteGeometry(const TopoDS_Shape& shape, Standard_OStream& output) const override
  {
    const auto precision = output.precision(std::numeric_limits<double>::max_digits10);
    BRepTools_ShapeSet::WriteGeometry(shape, output);
    output.precision(precision);
  }
};

const BRep_TImplicitFace& primitive(const TopoDS_Face& face)
{
  const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  if (data.IsNull() || typeid(*data) != typeid(BRep_TImplicitFace)
      || typeid(*data->Field()) != typeid(Geom_ImplicitGyroid) || !face.Location().IsIdentity()
      || !data->Location().IsIdentity()
      || (face.Orientation() != TopAbs_FORWARD && face.Orientation() != TopAbs_REVERSED))
  {
    std::ostringstream message;
    message << "Topology I/O requires an identity-placed native gyroid/box boundary: native="
            << !data.IsNull() << ", instance_identity=" << face.Location().IsIdentity()
            << ", orientation=" << face.Orientation();
    if (!data.IsNull())
      message << ", geometry_identity=" << data->Location().IsIdentity()
              << ", face_type=" << typeid(*data).name()
              << ", field_type=" << typeid(*data->Field()).name();
    throw Standard_NotImplemented(message.str().c_str());
  }
  return *data;
}

double distance(const gp_Pnt& a, const gp_Pnt& b)
{
  return std::hypot(a.X() - b.X(), a.Y() - b.Y(), a.Z() - b.Z());
}

void resetTransientFlags(const TopoDS_Shape& shape)
{
  Map parts;
  TopExp::MapShapes(shape, parts);
  for (auto part : parts)
  {
    part.Locked(false);
    part.Checked(false);
    part.Modified(true);
    part.Convex(false);
  }
}

// Validate the stored graph, NOT the original numerical tracing certificate.
// The native implicit face must never be passed to UV-based BRepCheck.
void validate(const Boundary& boundary)
{
  require(!boundary.Shape.IsNull() && boundary.Shape.Location().IsIdentity()
          && boundary.Shape.Orientation() == TopAbs_FORWARD && !boundary.Faces.empty()
          && boundary.Faces.size() <= 128 && boundary.Faces.size() == boundary.Sources.size());
  const auto& implicit = boundary.Faces.back();
  const auto& data     = primitive(implicit);
  require(boundary.Sources.back() == -1 && std::isfinite(boundary.GeometricBound)
          && boundary.GeometricBound >= 0 && boundary.GeometricBound <= data.Tolerance());
  size_t faceIndex = 0;
  for (TopoDS_Iterator it(boundary.Shape); it.More(); it.Next(), ++faceIndex)
    require(faceIndex < boundary.Faces.size() && it.Value().IsEqual(boundary.Faces[faceIndex]));
  require(faceIndex == boundary.Faces.size());
  double magnitude = 1;
  for (int axis = 1; axis <= 3; ++axis)
    magnitude = std::max(
      {magnitude, std::abs(data.Minimum().Coord(axis)), std::abs(data.Maximum().Coord(axis))});
  const double reserve   = 1024 * std::numeric_limits<double>::epsilon() * magnitude;
  const double precision = std::max(reserve, std::min(1.e-7, data.Tolerance() / 4096));
  Map          allEdges, nativeEdges, nativeVertices;
  TopExp::MapShapes(boundary.Shape, TopAbs_EDGE, allEdges);
  TopExp::MapShapes(implicit, TopAbs_EDGE, nativeEdges);
  TopExp::MapShapes(implicit, TopAbs_VERTEX, nativeVertices);
  require(allEdges.Extent() <= 500000 && nativeEdges.Extent() == boundary.SeamCount);
  const bool regularized = boundary.Policy == BRepClass3d_ImplicitBox::CornerPolicy::Regularized;
  require(regularized || boundary.Policy == BRepClass3d_ImplicitBox::CornerPolicy::Reject);
  BRepClass3d_ImplicitBox::CornerModel corners;
  std::array<int, 8> cornerVertices{}, cornerEnds{};
  if (regularized)
  {
    // Recheck the local mathematical identities, never retrace or recertify
    // stored seam bounds. Geometry and topology remain the stored graph.
    corners = BRepClass3d_ImplicitBox::CheckCornerGerms(implicit, data.Tolerance());
    require(corners.Resolved && boundary.ExactCornerCount == corners.ExactZeroCount);
    const auto exclusions = BRepClass3d_ImplicitBox::SurfaceCornerExclusions(implicit);
    require(exclusions.size() == data.ExcludedDomains().size());
    for (size_t i = 0; i < exclusions.size(); ++i)
      for (int end = 0; end < 2; ++end)
        require(distance(exclusions[i][size_t(end)], data.ExcludedDomains()[i][size_t(end)]) == 0);
    std::vector<gp_Pnt> excluded;
    int contacts = 0;
    for (const auto& corner : corners.Corners)
      if (corner.ExactZero)
      {
        using Kind = BRepClass3d_ImplicitBox::CornerKind;
        if (corner.Kind == (implicit.Orientation() == TopAbs_FORWARD
                             ? Kind::PositiveCone : Kind::NegativeCone))
          excluded.push_back(corner.Anchor);
        for (int axis = 0; axis < 3; ++axis)
          contacts += (corner.InwardDerivativeLower[size_t((axis + 1) % 3)] > 0)
                      == (corner.InwardDerivativeLower[size_t((axis + 2) % 3)] > 0);
      }
    require(boundary.PlanarCornerContacts == contacts
            && boundary.ExcludedCornerContacts.size() == excluded.size());
    for (size_t i = 0; i < excluded.size(); ++i)
      require(distance(boundary.ExcludedCornerContacts[i], excluded[i]) == 0);
    Map vertices;
    TopExp::MapShapes(boundary.Shape, TopAbs_VERTEX, vertices);
    for (const auto& vertex : vertices)
      for (size_t i = 0; i < 8; ++i)
        if (distance(BRep_Tool::Pnt(TopoDS::Vertex(vertex)), corners.Corners[i].Anchor) == 0)
          ++cornerVertices[i];
    for (size_t i = 0; i < 8; ++i)
    {
      const auto membership = BRepClass3d_ImplicitBox::ClassifyCornerGerm(
        implicit, int(i), corners.Corners[i].Anchor,
        implicit.Orientation() == TopAbs_FORWARD ? BRepClass3d_ImplicitBox::Operation::Common
                                                 : BRepClass3d_ImplicitBox::Operation::BoxMinusField,
        data.Tolerance());
      require(membership.Resolved && cornerVertices[i] == int(membership.RegularizedRetained));
    }
    for (const auto& edge : nativeEdges)
    {
      TopoDS_Vertex a, b;
      TopExp::Vertices(TopoDS::Edge(edge), a, b);
      for (const auto& vertex : {a, b})
        for (size_t i = 0; i < 8; ++i)
          if (distance(BRep_Tool::Pnt(vertex), corners.Corners[i].Anchor) == 0)
            ++cornerEnds[i];
    }
    int ends = 0;
    for (size_t i = 0; i < 8; ++i)
    {
      require(cornerEnds[i] == (corners.Corners[i].Kind
                                == BRepClass3d_ImplicitBox::CornerKind::SplitCone ? 2 : 0));
      ends += cornerEnds[i];
    }
    require(boundary.CornerContourEnds == ends);
  }
  else
    require(boundary.ExactCornerCount == 0 && boundary.CornerContourEnds == 0
            && boundary.PlanarCornerContacts == 0 && boundary.ExcludedCornerContacts.empty()
            && data.ExcludedDomains().empty());
  int rootCount = 0;
  for (const auto& vertex : nativeVertices)
  {
    const auto p     = BRep_Tool::Pnt(TopoDS::Vertex(vertex));
    int        sides = 0;
    for (int axis = 1; axis <= 3; ++axis)
      sides += std::min(std::abs(p.Coord(axis) - data.Minimum().Coord(axis)),
                        std::abs(p.Coord(axis) - data.Maximum().Coord(axis)))
               <= precision;
    require(sides >= 1 && sides <= (regularized ? 3 : 2));
    if (sides == 3)
    {
      bool matched = false;
      for (const auto& corner : corners.Corners)
        matched |= corner.ExactZero
                   && corner.Kind == BRepClass3d_ImplicitBox::CornerKind::SplitCone
                   && distance(p, corner.Anchor) == 0;
      require(matched);
    }
    rootCount += sides == 2;
  }
  require(rootCount == boundary.BoxEdgeRootCount);
  for (const auto& edgeShape : allEdges)
  {
    const auto      edge = TopoDS::Edge(edgeShape);
    TopLoc_Location location;
    double          first, last;
    const auto      curve  = BRep_Tool::Curve(edge, location, first, last);
    const auto      spline = occ::down_cast<Geom_BSplineCurve>(curve);
    require(location.IsIdentity() && !spline.IsNull()
            && typeid(*spline) == typeid(Geom_BSplineCurve) && spline->Degree() == 1
            && !spline->IsRational() && !spline->IsPeriodic() && spline->NbPoles() >= 2
            && spline->NbPoles() <= 300000 && spline->NbKnots() == spline->NbPoles() && first == 0
            && last == spline->NbPoles() - 1);
    for (int i = 1; i <= spline->NbPoles(); ++i)
    {
      require(spline->Knot(i) == i - 1
              && spline->Multiplicity(i) == ((i == 1 || i == spline->NbPoles()) ? 2 : 1));
      const auto p = spline->Pole(i);
      for (int axis = 1; axis <= 3; ++axis)
        require(std::isfinite(p.Coord(axis))
                && p.Coord(axis) >= data.Minimum().Coord(axis) - precision
                && p.Coord(axis) <= data.Maximum().Coord(axis) + precision);
    }
    TopoDS_Vertex a, b;
    TopExp::Vertices(edge, a, b);
    require(!a.IsNull() && !b.IsNull() && distance(BRep_Tool::Pnt(a), spline->Pole(1)) <= precision
            && distance(BRep_Tool::Pnt(b), spline->Pole(spline->NbPoles())) <= precision);
  }
  std::vector<int> uses(size_t(allEdges.Extent()), 0), balance(uses.size(), 0),
    capUses(uses.size(), 0);
  for (size_t i = 0; i < boundary.Faces.size(); ++i)
  {
    const auto& face = boundary.Faces[i];
    if (i + 1 < boundary.Faces.size())
    {
      require(boundary.Sources[i] >= 0 && boundary.Sources[i] < 6 && face.Location().IsIdentity()
              && face.Orientation() == TopAbs_FORWARD && !BRep_Tool::IsImplicit(face)
              && BRep_Tool::HasParametricSurface(face) && BRepCheck_Analyzer(face).IsValid());
      const BRepAdaptor_Surface adaptor(face);
      require(adaptor.GetType() == GeomAbs_Plane);
      const auto plane  = adaptor.Plane();
      const auto normal = plane.Axis().Direction();
      int        axis   = 0;
      for (int j = 1; j <= 3; ++j)
        if (std::abs(normal.Coord(j)) > 1 - 1.e-12)
          axis = j;
      require(axis != 0);
      const double coordinate =
        normal.Coord(axis) > 0 ? data.Maximum().Coord(axis) : data.Minimum().Coord(axis);
      require(std::abs(plane.Location().Coord(axis) - coordinate) <= precision);
      if (regularized)
      {
        require(plane.Location().Coord(axis) == coordinate);
        for (int j = 1; j <= 3; ++j)
          require(normal.Coord(j) == (j == axis ? (normal.Coord(axis) > 0 ? 1 : -1) : 0));
      }
    }
    for (TopoDS_Iterator wire(face); wire.More(); wire.Next())
    {
      require(wire.Value().ShapeType() == TopAbs_WIRE && wire.Value().Location().IsIdentity());
      TopoDS_Vertex start, previous;
      bool          nonempty = false;
      for (TopoDS_Iterator edge(wire.Value()); edge.More(); edge.Next())
      {
        require(edge.Value().ShapeType() == TopAbs_EDGE);
        TopoDS_Vertex a, b;
        TopExp::Vertices(TopoDS::Edge(edge.Value()), a, b, true);
        require(!a.IsNull() && !b.IsNull());
        if (!nonempty)
          start = a;
        else
          require(previous.IsSame(a));
        previous = b;
        nonempty = true;
      }
      require(nonempty && start.IsSame(previous));
    }
    for (TopExp_Explorer edge(face, TopAbs_EDGE); edge.More(); edge.Next())
    {
      const int index = allEdges.FindIndex(edge.Current()) - 1;
      require(index >= 0
              && (edge.Current().Orientation() == TopAbs_FORWARD
                  || edge.Current().Orientation() == TopAbs_REVERSED));
      ++uses[size_t(index)];
      balance[size_t(index)] += edge.Current().Orientation() == TopAbs_FORWARD ? 1 : -1;
      capUses[size_t(index)] += i + 1 < boundary.Faces.size();
    }
  }
  for (size_t i = 0; i < uses.size(); ++i)
    require(uses[i] == 2 && balance[i] == 0);
  for (const auto& edge : nativeEdges)
    require(capUses[size_t(allEdges.FindIndex(edge) - 1)] == 1);
}
} // namespace

void BRepClass3d_ImplicitBoxIO::Write(const Boundary& boundary, Standard_OStream& output)
{
  validate(boundary);
  const auto&  face = boundary.Faces.back();
  const auto&  data = primitive(face);
  BRep_Builder builder;
  TopoDS_Face  bare;
  builder.MakeImplicitFace(bare,
                           data.Field(),
                           data.Minimum(),
                           data.Maximum(),
                           TopLoc_Location(),
                           data.Tolerance());
  bare.Orientation(face.Orientation());
  bare.TShape(data.EmptyCopy());
  std::ostringstream math;
  BRep_ImplicitIO::Write(bare, math);
  TopoDS_Compound payload, wires;
  builder.MakeCompound(payload);
  builder.MakeCompound(wires);
  for (size_t i = 0; i + 1 < boundary.Faces.size(); ++i)
    builder.Add(payload, boundary.Faces[i]);
  for (TopoDS_Iterator it(face); it.More(); it.Next())
    builder.Add(wires, it.Value());
  builder.Add(payload, wires);
  // Strip generated meshes from a deep copy; never mutate the source caches.
  const auto clean = BRepBuilderAPI_Copy(payload, true, false).Shape();
  resetTransientFlags(clean);
  BRepTools::Clean(clean, true);
  std::ostringstream topology;
  topology.imbue(std::locale::classic());
  ExactShapeSet shapeSet;
  shapeSet.SetFormatNb(TopTools_FormatVersion_VERSION_3);
  shapeSet.Add(clean);
  shapeSet.Write(topology);
  shapeSet.Write(clean, topology);
  require(bool(topology));
  const auto         m = math.str(), t = topology.str();
  std::ostringstream record;
  record.imbue(std::locale::classic());
  record << std::setprecision(std::numeric_limits<double>::max_digits10)
         << "OCCT_IMPLICIT_BOX_TOPOLOGY "
         << (boundary.Policy == BRepClass3d_ImplicitBox::CornerPolicy::Regularized ? 2 : 1)
         << ' ' << m.size() << ' ' << t.size() << ' '
         << boundary.Faces.size() << ' ' << boundary.SeamCount << ' ' << boundary.BoxEdgeRootCount
         << ' ' << boundary.GeometricBound;
  if (boundary.Policy == BRepClass3d_ImplicitBox::CornerPolicy::Regularized)
    record << " 1 " << boundary.ExactCornerCount << ' ' << boundary.CornerContourEnds
           << ' ' << boundary.PlanarCornerContacts << ' ' << boundary.ExcludedCornerContacts.size();
  record << '\n';
  for (int source : boundary.Sources)
    record << source << ' ';
  record << '\n';
  if (boundary.Policy == BRepClass3d_ImplicitBox::CornerPolicy::Regularized)
    for (const auto& contact : boundary.ExcludedCornerContacts)
      record << contact.X() << ' ' << contact.Y() << ' ' << contact.Z() << '\n';
  record << m << t;
  const auto text = record.str();
  require(text.size() <= MaxBytes && bool(record));
  output.write(text.data(), std::streamsize(text.size()));
  require(bool(output));
}

void BRepClass3d_ImplicitBoxIO::Read(Standard_IStream& input, Boundary& destination)
{
  require(destination.Shape.IsNull() || !destination.Shape.Locked());
  std::string text;
  char        block[4096];
  while (input.read(block, sizeof(block)) || input.gcount())
  {
    require(text.size() + size_t(input.gcount()) <= MaxBytes);
    text.append(block, size_t(input.gcount()));
  }
  require(input.eof() && !input.bad());
  std::istringstream record(text);
  record.imbue(std::locale::classic());
  std::string magic;
  int         version;
  size_t      mathSize, topologySize, faceCount;
  Boundary    result;
  require(bool(record >> magic >> version >> mathSize >> topologySize >> faceCount
               >> result.SeamCount >> result.BoxEdgeRootCount >> result.GeometricBound));
  require(magic == "OCCT_IMPLICIT_BOX_TOPOLOGY" && (version == 1 || version == 2) && mathSize > 0 && mathSize <= 65536
          && topologySize > 0 && topologySize <= MaxBytes && faceCount > 0 && faceCount <= 128
          && result.SeamCount >= 0 && result.SeamCount <= 500000 && result.BoxEdgeRootCount >= 0
          && result.BoxEdgeRootCount <= 500000);
  if (version == 2)
  {
    int policy = 0, contacts = 0;
    require(bool(record >> policy >> result.ExactCornerCount >> result.CornerContourEnds
                        >> result.PlanarCornerContacts >> contacts));
    require(policy == 1 && result.ExactCornerCount >= 0 && result.ExactCornerCount <= 8
            && result.CornerContourEnds >= 0 && result.CornerContourEnds <= 24
            && result.PlanarCornerContacts >= 0 && result.PlanarCornerContacts <= 24
            && contacts >= 0 && contacts <= 8);
    result.Policy = BRepClass3d_ImplicitBox::CornerPolicy::Regularized;
    result.ExcludedCornerContacts.resize(size_t(contacts));
  }
  require(record.get() == '\n');
  result.Sources.resize(faceCount);
  for (int& source : result.Sources)
    require(bool(record >> source));
  require(record.get() == ' ' && record.get() == '\n');
  for (auto& contact : result.ExcludedCornerContacts)
  {
    double x, y, z;
    require(bool(record >> x >> y >> z) && std::isfinite(x) && std::isfinite(y) && std::isfinite(z)
            && record.get() == '\n');
    contact.SetCoord(x, y, z);
  }
  std::string math(mathSize, '\0'), topology(topologySize, '\0');
  require(bool(record.read(math.data(), std::streamsize(mathSize))));
  require(bool(record.read(topology.data(), std::streamsize(topologySize)))
          && record.peek() == EOF);
  std::istringstream mathInput(math), topologyInput(topology);
  mathInput.imbue(std::locale::classic());
  topologyInput.imbue(std::locale::classic());
  TopoDS_Face implicit;
  BRep_ImplicitIO::Read(mathInput, implicit);
  // Primitive I/O preserves transform values, not empty location-chain
  // representation. Canonicalize exact numerical identity ONLY, without
  // tolerance-snapping a nonidentity transform or changing the primitive API.
  const auto raw = occ::down_cast<BRep_TImplicitFace>(implicit.TShape());
  require(!raw.IsNull());
  for (const auto& location : {implicit.Location(), raw->Location()})
    for (int row = 1; row <= 3; ++row)
      for (int column = 1; column <= 4; ++column)
        require(location.Transformation().Value(row, column) == (row == column ? 1. : 0.));
  TopoDS_Face canonical;
  BRep_Builder().MakeImplicitFace(canonical,
                                  raw->Field(),
                                  raw->Minimum(),
                                  raw->Maximum(),
                                  TopLoc_Location(),
                                  raw->Tolerance());
  // Schema-2 mixed records written before explicit face restrictions stored
  // only the box. Migrate its local mathematical domain, not its trim topology
  // or saved seam provenance. New primitive schema-2 records retain exclusions.
  auto exclusions = raw->ExcludedDomains();
  if (result.Policy == BRepClass3d_ImplicitBox::CornerPolicy::Regularized && exclusions.empty())
    exclusions = BRepClass3d_ImplicitBox::SurfaceCornerExclusions(canonical);
  canonical.TShape(new BRep_TImplicitFace(raw->Field(), raw->Minimum(), raw->Maximum(),
                                         raw->Tolerance(), exclusions));
  canonical.Orientation(implicit.Orientation());
  implicit = canonical;
  primitive(implicit);
  TopoDS_Shape payload;
  BRep_Builder builder;
  BRepTools::Read(payload, topologyInput, builder);
  require(!payload.IsNull() && payload.ShapeType() == TopAbs_COMPOUND && !topologyInput.fail());
  topologyInput >> std::ws;
  require(topologyInput.eof() && payload.Location().IsIdentity()
          && payload.Orientation() == TopAbs_FORWARD);
  resetTransientFlags(payload);
  BRepTools::Clean(payload, true);
  size_t index = 0;
  for (TopoDS_Iterator it(payload); it.More(); it.Next(), ++index)
  {
    if (index + 1 < faceCount)
    {
      require(it.Value().ShapeType() == TopAbs_FACE);
      result.Faces.push_back(TopoDS::Face(it.Value()));
    }
    else
    {
      require(index + 1 == faceCount && it.Value().ShapeType() == TopAbs_COMPOUND
              && it.Value().Location().IsIdentity() && it.Value().Orientation() == TopAbs_FORWARD);
      for (TopoDS_Iterator wire(it.Value()); wire.More(); wire.Next())
      {
        require(wire.Value().ShapeType() == TopAbs_WIRE);
        builder.Add(implicit, wire.Value());
      }
    }
  }
  require(index == faceCount);
  result.Faces.push_back(implicit);
  builder.MakeCompound(result.Shape);
  for (const auto& face : result.Faces)
    builder.Add(result.Shape, face);
  validate(result);
  // Reconstruct preview curves from stored geometry, NOT by retracing/rebuilding.
  BRep_ImplicitPlaneSection::Result preview;
  Map                               edges;
  TopExp::MapShapes(implicit, TopAbs_EDGE, edges);
  for (const auto& shape : edges)
  {
    TopLoc_Location location;
    double          first, last;
    const auto      curve = occ::down_cast<Geom_BSplineCurve>(
      BRep_Tool::Curve(TopoDS::Edge(shape), location, first, last));
    BRep_ImplicitPlaneSection::Curve polyline;
    for (int i = 1; i <= curve->NbPoles(); ++i)
      polyline.Points.push_back(curve->Pole(i));
    polyline.Closed = polyline.Points.front().IsEqual(polyline.Points.back(), 0);
    preview.Curves.push_back(std::move(polyline));
  }
  result.Sections.push_back(std::move(preview));
  result.SeamBoundVerified = false;
  destination              = std::move(result);
}

#undef require
