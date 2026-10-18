// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitSolidIO.hxx>
#include <BRepClass3d_ImplicitBoxIO.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>
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
using Ops                 = BRepClass3d_ImplicitBox;
using Map                 = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
constexpr size_t MaxBytes = 32 * 1024 * 1024, MaxSourceBytes = 16 * 1024 * 1024;

void require(bool condition)
{
  if (!condition)
    throw Standard_DomainError("Invalid native solid topology record or unsupported structure");
}

std::vector<TopoDS_Shape> children(const TopoDS_Shape& shape)
{
  std::vector<TopoDS_Shape> result;
  for (TopoDS_Iterator it(shape, false, false); it.More(); it.Next())
    result.push_back(it.Value());
  return result;
}

const BRep_TImplicitFace& implicit(const TopoDS_Face& face)
{
  const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  require(!data.IsNull() && typeid(*data) == typeid(BRep_TImplicitFace) && !data->Field().IsNull()
          && typeid(*data->Field()) == typeid(Geom_ImplicitGyroid) && face.Location().IsIdentity()
          && data->Location().IsIdentity()
          && (face.Orientation() == TopAbs_FORWARD || face.Orientation() == TopAbs_REVERSED));
  return *data;
}

void identity(const TopoDS_Shape& shape, TopAbs_ShapeEnum type)
{
  require(!shape.IsNull() && shape.ShapeType() == type && shape.Location().IsIdentity()
          && shape.Orientation() == TopAbs_FORWARD);
}

// Structural graph checks only. Field topology, component separator coverage
// and geometric tolerance validity require the independent numerical pipeline.
void checkShell(const TopoDS_Shell& shell, const std::vector<TopoDS_Face>& faces)
{
  identity(shell, TopAbs_SHELL);
  const auto actual = children(shell);
  require(!actual.empty() && actual.size() == faces.size());
  Map edges;
  TopExp::MapShapes(shell, TopAbs_EDGE, edges);
  std::vector<std::vector<int>> owners(size_t(edges.Extent()));
  std::vector<int>              balance(size_t(edges.Extent()));
  for (size_t f = 0; f < faces.size(); ++f)
  {
    require(actual[f].IsEqual(faces[f]));
    for (TopExp_Explorer it(faces[f], TopAbs_EDGE); it.More(); it.Next())
    {
      require(it.Current().Orientation() == TopAbs_FORWARD
              || it.Current().Orientation() == TopAbs_REVERSED);
      const size_t e = size_t(edges.FindIndex(it.Current()) - 1);
      require(e < owners.size());
      owners[e].push_back(int(f));
      balance[e] += it.Current().Orientation() == TopAbs_FORWARD ? 1 : -1;
    }
  }
  std::vector<std::vector<int>> adjacency(faces.size());
  for (size_t e = 0; e < owners.size(); ++e)
  {
    require(owners[e].size() == 2 && owners[e][0] != owners[e][1] && balance[e] == 0);
    adjacency[size_t(owners[e][0])].push_back(owners[e][1]);
    adjacency[size_t(owners[e][1])].push_back(owners[e][0]);
  }
  std::vector<bool> visited(faces.size(), false);
  std::vector<int>  stack{0};
  visited[0] = true;
  while (!stack.empty())
  {
    const int f = stack.back();
    stack.pop_back();
    for (int next : adjacency[size_t(f)])
      if (!visited[size_t(next)])
      {
        visited[size_t(next)] = true;
        stack.push_back(next);
      }
  }
  require(std::all_of(visited.begin(), visited.end(), [](bool value) { return value; }));
}

struct Descriptor
{
  std::vector<std::vector<int>> Wires;
  double                        MaximumTolerance = 0;
};

Descriptor describe(const Ops::NativeSolids& value)
{
  require(value.Constructed && value.Boundary.ShellsConstructed && std::isfinite(value.Tolerance)
          && value.Tolerance > 0 && value.Solids.size() <= 128);
  const auto& b = value.Boundary;
  require(!b.Source.Faces.empty());
  const auto& original    = b.Source.Faces.back();
  const auto& fieldData   = implicit(original);
  const auto& field       = *static_cast<const Geom_ImplicitGyroid*>(fieldData.Field().get());
  const auto  sourceWires = children(original);
  Map         sourceWireMap;
  for (size_t w = 0; w < sourceWires.size(); ++w)
    require(sourceWires[w].ShapeType() == TopAbs_WIRE && !sourceWires[w].Free()
            && sourceWireMap.Add(sourceWires[w]) == int(w + 1));
  const size_t caps = b.Source.Faces.size() - 1, components = b.ComponentDomains.size();
  require(caps + components <= 128 && sourceWires.size() <= 500000
          && b.Faces.size() == caps + components && b.FaceSources.size() == b.Faces.size()
          && b.FaceComponents.size() == b.Faces.size() && b.FaceShells.size() == b.Faces.size()
          && b.Shells.size() <= 128 && b.ShellRegions.size() == b.Shells.size()
          && b.ShellEuler.size() == b.Shells.size());
  identity(b.Shape, TopAbs_COMPOUND);
  identity(value.Shape, TopAbs_COMPOUND);
  const auto shells = children(b.Shape), solids = children(value.Shape);
  require(shells.size() == b.Shells.size() && solids.size() == value.Solids.size());
  Descriptor result;
  result.Wires.resize(components);
  std::vector<int>                      used(sourceWires.size(), 0);
  std::vector<std::vector<TopoDS_Face>> groups(b.Shells.size());
  for (size_t f = 0; f < b.Faces.size(); ++f)
  {
    const auto& face = b.Faces[f];
    require(b.FaceShells[f] >= 0 && size_t(b.FaceShells[f]) < groups.size());
    groups[size_t(b.FaceShells[f])].push_back(face);
    if (f < caps)
      // Avoid changing manually freed children when the source codec builds
      // its temporary payload. Factory-owned caps/wires are already frozen.
      require(!face.Free() && face.IsEqual(b.Source.Faces[f]) && b.FaceSources[f] == int(f)
              && b.FaceComponents[f] == -1);
    else
    {
      const size_t c    = f - caps;
      const auto&  data = implicit(face);
      const auto&  g    = *static_cast<const Geom_ImplicitGyroid*>(data.Field().get());
      require(b.FaceSources[f] == int(caps) && b.FaceComponents[f] == int(c)
              && g.Period() == field.Period() && g.IsoValue() == field.IsoValue()
              && face.Orientation() == original.Orientation()
              && data.Tolerance() == fieldData.Tolerance());
      for (int a = 1; a <= 3; ++a)
      {
        const double lo = b.ComponentDomains[c][0].Coord(a), hi = b.ComponentDomains[c][1].Coord(a);
        require(std::isfinite(lo) && std::isfinite(hi) && lo < hi
                && lo >= fieldData.Minimum().Coord(a) && hi <= fieldData.Maximum().Coord(a)
                && data.Minimum().Coord(a) == lo && data.Maximum().Coord(a) == hi);
      }
      for (const auto& wire : children(face))
      {
        require(wire.ShapeType() == TopAbs_WIRE);
        const size_t index = size_t(sourceWireMap.FindIndex(wire) - 1);
        require(index < sourceWires.size() && sourceWires[index].IsEqual(wire));
        require(++used[index] == 1);
        result.Wires[c].push_back(int(index));
      }
    }
  }
  require(std::all_of(used.begin(), used.end(), [](int count) { return count == 1; }));
  for (size_t s = 0; s < b.Shells.size(); ++s)
  {
    require(shells[s].IsEqual(b.Shells[s]) && b.ShellRegions[s] >= 0
            && size_t(b.ShellRegions[s]) < value.Solids.size());
    checkShell(b.Shells[s], groups[s]);
  }
  for (size_t r = 0; r < value.Solids.size(); ++r)
  {
    identity(value.Solids[r], TopAbs_SOLID);
    require(solids[r].IsEqual(value.Solids[r]));
    const auto actual = children(value.Solids[r]);
    size_t     index  = 0;
    for (size_t s = 0; s < b.Shells.size(); ++s)
      if (b.ShellRegions[s] == int(r))
      {
        require(index < actual.size() && actual[index].IsEqual(b.Shells[s]));
        ++index;
      }
    require(index > 0 && index == actual.size());
  }
  for (const auto& face : b.Faces)
    result.MaximumTolerance = std::max(result.MaximumTolerance, BRep_Tool::Tolerance(face));
  for (TopExp_Explorer it(value.Shape, TopAbs_EDGE); it.More(); it.Next())
    result.MaximumTolerance =
      std::max(result.MaximumTolerance, BRep_Tool::Tolerance(TopoDS::Edge(it.Current())));
  for (TopExp_Explorer it(value.Shape, TopAbs_VERTEX); it.More(); it.Next())
    result.MaximumTolerance =
      std::max(result.MaximumTolerance, BRep_Tool::Tolerance(TopoDS::Vertex(it.Current())));
  require(std::isfinite(result.MaximumTolerance) && result.MaximumTolerance <= value.Tolerance);
  return result;
}

void label(std::istream& input, const char* expected)
{
  std::string actual;
  require(bool(input >> actual) && actual == expected && input.get() == '\n');
}
} // namespace

void BRepClass3d_ImplicitSolidIO::Write(const Ops::NativeSolids& value, Standard_OStream& output)
{
  const auto         descriptor = describe(value);
  std::ostringstream source;
  BRepClass3d_ImplicitBoxIO::Write(value.Boundary.Source, source);
  const auto         sourceText = source.str();
  const auto&        b          = value.Boundary;
  std::ostringstream record;
  record.imbue(std::locale::classic());
  record << std::setprecision(std::numeric_limits<double>::max_digits10)
         << "OCCT_IMPLICIT_SOLIDS 1 " << sourceText.size() << ' ' << b.Faces.size() << ' '
         << b.ComponentDomains.size() << ' ' << b.Shells.size() << ' ' << value.Solids.size() << ' '
         << value.Tolerance << "\nFACES\n";
  for (size_t f = 0; f < b.Faces.size(); ++f)
    record << b.FaceSources[f] << ' ' << b.FaceComponents[f] << ' ' << b.FaceShells[f] << '\n';
  record << "COMPONENTS\n";
  for (size_t c = 0; c < b.ComponentDomains.size(); ++c)
  {
    for (const auto& p : b.ComponentDomains[c])
      record << p.X() << ' ' << p.Y() << ' ' << p.Z() << ' ';
    record << descriptor.Wires[c].size();
    for (int wire : descriptor.Wires[c])
      record << ' ' << wire;
    record << '\n';
  }
  record << "SHELLS\n";
  for (size_t s = 0; s < b.Shells.size(); ++s)
    record << b.ShellRegions[s] << ' ' << b.ShellEuler[s] << '\n';
  record << "SOURCE\n" << sourceText;
  const auto text = record.str();
  require(bool(record) && sourceText.size() <= MaxSourceBytes && text.size() <= MaxBytes);
  output.write(text.data(), std::streamsize(text.size()));
  require(bool(output));
}

void BRepClass3d_ImplicitSolidIO::Read(Standard_IStream& input, Ops::NativeSolids& destination)
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
  std::string       magic;
  int               version;
  size_t            sourceSize, faces, components, shells, regions;
  Ops::NativeSolids result;
  require(bool(record >> magic >> version >> sourceSize >> faces >> components >> shells >> regions
               >> result.Tolerance));
  require(magic == "OCCT_IMPLICIT_SOLIDS" && version == 1 && sourceSize > 0
          && sourceSize <= MaxSourceBytes && faces <= 128 && components <= faces && shells <= 128
          && regions <= shells && std::isfinite(result.Tolerance) && result.Tolerance > 0
          && record.get() == '\n');
  auto& b = result.Boundary;
  b.FaceSources.resize(faces);
  b.FaceComponents.resize(faces);
  b.FaceShells.resize(faces);
  label(record, "FACES");
  for (size_t f = 0; f < faces; ++f)
    require(bool(record >> b.FaceSources[f] >> b.FaceComponents[f] >> b.FaceShells[f])
            && record.get() == '\n');
  label(record, "COMPONENTS");
  b.ComponentDomains.resize(components);
  std::vector<std::vector<int>> wires(components);
  size_t                        totalWires = 0;
  for (size_t c = 0; c < components; ++c)
  {
    for (auto& p : b.ComponentDomains[c])
    {
      double x, y, z;
      require(bool(record >> x >> y >> z) && std::isfinite(x) && std::isfinite(y)
              && std::isfinite(z));
      p.SetCoord(x, y, z);
    }
    size_t count;
    require(bool(record >> count) && count <= 500000 && (totalWires += count) <= 500000);
    wires[c].resize(count);
    for (int& wire : wires[c])
      require(bool(record >> wire) && wire >= 0 && wire < 500000);
    require(record.get() == '\n');
  }
  label(record, "SHELLS");
  b.ShellRegions.resize(shells);
  b.ShellEuler.resize(shells);
  for (size_t s = 0; s < shells; ++s)
    require(bool(record >> b.ShellRegions[s] >> b.ShellEuler[s]) && record.get() == '\n');
  label(record, "SOURCE");
  std::string source(sourceSize, '\0');
  require(bool(record.read(source.data(), std::streamsize(sourceSize))) && record.peek() == EOF);
  std::istringstream sourceInput(source);
  BRepClass3d_ImplicitBoxIO::Read(sourceInput, b.Source);
  require(!b.Source.Faces.empty() && faces == b.Source.Faces.size() - 1 + components);
  const auto&  original    = b.Source.Faces.back();
  const auto&  data        = implicit(original);
  const auto   sourceWires = children(original);
  BRep_Builder builder;
  b.Faces.assign(b.Source.Faces.begin(), b.Source.Faces.end() - 1);
  for (size_t c = 0; c < components; ++c)
  {
    TopoDS_Face face;
    builder.MakeImplicitFace(face,
                             data.Field(),
                             b.ComponentDomains[c][0],
                             b.ComponentDomains[c][1],
                             TopLoc_Location(),
                             data.Tolerance());
    for (int wire : wires[c])
    {
      require(size_t(wire) < sourceWires.size());
      builder.Add(face, sourceWires[size_t(wire)]);
    }
    // Assign local stored wires while the parent is forward. BRep_Builder::Add
    // composes parent orientation; setting a cut face reversed first would
    // invert the raw source-wire orientations a second time.
    face.Orientation(original.Orientation());
    face.Closed(false);
    b.Faces.push_back(face);
  }
  builder.MakeCompound(b.Shape);
  for (size_t s = 0; s < shells; ++s)
  {
    TopoDS_Shell shell;
    builder.MakeShell(shell);
    for (size_t f = 0; f < faces; ++f)
      if (b.FaceShells[f] == int(s))
        builder.Add(shell, b.Faces[f]);
    shell.Closed(false);
    b.Shells.push_back(shell);
    builder.Add(b.Shape, shell);
  }
  builder.MakeCompound(result.Shape);
  for (size_t r = 0; r < regions; ++r)
  {
    TopoDS_Solid solid;
    builder.MakeSolid(solid);
    for (size_t s = 0; s < shells; ++s)
      if (b.ShellRegions[s] == int(r))
        builder.Add(solid, b.Shells[s]);
    solid.Closed(false);
    result.Solids.push_back(solid);
    builder.Add(result.Shape, solid);
  }
  result.Constructed = b.ShellsConstructed = true;
  result.MaximumStoredTolerance            = describe(result).MaximumTolerance;
  result.Reason = b.Reason = "restored_native_solids_unverified";
  result.Geometry.Reason   = "not_rechecked_after_restore";
  destination              = std::move(result);
}
