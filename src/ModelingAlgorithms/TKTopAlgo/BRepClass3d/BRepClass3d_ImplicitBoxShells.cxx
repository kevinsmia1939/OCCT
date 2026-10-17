// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitBox.hxx>
#include <BRepClass3d_ImplicitBoxIO.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_BSplineCurve.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <sstream>

namespace
{
using Box = std::array<gp_Pnt, 2>;
using Map = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

void require(bool condition, const char* reason)
{
  if (!condition)
    throw Standard_NotImplemented(reason);
}

struct Partition
{
  const Geom_ImplicitGyroid&             Field;
  const Box                              Original;
  const std::vector<Box>&                Proposals;
  std::vector<Box>&                      Domains;
  BRepClass3d_ImplicitBox::NativeShells& Result;
  int                                    Resolution, MaximumDepth;
  size_t                                 MaximumCells;
  double                                 L, H, Reserve, FieldReserve;

  double linear(const gp_Pnt& p) const
  {
    int    index[3], order[3] = {0, 1, 2};
    double fraction[3];
    for (int a = 0; a < 3; ++a)
    {
      const double position = (p.Coord(a + 1) - Original[0].Coord(a + 1)) * Resolution
                              / (Original[1].Coord(a + 1) - Original[0].Coord(a + 1));
      index[a]              = std::clamp(int(std::floor(position)), 0, Resolution - 1);
      fraction[a]           = std::clamp(position - index[a], 0., 1.);
    }
    std::sort(order, order + 3, [&](int a, int b) { return fraction[a] > fraction[b]; });
    auto sample = [&]() {
      gp_Pnt q;
      for (int a = 0; a < 3; ++a)
        q.SetCoord(a + 1,
                   index[a] == Resolution
                     ? Original[1].Coord(a + 1)
                     : Original[0].Coord(a + 1)
                         + (Original[1].Coord(a + 1) - Original[0].Coord(a + 1)) * index[a]
                             / Resolution);
      return Field.Value(q);
    };
    double result = (1 - fraction[order[0]]) * sample();
    ++index[order[0]];
    result += (fraction[order[0]] - fraction[order[1]]) * sample();
    ++index[order[1]];
    result += (fraction[order[1]] - fraction[order[2]]) * sample();
    ++index[order[2]];
    return result + fraction[order[2]] * sample();
  }

  bool exclude(const Box& plane, int fixed, int depth, int& sign)
  {
    if (Result.SeparatorCells >= MaximumCells)
      return false;
    ++Result.SeparatorCells;
    const gp_Pnt p((plane[0].XYZ() + plane[1].XYZ()) * .5);
    auto         gradient = Field.Gradient(p);
    gradient.SetCoord(fixed + 1, 0);
    const double radius = plane[0].Distance(plane[1]) * .5 + 8 * Reserve;
    const double variation =
      std::min(L * radius, gradient.Magnitude() * radius + .5 * H * radius * radius)
        * (1 + 128 * std::numeric_limits<double>::epsilon())
      + FieldReserve;
    const double value = Field.Value(p);
    if (std::abs(value) > variation)
    {
      const int current = value > 0 ? 1 : -1;
      if (sign && sign != current)
        return false;
      sign = current;
      return true;
    }
    if (depth >= MaximumDepth)
      return false;
    int split = (fixed + 1) % 3;
    for (int a = 0; a < 3; ++a)
      if (a != fixed
          && plane[1].Coord(a + 1) - plane[0].Coord(a + 1)
               > plane[1].Coord(split + 1) - plane[0].Coord(split + 1))
        split = a;
    const double middle = p.Coord(split + 1);
    if (!(middle > plane[0].Coord(split + 1) && middle < plane[1].Coord(split + 1)))
      return false;
    auto left = plane, right = plane;
    left[1].SetCoord(split + 1, middle);
    right[0].SetCoord(split + 1, middle);
    return exclude(left, fixed, depth + 1, sign) && exclude(right, fixed, depth + 1, sign);
  }

  bool divide(std::vector<int> components, const Box& domain)
  {
    if (components.empty())
      return true;
    if (components.size() == 1)
    {
      Domains[size_t(components[0])] = domain;
      return true;
    }
    for (int axis = 0; axis < 3; ++axis)
    {
      std::sort(components.begin(), components.end(), [&](int a, int b) {
        return Proposals[size_t(a)][0].Coord(axis + 1) < Proposals[size_t(b)][0].Coord(axis + 1);
      });
      double previous = Proposals[size_t(components[0])][1].Coord(axis + 1);
      for (size_t split = 1; split < components.size(); ++split)
      {
        const double next = Proposals[size_t(components[split])][0].Coord(axis + 1);
        if (next - previous > 32 * Reserve)
          for (double fraction : {.5, .25, .75})
          {
            const double middle = previous + (next - previous) * fraction;
            Box          plane  = domain;
            plane[0].SetCoord(axis + 1, middle);
            plane[1].SetCoord(axis + 1, middle);
            const gp_Pnt p((plane[0].XYZ() + plane[1].XYZ()) * .5);
            const double v = linear(p);
            if (!(std::abs(v) > FieldReserve))
              continue;
            int sign = 0;
            if (exclude(plane, axis, 0, sign) && sign * v > FieldReserve)
            {
              // The proposal gap excludes ALL reference zero triangles in this
              // subdomain. L has a constant sign on the connected plane. F's
              // same strict sign excludes zeros throughout the whole homotopy.
              ++Result.Separators;
              auto left = domain, right = domain;
              left[1].SetCoord(axis + 1, middle);
              right[0].SetCoord(axis + 1, middle);
              return divide(std::vector<int>(components.begin(), components.begin() + split), left)
                     && divide(std::vector<int>(components.begin() + split, components.end()),
                               right);
            }
            if (Result.SeparatorCells >= MaximumCells)
              return false;
          }
        previous = std::max(previous, Proposals[size_t(components[split])][1].Coord(axis + 1));
      }
    }
    return false;
  }
};

void validateShell(const TopoDS_Shell&             shell,
                   const std::vector<int>&         faceIndices,
                   const std::vector<TopoDS_Face>& faces,
                   const std::vector<int>&         faceEuler,
                   int                             expectedEuler)
{
  Map edges, vertices;
  TopExp::MapShapes(shell, TopAbs_EDGE, edges);
  TopExp::MapShapes(shell, TopAbs_VERTEX, vertices);

  struct Use
  {
    int              Count = 0, Balance = 0;
    std::vector<int> Faces;
  };

  std::vector<Use>                             uses(size_t(edges.Extent()));
  std::vector<std::map<int, std::vector<int>>> links(size_t(vertices.Extent()));
  std::vector<std::vector<int>>                adjacent(faceIndices.size());
  int                                          euler = vertices.Extent() - edges.Extent();
  for (size_t f = 0; f < faceIndices.size(); ++f)
  {
    const int index = faceIndices[f];
    euler += faceEuler[size_t(index)];
    const auto& face = faces[size_t(index)];
    for (TopoDS_Iterator wire(face); wire.More(); wire.Next())
    {
      // Incoming/outgoing germs avoid assuming TopoDS child-list traversal is
      // geometric wire order. Closed edges have two distinct endpoint germs.
      std::map<int, std::pair<std::vector<int>, std::vector<int>>> corners;
      for (TopoDS_Iterator it(wire.Value()); it.More(); it.Next())
      {
        const auto edge = TopoDS::Edge(it.Value());
        const int  id   = edges.FindIndex(edge) - 1;
        require(id >= 0, "Native shell edge map incomplete");
        auto& use = uses[size_t(id)];
        ++use.Count;
        use.Faces.push_back(int(f));
        const bool forward = edge.Orientation() == TopAbs_FORWARD;
        require(forward || edge.Orientation() == TopAbs_REVERSED,
                "Native shell edge orientation unsupported");
        use.Balance += forward ? 1 : -1;
        TopoDS_Vertex ends[2];
        TopExp::Vertices(TopoDS::Edge(edge.Oriented(TopAbs_FORWARD)), ends[0], ends[1]);
        const int begin = forward ? 0 : 1, end = 1 - begin;
        const int a = vertices.FindIndex(ends[begin]) - 1, b = vertices.FindIndex(ends[end]) - 1;
        require(a >= 0 && b >= 0, "Native shell vertex map incomplete");
        corners[a].second.push_back(2 * id + begin);
        corners[b].first.push_back(2 * id + end);
      }
      for (const auto& entry : corners)
      {
        require(entry.second.first.size() == 1 && entry.second.second.size() == 1,
                "Native shell wire has ambiguous vertex corners");
        const int a = entry.second.first[0], b = entry.second.second[0];
        require(a != b, "Native shell has a collapsed corner germ");
        links[size_t(entry.first)][a].push_back(b);
        links[size_t(entry.first)][b].push_back(a);
      }
    }
  }
  require(euler == expectedEuler, "Native/reference shell Euler characteristics differ");
  for (const auto& use : uses)
  {
    require(use.Count == 2 && use.Balance == 0 && use.Faces[0] != use.Faces[1],
            "Native shell edges are not paired between opposite face uses");
    adjacent[size_t(use.Faces[0])].push_back(use.Faces[1]);
    adjacent[size_t(use.Faces[1])].push_back(use.Faces[0]);
  }
  for (const auto& link : links)
  {
    require(!link.empty(), "Native shell has an unused vertex");
    for (const auto& germ : link)
      require(germ.second.size() == 2, "Native shell vertex link is not cyclic");
    std::set<int>    visited;
    std::vector<int> todo{link.begin()->first};
    while (!todo.empty())
    {
      const int node = todo.back();
      todo.pop_back();
      if (!visited.insert(node).second)
        continue;
      for (int next : link.at(node))
        todo.push_back(next);
    }
    require(visited.size() == link.size(), "Native shell has a pinched vertex link");
  }
  require(!faceIndices.empty(), "Reference boundary shell has no native faces");
  std::set<int>    visited;
  std::vector<int> todo{0};
  while (!todo.empty())
  {
    const int node = todo.back();
    todo.pop_back();
    if (!visited.insert(node).second)
      continue;
    for (int next : adjacent[size_t(node)])
      todo.push_back(next);
  }
  require(visited.size() == faceIndices.size(), "Native shell faces are disconnected");
}
} // namespace

BRepClass3d_ImplicitBox::NativeShells BRepClass3d_ImplicitBox::BuildNativeShells(
  const Boundary& source,
  int             resolution,
  int             maximumDepth,
  size_t          maximumCells,
  size_t          maximumTests,
  size_t          maximumSeparatorCells)
{
  if (!maximumSeparatorCells || maximumSeparatorCells > 10000000)
    throw Standard_DomainError("Invalid native component separator budget");
  NativeShells       result;
  std::ostringstream output;
  BRepClass3d_ImplicitBoxIO::Write(source, output);
  std::istringstream input(output.str());
  Boundary           copy;
  BRepClass3d_ImplicitBoxIO::Read(input, copy);
  const auto attachment =
    CheckTrimAttachment(copy, resolution, maximumDepth, maximumCells, maximumTests);
  if (!attachment.Resolved)
  {
    result.Reason = "trim_attachment_unresolved: " + attachment.Reason;
    return result;
  }
  const auto&  retained   = attachment.RetainedCaps;
  const auto&  reference  = retained.Reference;
  const auto&  original   = copy.Faces.back();
  const auto   data       = occ::down_cast<BRep_TImplicitFace>(original.TShape());
  const auto&  field      = *static_cast<const Geom_ImplicitGyroid*>(data->Field().get());
  const int    components = attachment.ReferenceSurfaceComponentCount;
  const size_t caps       = copy.Faces.size() - 1;
  try
  {
    require(caps + size_t(components) <= 128, "Native component face budget exhausted");
    const Box         originalDomain{data->Minimum(), data->Maximum()};
    std::vector<Box>  proposals(size_t(components), originalDomain);
    std::vector<bool> initialized(size_t(components), false);
    auto              extend = [&](int c, const gp_Pnt& p) {
      auto& box = proposals[size_t(c)];
      if (!initialized[size_t(c)])
      {
        box                    = {p, p};
        initialized[size_t(c)] = true;
      }
      for (int axis = 1; axis <= 3; ++axis)
      {
        box[0].SetCoord(axis, std::min(box[0].Coord(axis), p.Coord(axis)));
        box[1].SetCoord(axis, std::max(box[1].Coord(axis), p.Coord(axis)));
      }
    };
    for (size_t i = 0; i < reference.Surface.Nodes.size(); ++i)
      extend(attachment.SurfaceNodeComponents[i], reference.Surface.Nodes[i]);
    std::vector<TopoDS_Shape> wires;
    for (TopoDS_Iterator it(original, false, false); it.More(); it.Next())
      wires.push_back(it.Value());
    require(wires.size() == attachment.WireSurfaceComponents.size(),
            "Native face wire component map incomplete");
    for (size_t i = 0; i < wires.size(); ++i)
      for (TopExp_Explorer it(wires[i], TopAbs_EDGE); it.More(); it.Next())
      {
        TopLoc_Location location;
        double          first, last;
        const auto      curve = occ::down_cast<Geom_BSplineCurve>(
          BRep_Tool::Curve(TopoDS::Edge(it.Current()), location, first, last));
        for (int pole = 1; pole <= curve->NbPoles(); ++pole)
          extend(attachment.WireSurfaceComponents[i], curve->Pole(pole));
      }
    double magnitude = 1;
    for (int axis = 1; axis <= 3; ++axis)
      magnitude = std::max({magnitude,
                            std::abs(originalDomain[0].Coord(axis)),
                            std::abs(originalDomain[1].Coord(axis))});
    const double reserve = 8192 * std::numeric_limits<double>::epsilon() * magnitude;
    const double k = 2 * std::acos(-1.) / field.Period(), l = field.LipschitzBound();
    result.ComponentDomains.resize(size_t(components), originalDomain);
    Partition        partition{field,
                               originalDomain,
                               proposals,
                               result.ComponentDomains,
                               result,
                               resolution,
                               maximumDepth,
                               maximumSeparatorCells,
                               l,
                               4 * k * k,
                               reserve,
                               4096 * std::numeric_limits<double>::epsilon()
                                   * (1 + std::abs(field.IsoValue()))
                                 + l * reserve};
    std::vector<int> ids;
    for (int c = 0; c < components; ++c)
      ids.push_back(c);
    require(partition.divide(std::move(ids), originalDomain),
            result.SeparatorCells >= maximumSeparatorCells
              ? "Component separator budget exhausted"
              : "Zero-free component separator unresolved");
    BRep_Builder     builder;
    std::vector<int> faceEuler;
    for (size_t cap = 0; cap < caps; ++cap)
    {
      result.Faces.push_back(copy.Faces[cap]);
      result.FaceSources.push_back(int(cap));
      result.FaceComponents.push_back(-1);
      result.FaceShells.push_back(retained.ReferenceShells[cap]);
      int boundaries = 0;
      for (TopoDS_Iterator it(copy.Faces[cap], false, false); it.More(); it.Next())
        ++boundaries;
      faceEuler.push_back(2 - boundaries);
    }
    for (int c = 0; c < components; ++c)
    {
      const auto& domain = result.ComponentDomains[size_t(c)];
      for (int axis = 1; axis <= 3; ++axis)
        require(proposals[size_t(c)][0].Coord(axis) >= domain[0].Coord(axis) - reserve
                  && proposals[size_t(c)][1].Coord(axis) <= domain[1].Coord(axis) + reserve,
                "Native trim geometry exceeds its component domain");
      TopoDS_Face face;
      builder.MakeImplicitFace(face,
                               data->Field(),
                               domain[0],
                               domain[1],
                               TopLoc_Location(),
                               data->Tolerance());
      int boundaries = 0;
      for (size_t i = 0; i < wires.size(); ++i)
        if (attachment.WireSurfaceComponents[i] == c)
        {
          builder.Add(face, wires[i]);
          ++boundaries;
        }
      face.Orientation(original.Orientation());
      face.Closed(boundaries == 0);
      result.Faces.push_back(face);
      result.FaceSources.push_back(int(caps));
      result.FaceComponents.push_back(c);
      result.FaceShells.push_back(attachment.SurfaceComponentShells[size_t(c)]);
      faceEuler.push_back(attachment.SurfaceComponentEuler[size_t(c)]);
    }
    result.FacesAttached = true;
    std::vector<std::vector<int>> groups(reference.ComponentEuler.size());
    for (size_t f = 0; f < result.Faces.size(); ++f)
    {
      require(result.FaceShells[f] >= 0 && size_t(result.FaceShells[f]) < groups.size(),
              "Native face shell attachment invalid");
      groups[size_t(result.FaceShells[f])].push_back(int(f));
    }
    std::vector<TopoDS_Shell> shells;
    for (size_t s = 0; s < groups.size(); ++s)
    {
      TopoDS_Shell shell;
      builder.MakeShell(shell);
      for (int f : groups[s])
        builder.Add(shell, result.Faces[size_t(f)]);
      validateShell(shell, groups[s], result.Faces, faceEuler, reference.ComponentEuler[s]);
      shell.Closed(true);
      shells.push_back(shell);
    }
    builder.MakeCompound(result.Shape);
    for (const auto& shell : shells)
      builder.Add(result.Shape, shell);
    result.Shells            = std::move(shells);
    result.ShellRegions      = reference.ComponentRegions;
    result.ShellEuler        = reference.ComponentEuler;
    result.Source            = std::move(copy);
    result.TopologyValidated = true;
    result.ShellsConstructed = true;
    result.Reason            = "native_shells_constructed";
  }
  catch (const Standard_NotImplemented& error)
  {
    result.Reason = error.what();
  }
  return result;
}
