// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepMesh_ImplicitBoundary.hxx>
#include <BRepClass3d_ImplicitBoxIO.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>

namespace
{
using Triangle = std::array<int, 3>;
using Mesh     = BRepMesh_ImplicitBoundary::Result;
using Report   = BRepMesh_ImplicitBoundary::Coverage;

struct Box
{
  gp_Pnt Low, High;
};

gp_Pnt midpoint(const gp_Pnt& a, const gp_Pnt& b)
{
  return gp_Pnt((a.XYZ() + b.XYZ()) * .5);
}

double radius(const Box& b)
{
  return gp_Vec(b.Low, b.High).Magnitude() * .5;
}

// Upper distances always refer to an explicit convex point on a triangle.
// A missed/ill-conditioned projection merely gives a looser upper bound.
double triangleUpper(const gp_Pnt& p, const gp_Pnt& a, const gp_Pnt& b, const gp_Pnt& c)
{
  double       best      = std::min({p.Distance(a), p.Distance(b), p.Distance(c)});
  const gp_Pnt points[3] = {a, b, c};
  for (int i = 0; i < 3; ++i)
  {
    const auto&  x = points[i];
    const auto&  y = points[(i + 1) % 3];
    const gp_Vec edge(x, y);
    if (edge.SquareMagnitude() > 0)
    {
      const double t = std::clamp(gp_Vec(x, p).Dot(edge) / edge.SquareMagnitude(), 0., 1.);
      best           = std::min(best, p.Distance(x.Translated(edge * t)));
    }
  }
  const gp_Vec ab(a, b), ac(a, c), ap(a, p), n = ab.Crossed(ac);
  const double square = n.SquareMagnitude();
  if (square > 0 && std::isfinite(square))
  {
    const double u = ap.Crossed(ac).Dot(n) / square, v = ab.Crossed(ap).Dot(n) / square;
    if (std::isfinite(u) && std::isfinite(v) && u >= 0 && v >= 0 && u + v <= 1)
      best = std::min(best, p.Distance(a.Translated(ab * u + ac * v)));
  }
  return best;
}

class Tree
{
  struct Node
  {
    Box    Bounds;
    int    Left = -1, Right = -1;
    size_t First = 0, Last = 0;
  };

  const Mesh&           myMesh;
  std::vector<Triangle> myTriangles;
  std::vector<Node>     myNodes;
  double                myReserve;

  int build(size_t first, size_t last)
  {
    Node node;
    node.First      = first;
    node.Last       = last;
    node.Bounds.Low = node.Bounds.High = myMesh.Nodes[size_t(myTriangles[first][0])];
    for (size_t i = first; i < last; ++i)
      for (int index : myTriangles[i])
        for (int axis = 1; axis <= 3; ++axis)
        {
          const double v = myMesh.Nodes[size_t(index)].Coord(axis);
          node.Bounds.Low.SetCoord(axis, std::min(node.Bounds.Low.Coord(axis), v));
          node.Bounds.High.SetCoord(axis, std::max(node.Bounds.High.Coord(axis), v));
        }
    const int id = int(myNodes.size());
    myNodes.push_back(node);
    if (last - first > 8)
    {
      int axis = 1;
      for (int j = 2; j <= 3; ++j)
        if (node.Bounds.High.Coord(j) - node.Bounds.Low.Coord(j)
            > node.Bounds.High.Coord(axis) - node.Bounds.Low.Coord(axis))
          axis = j;
      const size_t mid    = first + (last - first) / 2;
      const auto   center = [&](const Triangle& t) {
        return myMesh.Nodes[size_t(t[0])].Coord(axis) + myMesh.Nodes[size_t(t[1])].Coord(axis)
               + myMesh.Nodes[size_t(t[2])].Coord(axis);
      };
      std::nth_element(myTriangles.begin() + first,
                       myTriangles.begin() + mid,
                       myTriangles.begin() + last,
                       [&](const Triangle& a, const Triangle& b) { return center(a) < center(b); });
      myNodes[size_t(id)].Left  = build(first, mid);
      myNodes[size_t(id)].Right = build(mid, last);
    }
    return id;
  }

  double lower(const gp_Pnt& p, const Box& b) const
  {
    double d[3];
    for (int j = 1; j <= 3; ++j)
      d[j - 1] = std::max({0., b.Low.Coord(j) - p.Coord(j), p.Coord(j) - b.High.Coord(j)});
    return std::max(0., std::hypot(d[0], d[1], d[2]) - myReserve);
  }

  void visit(int index, const gp_Pnt& p, double& best) const
  {
    const auto& node = myNodes[size_t(index)];
    if (lower(p, node.Bounds) > best)
      return;
    if (node.Left < 0)
    {
      for (size_t i = node.First; i < node.Last; ++i)
      {
        const auto& t = myTriangles[i];
        best          = std::min(best,
                                 triangleUpper(p,
                                               myMesh.Nodes[size_t(t[0])],
                                               myMesh.Nodes[size_t(t[1])],
                                               myMesh.Nodes[size_t(t[2])])
                                   + myReserve);
      }
    }
    else
    {
      int a = node.Left, b = node.Right;
      if (lower(p, myNodes[size_t(a)].Bounds) > lower(p, myNodes[size_t(b)].Bounds))
        std::swap(a, b);
      visit(a, p, best);
      visit(b, p, best);
    }
  }

public:
  Tree(const Mesh& mesh, std::vector<Triangle> triangles, double reserve)
      : myMesh(mesh),
        myTriangles(std::move(triangles)),
        myReserve(reserve)
  {
    if (!myTriangles.empty())
      build(0, myTriangles.size());
  }

  double upper(const gp_Pnt& p) const
  {
    double best = std::numeric_limits<double>::infinity();
    if (!myNodes.empty())
      visit(0, p, best);
    return best;
  }
};

struct Checker
{
  const BRepClass3d_ImplicitBox::Boundary& Boundary;
  const Geom_ImplicitGyroid&               Field;
  gp_Pnt                                   Low, High;
  double Tolerance, Reserve, FieldReserve, GradientReserve, H, L, Direction;
  int    MaximumDepth;
  size_t MaximumCells, MaximumPatches;
  Report Result;

  double capWitness(const gp_Pnt& p, double available, int requestedSide) const
  {
    double best = std::numeric_limits<double>::infinity();
    for (int side = 0; side < 6; ++side)
    {
      if (requestedSide >= 0 && requestedSide != side)
        continue;
      const int axis = side / 2 + 1;
      gp_Pnt    base = p;
      for (int j = 1; j <= 3; ++j)
        base.SetCoord(j, std::clamp(p.Coord(j), Low.Coord(j), High.Coord(j)));
      base.SetCoord(axis, side % 2 ? High.Coord(axis) : Low.Coord(axis));
      if (p.Distance(base) + 8 * Reserve > available)
        continue;
      const auto test = [&](const gp_Pnt& q) {
        const double distance = p.Distance(q) + 8 * Reserve;
        if (distance <= available && Direction * Field.Value(q) < -FieldReserve)
          best = std::min(best, distance);
      };
      test(base);
      auto gradient = Field.Gradient(base);
      gradient.SetCoord(axis, 0);
      std::vector<gp_Vec> directions;
      if (gradient.Magnitude() > GradientReserve)
        directions.push_back(gradient / gradient.Magnitude());
      for (int j = 1; j <= 3; ++j)
        if (j != axis)
        {
          gp_Vec direction(0, 0, 0);
          direction.SetCoord(j, 1);
          directions.push_back(direction);
        }
      for (const auto& direction : directions)
        for (double factor : {-.25, .25, -.5, .5, -.9, .9})
        {
          auto q = base.Translated(direction * (available * factor));
          for (int j = 1; j <= 3; ++j)
            q.SetCoord(j, std::clamp(q.Coord(j), Low.Coord(j), High.Coord(j)));
          test(q);
        }
    }
    return best;
  }

  double variation(const gp_Pnt& p, double r) const
  {
    return std::min(L * r, (Field.Gradient(p).Magnitude() + GradientReserve) * r + H * r * r * .5)
             * (1 + 64 * std::numeric_limits<double>::epsilon())
           + FieldReserve;
  }

  bool coverBox(const Box& box, const Tree& tree, int depth, bool cap)
  {
    Result.UnresolvedMinimum = box.Low;
    Result.UnresolvedMaximum = box.High;
    if (Result.SpatialCells >= MaximumCells)
      return false;
    ++Result.SpatialCells;
    const auto   p = midpoint(box.Low, box.High);
    const double r = radius(box) + 4 * Reserve, value = Field.Value(p), error = variation(p, r);
    if ((!cap && std::abs(value) > error) || (cap && Direction * value > error))
      return true;
    const double upper = tree.upper(p) + r;
    if (upper <= Tolerance)
    {
      Result.BoundaryToMeshBound = std::max(Result.BoundaryToMeshBound, upper);
      return true;
    }
    if (depth >= MaximumDepth)
      return false;
    std::vector<Box> children{box};
    for (int axis = 1; axis <= 3; ++axis)
    {
      if (box.Low.Coord(axis) == box.High.Coord(axis))
        continue;
      const double mid = p.Coord(axis);
      if (!(mid > box.Low.Coord(axis) && mid < box.High.Coord(axis)))
        return false;
      const size_t count = children.size();
      for (size_t i = 0; i < count; ++i)
      {
        auto other = children[i];
        children[i].High.SetCoord(axis, mid);
        other.Low.SetCoord(axis, mid);
        children.push_back(other);
      }
    }
    for (const auto& child : children)
      if (!coverBox(child, tree, depth + 1, cap))
        return false;
    return true;
  }

  bool coverTriangle(const gp_Pnt& a, const gp_Pnt& b, const gp_Pnt& c, int depth, int capSide)
  {
    for (int axis = 1; axis <= 3; ++axis)
    {
      Result.UnresolvedMinimum.SetCoord(axis,
                                        std::min({a.Coord(axis), b.Coord(axis), c.Coord(axis)}));
      Result.UnresolvedMaximum.SetCoord(axis,
                                        std::max({a.Coord(axis), b.Coord(axis), c.Coord(axis)}));
    }
    if (Result.TrianglePatches >= MaximumPatches)
      return false;
    ++Result.TrianglePatches;
    const gp_Pnt p((a.XYZ() + b.XYZ() + c.XYZ()) / 3);
    const double r = std::max({p.Distance(a), p.Distance(b), p.Distance(c)}) + 4 * Reserve;
    if (capSide >= 0 && Direction * Field.Value(p) + variation(p, r) < 0)
    {
      Result.MeshToBoundaryBound = std::max(Result.MeshToBoundaryBound, 8 * Reserve);
      return true;
    }
    // A normal-line witness for the whole patch is clipped only when its
    // entire normal slab lies in the box. Never reuse the untrimmed flag.
    const auto   g        = Field.Gradient(p);
    const double lower    = g.Magnitude() - H * r - GradientReserve;
    const double diameter = std::max({a.Distance(b), b.Distance(c), c.Distance(a)});
    const double error =
      std::max({std::abs(Field.Value(a)), std::abs(Field.Value(b)), std::abs(Field.Value(c))})
      + H * diameter * diameter / 6 + FieldReserve;
    const double discriminant = lower * lower - 2 * H * error;
    if (std::isfinite(lower) && std::isfinite(discriminant) && lower > 0 && discriminant > 0)
    {
      const double distance = 2 * error / (lower + std::sqrt(discriminant)) + 8 * Reserve;
      bool         inside   = distance <= Tolerance;
      for (int axis = 1; axis <= 3; ++axis)
      {
        const double reach = distance * std::abs(g.Coord(axis)) / g.Magnitude() + 8 * Reserve;
        inside &=
          std::min({a.Coord(axis), b.Coord(axis), c.Coord(axis)}) - reach >= Low.Coord(axis)
          && std::max({a.Coord(axis), b.Coord(axis), c.Coord(axis)}) + reach <= High.Coord(axis);
      }
      if (inside)
      {
        Result.MeshToBoundaryBound = std::max(Result.MeshToBoundaryBound, distance);
        return true;
      }
    }
    if (r + 8192 * Reserve < Tolerance)
    {
      const double capDistance = capWitness(p, Tolerance - r - 8 * Reserve, capSide);
      if (capDistance + r + 8 * Reserve <= Tolerance)
      {
        Result.MeshToBoundaryBound =
          std::max(Result.MeshToBoundaryBound, capDistance + r + 8 * Reserve);
        return true;
      }
      const auto witness =
        BRepClass3d_ImplicitBox::Classify(Boundary, p, Tolerance - r - 8 * Reserve);
      if (witness.HasBoundaryWitness
          && witness.WitnessDistanceUpperBound + r + 8 * Reserve <= Tolerance)
      {
        Result.MeshToBoundaryBound =
          std::max(Result.MeshToBoundaryBound, witness.WitnessDistanceUpperBound + r + 8 * Reserve);
        return true;
      }
    }
    if (depth >= MaximumDepth)
      return false;
    const auto ab = midpoint(a, b), bc = midpoint(b, c), ca = midpoint(c, a);
    return coverTriangle(a, ab, ca, depth + 1, capSide)
           && coverTriangle(ab, b, bc, depth + 1, capSide)
           && coverTriangle(ca, bc, c, depth + 1, capSide)
           && coverTriangle(ab, bc, ca, depth + 1, capSide);
  }
};
} // namespace

BRepMesh_ImplicitBoundary::Coverage BRepMesh_ImplicitBoundary::CheckCoverage(
  const BRepClass3d_ImplicitBox::Boundary& source,
  const Result&                            mesh,
  double                                   tolerance,
  int                                      maximumDepth,
  size_t                                   maximumCells,
  size_t                                   maximumPatches)
{
  if (source.Policy != BRepClass3d_ImplicitBox::CornerPolicy::Reject)
    throw Standard_NotImplemented("Corner-aware clipped coverage is not implemented");
  if (!std::isfinite(tolerance) || tolerance <= 0 || maximumDepth < 1 || maximumDepth > 20
      || !maximumCells || maximumCells > 10000000 || !maximumPatches || maximumPatches > 4000000
      || mesh.Nodes.size() > 1000000
      || mesh.CapTriangles.size() + mesh.ImplicitTriangles.size() > 2000000)
    throw Standard_DomainError("Invalid clipped coverage settings");
  std::ostringstream record;
  BRepClass3d_ImplicitBoxIO::Write(source, record);
  std::istringstream                input(record.str());
  BRepClass3d_ImplicitBox::Boundary copy;
  BRepClass3d_ImplicitBoxIO::Read(input, copy);
  const auto   face  = occ::down_cast<BRep_TImplicitFace>(copy.Faces.back().TShape());
  const auto&  field = *static_cast<const Geom_ImplicitGyroid*>(face->Field().get());
  const auto   low = face->Minimum(), high = face->Maximum();
  const double k = 2 * std::acos(-1.) / field.Period(), h = 4 * k * k, l = field.LipschitzBound();
  double       magnitude = 1;
  for (int axis = 1; axis <= 3; ++axis)
  {
    magnitude = std::max({magnitude, std::abs(low.Coord(axis)), std::abs(high.Coord(axis))});
    if (!std::isfinite(k * low.Coord(axis)) || !std::isfinite(k * high.Coord(axis))
        || std::abs(k * low.Coord(axis)) >= 1.e6 || std::abs(k * high.Coord(axis)) >= 1.e6)
      throw Standard_DomainError("Coverage trigonometric precision exhausted");
  }
  const double reserve = 4096 * std::numeric_limits<double>::epsilon() * magnitude;
  if (magnitude > 1.e50 || !std::isfinite(h) || !(h > 0) || !std::isfinite(l)
      || tolerance <= 4096 * reserve)
    throw Standard_DomainError("Coverage coordinate/derivative precision exhausted");
  for (const auto& p : mesh.Nodes)
    for (int axis = 1; axis <= 3; ++axis)
      if (!std::isfinite(p.Coord(axis)) || p.Coord(axis) < low.Coord(axis) - reserve
          || p.Coord(axis) > high.Coord(axis) + reserve)
        throw Standard_DomainError("Coverage mesh point outside supported box");
  const auto validate = [&](const std::array<int, 3>& t) {
    for (int n : t)
      if (n < 0 || size_t(n) >= mesh.Nodes.size())
        throw Standard_DomainError("Invalid coverage mesh index");
  };
  for (const auto& t : mesh.ImplicitTriangles)
    validate(t);
  std::vector<int> sides(copy.Faces.size() - 1);
  for (size_t i = 0; i < sides.size(); ++i)
  {
    const auto plane = BRepAdaptor_Surface(copy.Faces[i]).Plane();
    int        axis  = 1;
    for (int j = 2; j <= 3; ++j)
      if (std::abs(plane.Axis().Direction().Coord(j))
          > std::abs(plane.Axis().Direction().Coord(axis)))
        axis = j;
    sides[i] = (axis - 1) * 2
               + (std::abs(plane.Location().Coord(axis) - high.Coord(axis))
                      < std::abs(plane.Location().Coord(axis) - low.Coord(axis))
                    ? 1
                    : 0);
  }
  std::array<std::vector<std::array<int, 3>>, 6> caps;
  for (const auto& t : mesh.CapTriangles)
  {
    validate(t.Nodes);
    if (t.Cap < 0 || size_t(t.Cap) >= sides.size())
      throw Standard_DomainError("Invalid coverage cap index");
    const int    side = sides[size_t(t.Cap)], axis = side / 2 + 1;
    const double position = side % 2 ? high.Coord(axis) : low.Coord(axis);
    for (int n : t.Nodes)
      if (std::abs(mesh.Nodes[size_t(n)].Coord(axis) - position) > reserve)
        throw Standard_DomainError("Coverage cap does not lie on its native plane");
    caps[size_t(side)].push_back(t.Nodes);
  }
  Checker checker{copy,
                  field,
                  low,
                  high,
                  tolerance,
                  reserve,
                  1024 * std::numeric_limits<double>::epsilon() * (1 + std::abs(field.IsoValue()))
                    + l * reserve,
                  h * reserve + 1024 * std::numeric_limits<double>::epsilon() * l,
                  h,
                  l,
                  copy.Faces.back().Orientation() == TopAbs_FORWARD ? 1. : -1.,
                  maximumDepth,
                  maximumCells,
                  maximumPatches,
                  {}};
  int     index = 0;
  for (const auto& t : mesh.CapTriangles)
  {
    checker.Result.UnresolvedTriangle = index++;
    if (!checker.coverTriangle(mesh.Nodes[size_t(t.Nodes[0])],
                               mesh.Nodes[size_t(t.Nodes[1])],
                               mesh.Nodes[size_t(t.Nodes[2])],
                               0,
                               sides[size_t(t.Cap)]))
      return checker.Result;
  }
  for (const auto& t : mesh.ImplicitTriangles)
  {
    checker.Result.UnresolvedTriangle = index++;
    if (!checker.coverTriangle(mesh.Nodes[size_t(t[0])],
                               mesh.Nodes[size_t(t[1])],
                               mesh.Nodes[size_t(t[2])],
                               0,
                               -1))
      return checker.Result;
  }
  checker.Result.UnresolvedTriangle = -1;
  const Tree implicit(mesh, mesh.ImplicitTriangles, reserve);
  if (!checker.coverBox({low, high}, implicit, 0, false))
    return checker.Result;
  for (int side = 0; side < 6; ++side)
  {
    Box       box{low, high};
    const int axis = side / 2 + 1;
    box.Low.SetCoord(axis, side % 2 ? high.Coord(axis) : low.Coord(axis));
    box.High.SetCoord(axis, box.Low.Coord(axis));
    const Tree cap(mesh, std::move(caps[size_t(side)]), reserve);
    if (!checker.coverBox(box, cap, 0, true))
      return checker.Result;
  }
  checker.Result.Resolved = true;
  return checker.Result;
}
