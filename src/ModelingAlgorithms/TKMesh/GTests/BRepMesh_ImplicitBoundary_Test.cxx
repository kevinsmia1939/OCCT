// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepMesh_ImplicitBoundary.hxx>
#include <BRepClass3d_ImplicitBoxIO.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_BSplineCurve.hxx>
#include <BRep_TImplicitFace.hxx>
#include <Poly_Triangulation.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <sstream>

namespace
{
using Ops = BRepClass3d_ImplicitBox;

double referenceMeshDistance(const gp_Pnt&                            p,
                             const BRepMesh_ImplicitBoundary::Result& mesh,
                             bool                                     caps)
{
  double     best     = std::numeric_limits<double>::infinity();
  const auto triangle = [&](const std::array<int, 3>& t) {
    gp_Pnt v[3] = {mesh.Nodes[size_t(t[0])], mesh.Nodes[size_t(t[1])], mesh.Nodes[size_t(t[2])]};
    auto   n    = gp_Vec(v[0], v[1]).Crossed(gp_Vec(v[0], v[2]));
    if (!(n.Magnitude() > 0))
      return;
    n.Normalize();
    const double height = gp_Vec(v[0], p).Dot(n);
    const gp_Pnt q      = p.Translated(n * (-height));
    bool         inside = true;
    for (int i = 0; i < 3; ++i)
    {
      const gp_Vec edge(v[i], v[(i + 1) % 3]);
      inside &= edge.Crossed(gp_Vec(v[i], q)).Dot(n) >= -1.e-12;
      const double u = std::clamp(gp_Vec(v[i], p).Dot(edge) / edge.SquareMagnitude(), 0., 1.);
      best           = std::min(best, p.Distance(v[i].Translated(edge * u)));
    }
    if (inside)
      best = std::min(best, std::abs(height));
  };
  if (caps)
    for (const auto& t : mesh.CapTriangles)
      triangle(t.Nodes);
  else
    for (const auto& t : mesh.ImplicitTriangles)
      triangle(t);
  return best;
}

TEST(BRepMesh_ImplicitBoundary, ExactTriangleAuditRejectsImproperContactsAndHandlesAdjacency)
{
  using Mesher = BRepMesh_ImplicitBoundary;
  using State  = Mesher::MeshAudit::State;
  Mesher::Result mesh;
  mesh.Nodes             = {gp_Pnt(0, 0, 0),
                            gp_Pnt(2, 0, 0),
                            gp_Pnt(0, 2, 0),
                            gp_Pnt(2, 2, 0),
                            gp_Pnt(.5, .5, -1),
                            gp_Pnt(.5, .5, 1),
                            gp_Pnt(3, 3, 0)};
  mesh.ImplicitTriangles = {{0, 1, 2}, {1, 3, 2}};
  EXPECT_EQ(Mesher::Audit(mesh).Status, State::Clear); // shared edge, coplanar
  mesh.ImplicitTriangles = {{0, 1, 2}, {0, 4, 5}};
  EXPECT_EQ(Mesher::Audit(mesh).Status, State::Contact); // shared vertex does NOT excuse crossing
  mesh.ImplicitTriangles = {{0, 1, 2}, {1, 4, 5}};
  EXPECT_EQ(Mesher::Audit(mesh).Status, State::Contact); // crosses beyond common vertex
  mesh.ImplicitTriangles = {{0, 1, 2}, {0, 1, 6}};
  EXPECT_EQ(Mesher::Audit(mesh).Status, State::Contact); // shared edge, same-side overlap
  mesh.ImplicitTriangles = {{0, 1, 2}, {0, 1, 4}};
  EXPECT_EQ(Mesher::Audit(mesh).Status, State::Clear); // shared edge, noncoplanar
  mesh.ImplicitTriangles = {{0, 1, 2}, {2, 1, 0}};
  EXPECT_EQ(Mesher::Audit(mesh).Status, State::Contact); // duplicate face
  mesh.Nodes.push_back(mesh.Nodes[0]);
  mesh.ImplicitTriangles = {{0, 1, 2}, {7, 4, 5}};
  EXPECT_EQ(Mesher::Audit(mesh).Status, State::Contact); // unrelated coincident node
  mesh.Nodes             = {gp_Pnt(0, 0, 0),
                            gp_Pnt(1, 0, 0),
                            gp_Pnt(0, 1, 0),
                            gp_Pnt(0, 0, 1.e-100),
                            gp_Pnt(1, 0, 1.e-100),
                            gp_Pnt(0, 1, 1.e-100)};
  mesh.ImplicitTriangles = {{0, 1, 2}, {3, 4, 5}};
  EXPECT_EQ(Mesher::Audit(mesh).Status, State::Clear); // no epsilon-based contact
  mesh.Nodes[3].SetZ(1.e-300);
  EXPECT_EQ(Mesher::Audit(mesh).Status, State::Limit); // bounded arithmetic, not false clear
  mesh.Nodes[3].SetZ(1.e-100);
  EXPECT_THROW(Mesher::Audit(mesh, 0), Standard_DomainError);
  // Near-parallel planes whose 3D boxes overlap. Exact predicates must retain
  // the tiny gap (and detect a tiny crossing), not rely on the box or epsilon.
  const double tiny = 1.e-100;
  mesh.Nodes        = {gp_Pnt(0, 0, 0),
                       gp_Pnt(1, 0, tiny),
                       gp_Pnt(0, 1, tiny),
                       gp_Pnt(0, 0, tiny),
                       gp_Pnt(1, 0, 2 * tiny),
                       gp_Pnt(0, 1, 2 * tiny)};
  EXPECT_EQ(Mesher::Audit(mesh).Status, State::Clear);
  EXPECT_EQ(Mesher::Audit(mesh).PairTests, 1u);
  mesh.Nodes[3].SetZ(-tiny);
  EXPECT_EQ(Mesher::Audit(mesh).Status, State::Contact);
  // Coplanar triangles have overlapping AABBs but an exact oblique separator.
  // Exercise large/small binary scales and integer carries, not just box culls.
  for (int exponent : {-100, 0, 100})
  {
    const double s         = std::ldexp(1., exponent);
    mesh.Nodes             = {gp_Pnt(0, 0, 0),
                              gp_Pnt(2 * s, 0, 0),
                              gp_Pnt(0, 2 * s, 0),
                              gp_Pnt(.75 * s, 1.5 * s, 0),
                              gp_Pnt(1.5 * s, .75 * s, 0),
                              gp_Pnt(1.5 * s, 1.5 * s, 0)};
    mesh.ImplicitTriangles = {{0, 1, 2}, {3, 4, 5}};
    const auto separated   = Mesher::Audit(mesh);
    EXPECT_EQ(separated.Status, State::Clear);
    EXPECT_EQ(separated.PairTests, 1u);
    mesh.Nodes[3].SetX(.5 * s); // exact touching point on x+y=2*s
    EXPECT_EQ(Mesher::Audit(mesh).Status, State::Contact);
    std::swap(mesh.ImplicitTriangles[0][0], mesh.ImplicitTriangles[0][2]);
    EXPECT_EQ(Mesher::Audit(mesh).Status, State::Contact);
  }
  mesh.Nodes             = {gp_Pnt(0, 0, 0),
                            gp_Pnt(2, 0, 0),
                            gp_Pnt(0, 2, 0),
                            gp_Pnt(0, 0, 2),
                            gp_Pnt(0, -2, -1)};
  mesh.ImplicitTriangles = {{0, 1, 2}, {0, 1, 3}, {0, 1, 4}};
  EXPECT_EQ(Mesher::Audit(mesh).Status, State::Clear);
  EXPECT_EQ(Mesher::Audit(mesh, 1).Status, State::Limit);
  EXPECT_EQ(Mesher::Audit(mesh, 100, 1).Status, State::Limit);
  mesh.ImplicitTriangles[0][1] = 99;
  EXPECT_THROW(Mesher::Audit(mesh), Standard_DomainError);
  mesh.Nodes[4].SetZ(1.e-300);
  EXPECT_THROW(Mesher::Audit(mesh), Standard_DomainError); // index checks precede arithmetic limits
}

TEST(BRepMesh_ImplicitBoundary, SpatialAuditCullsCoplanarSparsePairsWithoutHidingContacts)
{
  using Mesher = BRepMesh_ImplicitBoundary;
  Mesher::Result mesh;
  for (int i = 0; i < 4096; ++i)
  {
    const int n = int(mesh.Nodes.size());
    mesh.Nodes.push_back(gp_Pnt(0, 3 * i, 0));
    mesh.Nodes.push_back(gp_Pnt(0, 3 * i + 1, 0));
    mesh.Nodes.push_back(gp_Pnt(0, 3 * i, 1));
    mesh.ImplicitTriangles.push_back({n, n + 1, n + 2});
  }
  const auto separated = Mesher::Audit(mesh, 1000, 100000);
  EXPECT_EQ(separated.Status, Mesher::MeshAudit::State::Clear);
  EXPECT_EQ(separated.PairTests, 0u);
  EXPECT_LT(separated.BroadTests, 100000u);
  const size_t last = mesh.Nodes.size() - 3;
  for (size_t i = 0; i < 3; ++i) mesh.Nodes[last + i] = mesh.Nodes[i];
  const auto contact = Mesher::Audit(mesh, 1000, 100000);
  EXPECT_EQ(contact.Status, Mesher::MeshAudit::State::Contact);
  EXPECT_EQ(contact.FirstTriangle, 0);
  EXPECT_EQ(contact.SecondTriangle, 4095);
}

TEST(BRepMesh_ImplicitBoundary, NativeCandidateContactAuditDiagnostic)
{
  using Mesher = BRepMesh_ImplicitBoundary;
  TopoDS_Face field;
  BRep_Builder().MakeImplicitFace(field,
                                  new Geom_ImplicitGyroid(10, 0),
                                  gp_Pnt(-5, -5, -5),
                                  gp_Pnt(16, 16, 6),
                                  TopLoc_Location(),
                                  .1);
  const auto box =
    BRepPrimAPI_MakeBox(gp_Pnt(-3.83, -4.07, -2.67), gp_Pnt(4.17, 3.93, 3.37)).Shape();
  const auto source = Ops::Build(field, box, Ops::Operation::Common);
  const auto mesh   = Mesher::Tessellate(source);
  const auto audit  = Mesher::Audit(mesh);
  std::cout << "native contact audit state=" << int(audit.Status) << " pair=" << audit.FirstTriangle
            << "," << audit.SecondTriangle << " tests=" << audit.PairTests << std::endl;
  EXPECT_EQ(audit.Status, Mesher::MeshAudit::State::Clear);
}

Ops::Boundary boundary(bool crossing, Ops::Operation operation)
{
  TopoDS_Face field;
  BRep_Builder().MakeImplicitFace(field,
                                  new Geom_ImplicitGyroid(10, crossing ? 0 : 1.45),
                                  gp_Pnt(-5, -5, -5),
                                  gp_Pnt(16, 16, 6),
                                  TopLoc_Location(),
                                  .1);
  const auto box =
    crossing ? BRepPrimAPI_MakeBox(gp_Pnt(-3.83, -4.07, -2.67), gp_Pnt(4.17, 3.93, 3.37)).Shape()
             : BRepPrimAPI_MakeBox(gp_Pnt(-.31, -.37, 1.25), gp_Pnt(2.93, 2.91, 4)).Shape();
  return Ops::Build(field, box, operation);
}

struct Use
{
  int Count = 0, Balance = 0;
};

using Key = std::pair<int, int>;

void add(std::map<Key, Use>& uses, int a, int b)
{
  ASSERT_NE(a, b);
  auto& value = uses[std::minmax(a, b)];
  ++value.Count;
  value.Balance += a < b ? 1 : -1;
}
} // namespace

Ops::Boundary cornerBoundary(Ops::Operation operation)
{
  TopoDS_Face field;
  BRep_Builder().MakeImplicitFace(field, new Geom_ImplicitGyroid(10, 0), gp_Pnt(),
                                  gp_Pnt(10, 10, 10), TopLoc_Location(), .1);
  return Ops::Build(field, BRepPrimAPI_MakeBox(10, 10, 10).Shape(), operation,
                    .1, 512, Ops::CornerPolicy::Regularized);
}

TEST(BRepMesh_ImplicitBoundary, CornerCapConstraintsShareOnlyNativeCornerIdentities)
{
  using Mesher = BRepMesh_ImplicitBoundary;
  for (const auto operation : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
  {
    const auto source = cornerBoundary(operation);
    std::ostringstream before;
    BRepClass3d_ImplicitBoxIO::Write(source, before);
    const auto prepared = Mesher::Prepare(source);
    EXPECT_EQ(prepared.Policy, Ops::CornerPolicy::Regularized);
    EXPECT_EQ(prepared.FixedCornerCount, 8);
    ASSERT_EQ(prepared.Wires.size(), 1u);
    EXPECT_EQ(prepared.Seams.size(), 18u);
    std::set<int> corners;
    for (int node : prepared.Wires[0])
    {
      const auto p = prepared.Nodes[size_t(node)];
      if ((p.X() == 0 || p.X() == 10) && (p.Y() == 0 || p.Y() == 10)
          && (p.Z() == 0 || p.Z() == 10)) corners.insert(node);
    }
    EXPECT_EQ(corners.size(), 6u);
    EXPECT_FALSE(prepared.NativeSurfaceMeshed || prepared.SolidValidated);
    for (const auto& face : source.Faces)
    {
      TopLoc_Location location;
      EXPECT_TRUE(BRep_Tool::Triangulation(face, location).IsNull());
    }
    std::ostringstream after;
    BRepClass3d_ImplicitBoxIO::Write(source, after);
    EXPECT_EQ(before.str(), after.str());
  }
}

TEST(BRepMesh_ImplicitBoundary, CornerFieldCandidateRetainsTrimCornersAndExcludesContacts)
{
  using Mesher = BRepMesh_ImplicitBoundary;
  for (const auto operation : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
  {
    const auto source = cornerBoundary(operation);
    const auto mesh = Mesher::Tessellate(source);
    EXPECT_TRUE(mesh.NativeSurfaceMeshed && mesh.MeshClosed);
    EXPECT_FALSE(mesh.SurfaceCoverageVerified || mesh.SolidValidated);
    EXPECT_EQ(mesh.Policy, Ops::CornerPolicy::Regularized);
    EXPECT_EQ(mesh.FixedCornerCount, 8);
    EXPECT_EQ(mesh.Cells % 2, 1);
    EXPECT_LE(mesh.FieldDistanceBound, .1);
    std::set<int> used;
    for (const auto& t : mesh.ImplicitTriangles) for (int n : t) used.insert(n);
    std::set<int> corners;
    for (int n : used)
    {
      const auto p = mesh.Nodes[size_t(n)];
      EXPECT_GT(p.Distance(gp_Pnt()), .1);
      EXPECT_GT(p.Distance(gp_Pnt(10, 10, 10)), .1);
      if ((p.X() == 0 || p.X() == 10) && (p.Y() == 0 || p.Y() == 10)
          && (p.Z() == 0 || p.Z() == 10)) corners.insert(n);
    }
    EXPECT_EQ(corners.size(), 6u);
    const auto audit = Mesher::Audit(mesh);
    EXPECT_EQ(audit.Status, Mesher::MeshAudit::State::Clear)
      << "broad=" << audit.BroadTests << " exact=" << audit.PairTests
      << " caps=" << mesh.CapTriangles.size() << " field=" << mesh.ImplicitTriangles.size();
  }
}

TEST(BRepMesh_ImplicitBoundary, CapsAndNativeWiresShareCanonicalSamplesAndOppositeIncidences)
{
  for (bool crossing : {false, true})
    for (auto operation : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
    {
      const auto source   = boundary(crossing, operation);
      const auto prepared = BRepMesh_ImplicitBoundary::Prepare(source);
      ASSERT_FALSE(prepared.CapTriangles.empty());
      ASSERT_EQ(prepared.Seams.size(), size_t(source.SeamCount));
      ASSERT_FALSE(prepared.Wires.empty());
      EXPECT_FALSE(prepared.NativeSurfaceMeshed);
      EXPECT_FALSE(prepared.SolidValidated);
      EXPECT_TRUE(prepared.SeamBoundVerified);
      EXPECT_DOUBLE_EQ(prepared.SeamGeometricBound, source.GeometricBound);
      std::map<Key, Use> caps, wires;
      for (const auto& triangle : prepared.CapTriangles)
      {
        const auto& n = triangle.Nodes;
        for (int index : n)
          ASSERT_LT(size_t(index), prepared.Nodes.size());
        const auto normal =
          BRepAdaptor_Surface(source.Faces.at(size_t(triangle.Cap))).Plane().Axis().Direction();
        const auto area =
          gp_Vec(prepared.Nodes[size_t(n[0])], prepared.Nodes[size_t(n[1])])
            .Crossed(gp_Vec(prepared.Nodes[size_t(n[0])], prepared.Nodes[size_t(n[2])]));
        EXPECT_GT(area.Dot(gp_Vec(normal)), 0);
        add(caps, n[0], n[1]);
        add(caps, n[1], n[2]);
        add(caps, n[2], n[0]);
      }
      for (const auto& wire : prepared.Wires)
      {
        ASSERT_GE(wire.size(), 4u);
        EXPECT_EQ(wire.front(), wire.back());
        for (size_t i = 1; i < wire.size(); ++i)
          add(wires, wire[i - 1], wire[i]);
      }
      for (const auto& entry : wires)
      {
        ASSERT_EQ(caps.count(entry.first), 1u);
        EXPECT_EQ(entry.second.Count, 1);
        EXPECT_EQ(caps.at(entry.first).Count, 1);
        EXPECT_EQ(caps.at(entry.first).Balance + entry.second.Balance, 0);
      }
      for (const auto& entry : caps)
        if (!wires.count(entry.first))
        {
          EXPECT_EQ(entry.second.Count, 2);
          EXPECT_EQ(entry.second.Balance, 0);
        }
      NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> edges;
      TopExp::MapShapes(source.Faces.back(), TopAbs_EDGE, edges);
      for (int index = 1; index <= edges.Extent(); ++index)
      {
        TopLoc_Location location;
        double          first, last;
        const auto      curve = occ::down_cast<Geom_BSplineCurve>(
          BRep_Tool::Curve(TopoDS::Edge(edges.FindKey(index)), location, first, last));
        const auto& samples = prepared.Seams[size_t(index - 1)];
        ASSERT_GE(samples.size(), size_t(curve->NbPoles()));
        for (int pole = 1; pole <= curve->NbPoles(); ++pole)
        {
          double distance = std::numeric_limits<double>::infinity();
          for (int node : samples)
            distance =
              std::min(distance, prepared.Nodes.at(size_t(node)).Distance(curve->Pole(pole)));
          EXPECT_LT(distance, 1.e-9);
        }
      }
      TopLoc_Location location;
      EXPECT_TRUE(BRep_Tool::Triangulation(source.Faces.back(), location).IsNull());
    }
}

TEST(BRepMesh_ImplicitBoundary, PrivatePreparationPreservesCachesAndRejectsBudgets)
{
  auto source = boundary(false, Ops::Operation::Common);
  auto mesh   = occ::handle<Poly_Triangulation>(new Poly_Triangulation(3, 1, false));
  mesh->SetNode(1, gp_Pnt(-.31, -.37, 1.25));
  mesh->SetNode(2, gp_Pnt(-.31, 2.91, 1.25));
  mesh->SetNode(3, gp_Pnt(-.31, -.37, 4));
  mesh->SetTriangle(1, Poly_Triangle(1, 2, 3));
  BRep_Builder().UpdateFace(source.Faces.front(), mesh);
  const bool free = source.Shape.Free();
  EXPECT_FALSE(BRepMesh_ImplicitBoundary::Prepare(source).CapTriangles.empty());
  EXPECT_THROW(BRepMesh_ImplicitBoundary::Prepare(source, .1, 3), Standard_NotImplemented);
  EXPECT_THROW(BRepMesh_ImplicitBoundary::Prepare(source, .1, 200000, 1), Standard_NotImplemented);
  EXPECT_THROW(BRepMesh_ImplicitBoundary::Prepare(source, 0), Standard_DomainError);
  EXPECT_THROW(BRepMesh_ImplicitBoundary::Prepare(source, .1, 1000001), Standard_DomainError);
  TopLoc_Location location;
  EXPECT_EQ(BRep_Tool::Triangulation(source.Faces.front(), location), mesh);
  EXPECT_TRUE(BRep_Tool::Triangulation(source.Faces.back(), location).IsNull());
  EXPECT_EQ(source.Shape.Free(), free);
  std::ostringstream record;
  BRepClass3d_ImplicitBoxIO::Write(source, record);
  std::istringstream input(record.str());
  Ops::Boundary      restored;
  BRepClass3d_ImplicitBoxIO::Read(input, restored);
  const auto prepared = BRepMesh_ImplicitBoundary::Prepare(restored);
  EXPECT_FALSE(prepared.SeamBoundVerified);
  EXPECT_DOUBLE_EQ(prepared.SeamGeometricBound, source.GeometricBound);
}

TEST(BRepMesh_ImplicitBoundary, TrimConformingSurfaceCandidateBothOperationsAndOpenContours)
{
  for (bool crossing : {false, true})
    for (auto operation : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
    {
      const auto source   = boundary(crossing, operation);
      const auto mesh     = BRepMesh_ImplicitBoundary::Tessellate(source);
      const auto prepared = BRepMesh_ImplicitBoundary::Prepare(source);
      ASSERT_GE(mesh.Nodes.size(), prepared.Nodes.size());
      EXPECT_EQ(mesh.Seams, prepared.Seams);
      EXPECT_EQ(mesh.Wires, prepared.Wires);
      for (size_t i = 0; i < prepared.Nodes.size(); ++i)
        EXPECT_DOUBLE_EQ(mesh.Nodes[i].Distance(prepared.Nodes[i]), 0);
      EXPECT_TRUE(mesh.NativeSurfaceMeshed);
      EXPECT_TRUE(mesh.MeshClosed);
      EXPECT_FALSE(mesh.SurfaceCoverageVerified);
      EXPECT_FALSE(mesh.SolidValidated);
      EXPECT_GT(mesh.FieldDistanceBound, 0);
      EXPECT_LE(mesh.FieldDistanceBound, .1);
      ASSERT_FALSE(mesh.ImplicitTriangles.empty());
      const auto field = occ::down_cast<BRep_TImplicitFace>(source.Faces.back().TShape())->Field();
      const double       direction = operation == Ops::Operation::Common ? 1 : -1;
      std::map<Key, Use> uses;
      for (const auto& t : mesh.ImplicitTriangles)
      {
        for (int i = 0; i < 3; ++i)
        {
          ASSERT_GE(t[i], 0);
          ASSERT_LT(size_t(t[i]), mesh.Nodes.size());
          add(uses, t[i], t[(i + 1) % 3]);
        }
        const auto&  a = mesh.Nodes[size_t(t[0])];
        const auto&  b = mesh.Nodes[size_t(t[1])];
        const auto&  c = mesh.Nodes[size_t(t[2])];
        const gp_Pnt center((a.XYZ() + b.XYZ() + c.XYZ()) / 3);
        EXPECT_GT(direction * gp_Vec(a, b).Crossed(gp_Vec(a, c)).Dot(field->Gradient(center)), 0);
        // Independent numerical witnesses, not the production bound formula.
        for (const auto& p : {a, center, gp_Pnt((a.XYZ() + b.XYZ()) * .5)})
        {
          auto normal = field->Gradient(p);
          ASSERT_GT(normal.Magnitude(), 0);
          normal.Normalize();
          EXPECT_LE(field->Value(p.Translated(normal * (-mesh.FieldDistanceBound))), 1.e-10);
          EXPECT_GE(field->Value(p.Translated(normal * mesh.FieldDistanceBound)), -1.e-10);
        }
      }
      for (const auto& t : mesh.CapTriangles)
        for (int i = 0; i < 3; ++i)
          add(uses, t.Nodes[i], t.Nodes[(i + 1) % 3]);
      for (const auto& use : uses)
      {
        EXPECT_EQ(use.second.Count, 2);
        EXPECT_EQ(use.second.Balance, 0);
      }
    }
}

TEST(BRepMesh_ImplicitBoundary, SurfaceLimitsAtomicityAndRestoredProvenance)
{
  const auto         source = boundary(true, Ops::Operation::Common);
  std::ostringstream before;
  BRepClass3d_ImplicitBoxIO::Write(source, before);
  EXPECT_THROW(BRepMesh_ImplicitBoundary::Tessellate(source, 0), Standard_DomainError);
  EXPECT_THROW(BRepMesh_ImplicitBoundary::Tessellate(source, .1, 257), Standard_DomainError);
  EXPECT_THROW(BRepMesh_ImplicitBoundary::Tessellate(source, .1, 7), Standard_DomainError);
  EXPECT_THROW(BRepMesh_ImplicitBoundary::Tessellate(source, .1, 8), Standard_NotImplemented);
  EXPECT_THROW(BRepMesh_ImplicitBoundary::Tessellate(source, .1, 128, 3), Standard_NotImplemented);
  EXPECT_THROW(BRepMesh_ImplicitBoundary::Tessellate(source, .1, 128, 200000, 1),
               Standard_NotImplemented);
  Ops::Boundary      restored;
  std::istringstream input(before.str());
  BRepClass3d_ImplicitBoxIO::Read(input, restored);
  const auto mesh = BRepMesh_ImplicitBoundary::Tessellate(restored);
  EXPECT_FALSE(mesh.SeamBoundVerified);
  EXPECT_FALSE(mesh.SurfaceCoverageVerified);
  EXPECT_FALSE(mesh.SolidValidated);
  TopLoc_Location location;
  for (const auto& face : source.Faces)
    EXPECT_TRUE(BRep_Tool::Triangulation(face, location).IsNull());
  std::ostringstream after;
  BRepClass3d_ImplicitBoxIO::Write(source, after);
  EXPECT_EQ(before.str(), after.str());
}

TEST(BRepMesh_ImplicitBoundary, SurfaceMultipleLoopsAndTranslatedPhysicalScale)
{
  TopoDS_Face field;
  BRep_Builder().MakeImplicitFace(field,
                                  new Geom_ImplicitGyroid(10, 1.45),
                                  gp_Pnt(-5, -5, -5),
                                  gp_Pnt(16, 16, 6),
                                  TopLoc_Location(),
                                  .1);
  const auto box = BRepPrimAPI_MakeBox(gp_Pnt(-.31, -.37, 1.25), gp_Pnt(12.93, 2.91, 4)).Shape();
  for (auto operation : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
  {
    const auto source = Ops::Build(field, box, operation);
    const auto mesh   = BRepMesh_ImplicitBoundary::Tessellate(source);
    EXPECT_EQ(mesh.Wires.size(), 2u);
    EXPECT_TRUE(mesh.MeshClosed);
    EXPECT_LE(mesh.FieldDistanceBound, .1);
  }
  // Periodic translations retain the same phase, without placing the native
  // face or relaxing the identity-operand Boolean contract.
  for (double scale : {.5, 2.})
  {
    const gp_Vec shift(20 * scale, -10 * scale, 30 * scale);
    const auto   placed = [&](double x, double y, double z) {
      return gp_Pnt(x * scale, y * scale, z * scale).Translated(shift);
    };
    TopoDS_Face scaled;
    BRep_Builder().MakeImplicitFace(scaled,
                                    new Geom_ImplicitGyroid(10 * scale, 1.45),
                                    placed(-5, -5, -5),
                                    placed(16, 16, 6),
                                    TopLoc_Location(),
                                    .1);
    const auto solid = BRepPrimAPI_MakeBox(placed(-.31, -.37, 1.25), placed(2.93, 2.91, 4)).Shape();
    const auto source = Ops::Build(scaled, solid, Ops::Operation::Common);
    const auto mesh   = BRepMesh_ImplicitBoundary::Tessellate(source);
    EXPECT_TRUE(mesh.MeshClosed);
    EXPECT_LE(mesh.FieldDistanceBound, .1);
  }
}

TEST(BRepMesh_ImplicitBoundary, ClippedCoverageBothDirectionsAndBooleanOperations)
{
  for (bool crossing : {false, true})
    for (auto operation : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
    {
      const auto source = boundary(crossing, operation);
      const auto mesh   = BRepMesh_ImplicitBoundary::Tessellate(source);
      const auto report = BRepMesh_ImplicitBoundary::CheckCoverage(source, mesh);
      std::cout << "coverage " << crossing << " " << int(operation)
                << " resolved=" << report.Resolved << " bounds=" << report.BoundaryToMeshBound
                << "," << report.MeshToBoundaryBound << " cells=" << report.SpatialCells
                << " patches=" << report.TrianglePatches
                << " triangle=" << report.UnresolvedTriangle
                << " box=" << report.UnresolvedMinimum.X() << "," << report.UnresolvedMinimum.Y()
                << "," << report.UnresolvedMinimum.Z() << std::endl;
      EXPECT_TRUE(report.Resolved);
      EXPECT_GT(report.BoundaryToMeshBound, 0);
      EXPECT_LE(report.BoundaryToMeshBound, .1);
      EXPECT_LE(report.MeshToBoundaryBound, .1);
      // Independent harmonic z roots and raw-equation cap samples, not the
      // production octree, BVH or classifier witness implementation.
      const auto   data = occ::down_cast<BRep_TImplicitFace>(source.Faces.back().TShape());
      const auto   low = data->Minimum(), high = data->Maximum();
      const double k = 2 * std::acos(-1.) / 10, iso = crossing ? 0 : 1.45;
      int          roots = 0;
      for (int i = 1; i < 6; ++i)
        for (int j = 1; j < 6; ++j)
        {
          const double x = low.X() + (high.X() - low.X()) * (double(i) / 6);
          const double y = low.Y() + (high.Y() - low.Y()) * (double(j) / 6);
          const double a = std::sin(k * y), b = std::cos(k * x),
                       d = iso - std::sin(k * x) * std::cos(k * y);
          if (std::abs(d) >= std::hypot(a, b))
            continue;
          const double angle = std::asin(d / std::hypot(a, b)), phase = std::atan2(a, b);
          for (double branch : {angle, std::acos(-1.) - angle})
            for (int turn = -2; turn <= 2; ++turn)
            {
              const double z = (branch - phase + 2 * std::acos(-1.) * turn) / k;
              if (z < low.Z() || z > high.Z())
                continue;
              ++roots;
              EXPECT_LE(referenceMeshDistance(gp_Pnt(x, y, z), mesh, false),
                        report.BoundaryToMeshBound + 1.e-8);
            }
        }
      EXPECT_GT(roots, 0);
      for (int side = 0; side < 6; ++side)
        for (int i = 1; i < 4; ++i)
          for (int j = 1; j < 4; ++j)
          {
            const int axis = side / 2 + 1, u = axis % 3 + 1, v = u % 3 + 1;
            gp_Pnt    p;
            p.SetCoord(axis, side % 2 ? high.Coord(axis) : low.Coord(axis));
            p.SetCoord(u, low.Coord(u) + (high.Coord(u) - low.Coord(u)) * (double(i) / 4));
            p.SetCoord(v, low.Coord(v) + (high.Coord(v) - low.Coord(v)) * (double(j) / 4));
            const double value = std::sin(k * p.X()) * std::cos(k * p.Y())
                                 + std::sin(k * p.Y()) * std::cos(k * p.Z())
                                 + std::sin(k * p.Z()) * std::cos(k * p.X()) - iso;
            if ((operation == Ops::Operation::Common ? value : -value) < -1.e-8)
              EXPECT_LE(referenceMeshDistance(p, mesh, true), report.BoundaryToMeshBound + 1.e-8);
          }
    }
}

TEST(BRepMesh_ImplicitBoundary, CoverageRefusesMissingGeometryAndExhaustedBudgets)
{
  const auto source = boundary(false, Ops::Operation::Common);
  auto       mesh   = BRepMesh_ImplicitBoundary::Tessellate(source);
  EXPECT_FALSE(BRepMesh_ImplicitBoundary::CheckCoverage(source, mesh, .1, 1).Resolved);
  EXPECT_FALSE(BRepMesh_ImplicitBoundary::CheckCoverage(source, mesh, .1, 14, 1).Resolved);
  EXPECT_FALSE(BRepMesh_ImplicitBoundary::CheckCoverage(source, mesh, .1, 14, 2000000, 1).Resolved);
  EXPECT_THROW(BRepMesh_ImplicitBoundary::CheckCoverage(source, mesh, 0), Standard_DomainError);
  EXPECT_THROW(BRepMesh_ImplicitBoundary::CheckCoverage(source, mesh, .1, 21),
               Standard_DomainError);
  mesh.SurfaceCoverageVerified = true; // cannot bypass geometry by forging metadata
  mesh.ImplicitTriangles.clear();
  EXPECT_FALSE(BRepMesh_ImplicitBoundary::CheckCoverage(source, mesh).Resolved);
  mesh = BRepMesh_ImplicitBoundary::Tessellate(source);
  mesh.CapTriangles.clear();
  EXPECT_FALSE(BRepMesh_ImplicitBoundary::CheckCoverage(source, mesh).Resolved);
}

TEST(BRepMesh_ImplicitBoundary, CoverageRejectsOneMissingComponentAndExtraneousInteriorTriangles)
{
  TopoDS_Face field;
  BRep_Builder().MakeImplicitFace(field,
                                  new Geom_ImplicitGyroid(10, 1.45),
                                  gp_Pnt(-5, -5, -5),
                                  gp_Pnt(16, 16, 6),
                                  TopLoc_Location(),
                                  .1);
  const auto box    = BRepPrimAPI_MakeBox(gp_Pnt(-.31, -.37, 1.25), gp_Pnt(12.93, 2.91, 4)).Shape();
  const auto source = Ops::Build(field, box, Ops::Operation::Common);
  auto       mesh   = BRepMesh_ImplicitBoundary::Tessellate(source);
  ASSERT_TRUE(BRepMesh_ImplicitBoundary::CheckCoverage(source, mesh).Resolved);
  const auto count = mesh.ImplicitTriangles.size();
  mesh.ImplicitTriangles.erase(
    std::remove_if(mesh.ImplicitTriangles.begin(),
                   mesh.ImplicitTriangles.end(),
                   [&](const auto& t) { return mesh.Nodes[size_t(t[0])].X() > 5; }),
    mesh.ImplicitTriangles.end());
  ASSERT_FALSE(mesh.ImplicitTriangles.empty());
  ASSERT_LT(mesh.ImplicitTriangles.size(), count);
  EXPECT_FALSE(BRepMesh_ImplicitBoundary::CheckCoverage(source, mesh).Resolved);

  const auto small = boundary(false, Ops::Operation::Common);
  mesh             = BRepMesh_ImplicitBoundary::Tessellate(small);
  const int index  = int(mesh.Nodes.size());
  mesh.Nodes.insert(mesh.Nodes.end(), {gp_Pnt(0, 0, 2), gp_Pnt(.01, 0, 2), gp_Pnt(0, .01, 2)});
  mesh.ImplicitTriangles.push_back({index, index + 1, index + 2});
  const auto bad = BRepMesh_ImplicitBoundary::CheckCoverage(small, mesh);
  EXPECT_FALSE(bad.Resolved);
  EXPECT_EQ(bad.UnresolvedTriangle,
            int(mesh.CapTriangles.size() + mesh.ImplicitTriangles.size() - 1));
  mesh.ImplicitTriangles.back()[0] = -1;
  EXPECT_THROW(BRepMesh_ImplicitBoundary::CheckCoverage(small, mesh), Standard_DomainError);
}
