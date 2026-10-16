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
#include <sstream>

namespace
{
using Ops = BRepClass3d_ImplicitBox;

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
