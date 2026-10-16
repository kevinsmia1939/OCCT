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
