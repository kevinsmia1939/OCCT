// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRep_ImplicitSphereOps.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_ImplicitIO.hxx>
#include <BRep_TImplicitSpherePatch.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Circle.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Iterator.hxx>
#include <Standard_NotImplemented.hxx>
#include <Standard_DomainError.hxx>
#include <gtest/gtest.h>
#include <map>
#include <sstream>

namespace
{
using Ops = BRep_ImplicitSphereOps;

TopoDS_Face sphere(gp_Pnt center = gp_Pnt(), double radius = 2)
{
  TopoDS_Face face;
  BRep_Builder().MakeImplicitFace(face,
                                  new Geom_ImplicitSphere(center, radius),
                                  center.Translated(gp_Vec(-radius - 1, -radius - 1, -radius - 1)),
                                  center.Translated(gp_Vec(radius + 1, radius + 1, radius + 1)),
                                  TopLoc_Location(),
                                  .1);
  return face;
}

TopoDS_Edge boundary(const TopoDS_Face& face)
{
  TopoDS_Iterator wires(face);
  EXPECT_TRUE(wires.More());
  TopoDS_Iterator edges(wires.Value());
  return TopoDS::Edge(edges.Value());
}
} // namespace

TEST(BRep_ImplicitSphereOps, TrimCircleDomainRaysAndCopy)
{
  const auto face  = Ops::Trim(sphere(), gp_Pln(gp_Pnt(0, 0, .5), gp_Dir(0, 0, 1)));
  auto       patch = occ::down_cast<BRep_TImplicitSpherePatch>(face.TShape());
  ASSERT_FALSE(patch.IsNull());
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(face, gp_Pnt(0, 0, -2)));
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(face, gp_Pnt(0, 0, 2)));
  auto hits = BRep_ImplicitTool::IntersectRay(face, gp_Pnt(0, 0, -5), gp_Vec(0, 0, 1), 0, 10);
  ASSERT_EQ(hits.size(), 1u);
  EXPECT_NEAR(hits[0].Point.Z(), -2, .05);
  const auto edge = boundary(face);
  double     first, last;
  const auto curve = BRep_Tool::Curve(edge, first, last);
  ASSERT_FALSE(occ::down_cast<Geom_Circle>(curve).IsNull());
  for (int i = 0; i <= 100; ++i)
  {
    gp_Pnt p = curve->Value(first + (last - first) * i / 100);
    EXPECT_NEAR(p.Z(), .5, 1.e-14);
    EXPECT_NEAR(BRep_ImplicitTool::Value(face, p), 0, 1.e-14);
  }
  auto copy = occ::down_cast<BRep_TImplicitSpherePatch>(patch->EmptyCopy());
  ASSERT_FALSE(copy.IsNull());
  EXPECT_DOUBLE_EQ(copy->Plane().Location().Z(), .5);
  TopoDS_Face empty;
  empty.TShape(copy);
  // Even an EmptyCopy with no wire must not lose its trim in schema 1.
  std::ostringstream out;
  EXPECT_THROW(BRep_ImplicitIO::Write(empty, out), Standard_NotImplemented);
  EXPECT_TRUE(out.str().empty());
  gp_Trsf move;
  move.SetTranslation(gp_Vec(3, 0, 10));
  TopoDS_Face placed = face;
  placed.Location(TopLoc_Location(move));
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(placed, gp_Pnt(3, 0, 12)));
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(placed, gp_Pnt(3, 0, 8)));
  patch->Location(TopLoc_Location(move));
  EXPECT_THROW(BRep_ImplicitTool::IsInDomain(face, gp_Pnt()), Standard_DomainError);
  EXPECT_THROW(Ops::Tessellate(face, .1), Standard_DomainError);
}

TEST(BRep_ImplicitSphereOps, TransverseBooleanBoundaryAndHistory)
{
  const auto a = sphere(), b = sphere(gp_Pnt(2, 0, 0));
  for (auto operation : {Ops::Operation::Common, Ops::Operation::Fuse, Ops::Operation::Cut})
  {
    const auto result = Ops::Boolean(a, b, operation);
    ASSERT_EQ(result.Faces.size(), 2u);
    EXPECT_EQ(result.Sources, (std::vector<int>{0, 1}));
    EXPECT_EQ(result.Boundary.ShapeType(), TopAbs_COMPOUND);
    const auto edgeA = boundary(result.Faces[0]), edgeB = boundary(result.Faces[1]);
    EXPECT_TRUE(edgeA.IsSame(edgeB));
    EXPECT_NE(edgeA.Orientation(), edgeB.Orientation());
    const bool common = operation == Ops::Operation::Common;
    EXPECT_EQ(BRep_ImplicitTool::IsInDomain(result.Faces[0], gp_Pnt(2, 0, 0)), common);
    EXPECT_EQ(BRep_ImplicitTool::IsInDomain(result.Faces[0], gp_Pnt(-2, 0, 0)), !common);
    if (operation == Ops::Operation::Cut)
    {
      EXPECT_EQ(result.Faces[1].Orientation(), TopAbs_REVERSED);
      EXPECT_NEAR(BRep_ImplicitTool::Normal(result.Faces[1], gp_Pnt(0, 0, 0)).X(), 1, 1.e-14);
    }
    else
      EXPECT_EQ(result.Faces[1].Orientation(), TopAbs_FORWARD);
  }
  EXPECT_DOUBLE_EQ(BRep_Tool::Tolerance(a), .1);
  EXPECT_FALSE(TopoDS_Iterator(a).More()); // Operands never mutated.
}

TEST(BRep_ImplicitSphereOps, SeparatedContainedAndUnsupportedCases)
{
  const auto a = sphere(), far = sphere(gp_Pnt(8, 0, 0)), small = sphere(gp_Pnt(.2, 0, 0), .5);
  EXPECT_TRUE(Ops::Boolean(a, far, Ops::Operation::Common).Faces.empty());
  EXPECT_EQ(Ops::Boolean(a, far, Ops::Operation::Fuse).Faces.size(), 2u);
  EXPECT_EQ(Ops::Boolean(a, far, Ops::Operation::Cut).Faces.size(), 1u);
  EXPECT_EQ(Ops::Boolean(a, small, Ops::Operation::Common).Sources, (std::vector<int>{1}));
  EXPECT_EQ(Ops::Boolean(a, small, Ops::Operation::Fuse).Sources, (std::vector<int>{0}));
  auto cavity = Ops::Boolean(a, small, Ops::Operation::Cut);
  ASSERT_EQ(cavity.Faces.size(), 2u);
  EXPECT_EQ(cavity.Faces[1].Orientation(), TopAbs_REVERSED);
  EXPECT_TRUE(Ops::Boolean(small, a, Ops::Operation::Cut).Faces.empty());
  EXPECT_THROW(Ops::Boolean(a, a, Ops::Operation::Common), Standard_DomainError);
  for (double d : {3.95, 4.0, 4.05})
    EXPECT_THROW(Ops::Boolean(a, sphere(gp_Pnt(d, 0, 0)), Ops::Operation::Fuse),
                 Standard_DomainError);
  EXPECT_THROW(Ops::Trim(a, gp_Pln(gp_Pnt(0, 0, 2), gp_Dir(0, 0, 1))), Standard_DomainError);
  auto trimmed = Ops::Trim(a, gp_Pln(gp_Pnt(), gp_Dir(1, 0, 0)));
  EXPECT_THROW(Ops::Trim(trimmed, gp_Pln()), Standard_NotImplemented);
  EXPECT_THROW(Ops::Boolean(a, trimmed, Ops::Operation::Common), Standard_NotImplemented);
  EXPECT_THROW(Ops::Trim(TopoDS_Face(), gp_Pln()), Standard_NotImplemented);
  EXPECT_THROW((void)BRep_TImplicitSpherePatch({}, gp_Pln(), .1), Standard_DomainError);
}

TEST(BRep_ImplicitSphereOps, MeshBoundWindingAndConnectivity)
{
  auto full    = sphere();
  auto trimmed = Ops::Trim(full, gp_Pln(gp_Pnt(0, 0, .7), gp_Dir(0, 0, 1)));
  for (auto face : {full, trimmed})
    for (bool reversed : {false, true})
    {
      if (reversed)
        face.Reverse();
      const auto result = Ops::Tessellate(face, .1);
      EXPECT_LE(result.GeometricBound, .1);
      const auto&                                        mesh = result.Triangulation;
      std::map<std::pair<int, int>, std::pair<int, int>> edges;
      for (int i = 1; i <= mesh->NbNodes(); ++i)
        EXPECT_NEAR(BRep_ImplicitTool::Value(face, mesh->Node(i)), 0, 2.e-14);
      for (int i = 1; i <= mesh->NbTriangles(); ++i)
      {
        int a, b, c;
        mesh->Triangle(i).Get(a, b, c);
        const auto   p = mesh->Node(a), q = mesh->Node(b), r = mesh->Node(c);
        const gp_Vec cross = gp_Vec(p, q).Crossed(gp_Vec(p, r));
        EXPECT_GT(cross.Dot(BRep_ImplicitTool::Normal(face, p)), 0);
        for (auto edge : {std::make_pair(a, b), std::make_pair(b, c), std::make_pair(c, a)})
        {
          int   sign  = edge.first < edge.second ? 1 : -1;
          auto& count = edges[std::minmax(edge.first, edge.second)];
          ++count.first;
          count.second += sign;
        }
        // Dense parameter-independent samples of affine triangle interiors.
        for (int u = 0; u <= 4; ++u)
          for (int v = 0; v <= 4 - u; ++v)
          {
            const gp_Pnt sample(p.XYZ() + (q.XYZ() - p.XYZ()) * (u / 4.0)
                                + (r.XYZ() - p.XYZ()) * (v / 4.0));
            EXPECT_LE(std::abs(BRep_ImplicitTool::Value(face, sample)),
                      result.GeometricBound + 1.e-13);
          }
      }
      int open = 0;
      for (const auto& item : edges)
      {
        EXPECT_LE(item.second.first, 2);
        if (item.second.first == 2)
          EXPECT_EQ(item.second.second, 0);
        else
          ++open;
      }
      EXPECT_EQ(open == 0, face.IsSame(full));
    }
  EXPECT_THROW(Ops::Tessellate(full, 0), Standard_DomainError);
  EXPECT_THROW(Ops::Tessellate(full, 1.e-15), Standard_DomainError);
  EXPECT_THROW(Ops::Tessellate(sphere(gp_Pnt(1.e15, 0, 0), 2), .1), Standard_DomainError);
  EXPECT_THROW(Ops::Tessellate(TopoDS_Face(), .1), Standard_NotImplemented);
}

TEST(BRep_ImplicitSphereOps, UnequalSphereClassificationAndPlacedInputs)
{
  auto a = sphere(gp_Pnt(0, 0, 0), 3), b = sphere(gp_Pnt(2.5, 0, 0), 2);
  for (auto operation : {Ops::Operation::Common, Ops::Operation::Fuse, Ops::Operation::Cut})
  {
    const auto result = Ops::Boolean(a, b, operation);
    ASSERT_EQ(result.Faces.size(), 2u);
    const auto edge = boundary(result.Faces[0]);
    double     first, last;
    const auto curve = BRep_Tool::Curve(edge, first, last);
    for (int i = 0; i < 100; ++i)
    {
      const auto p = curve->Value(first + (last - first) * i / 100);
      EXPECT_NEAR(BRep_ImplicitTool::Value(a, p), 0, 1.e-13);
      EXPECT_NEAR(BRep_ImplicitTool::Value(b, p), 0, 1.e-13);
    }
    for (int source = 0; source < 2; ++source)
    {
      const auto& face  = result.Faces[source];
      const auto  field = occ::down_cast<Geom_ImplicitSphere>(
        occ::down_cast<BRep_TImplicitFace>(face.TShape())->Field());
      for (int i = 0; i < 50; ++i)
      {
        const double angle = i * .137;
        const gp_Pnt p =
          field->Center().Translated(gp_Vec(std::cos(angle), std::sin(angle), 0) * field->Radius());
        const double other = BRep_ImplicitTool::Value(source == 0 ? b : a, p);
        const bool   inside =
          operation == Ops::Operation::Common || (operation == Ops::Operation::Cut && source == 1);
        if (std::abs(other) > 1.e-10)
          EXPECT_EQ(BRep_ImplicitTool::IsInDomain(face, p), inside ? other < 0 : other > 0);
      }
    }
  }
  auto    input = sphere();
  gp_Trsf geometry;
  geometry.SetTranslation(gp_Vec(1, 0, 0));
  occ::down_cast<BRep_TImplicitFace>(input.TShape())->Location(TopLoc_Location(geometry));
  gp_Trsf instance;
  instance.SetScale(gp_Pnt(), -1.5);
  input.Location(TopLoc_Location(instance), false);
  const auto trimmed = Ops::Trim(input, gp_Pln(gp_Pnt(), gp_Dir(0, 0, 1)));
  EXPECT_NEAR(BRep_ImplicitTool::Value(trimmed, gp_Pnt(-1.5, 0, -3)), 0, 1.e-14);
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(trimmed, gp_Pnt(-1.5, 0, 3)));
  auto mesh = Ops::Tessellate(trimmed, .1);
  EXPECT_LE(mesh.GeometricBound, .1);
  const auto* data = static_cast<const BRep_TImplicitFace*>(trimmed.TShape().get());
  EXPECT_TRUE(data->Location().IsIdentity());
  EXPECT_TRUE(trimmed.Location().IsIdentity());
}
