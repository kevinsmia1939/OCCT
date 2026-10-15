// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BinTools.hxx>
#include <Geom_ImplicitSphere.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_Plane.hxx>
#include <Poly_Triangulation.hxx>
#include <Poly_Triangle.hxx>
#include <Standard_NotImplemented.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_TypeMismatch.hxx>
#include <TopoDS.hxx>
#include <gp_Ax1.hxx>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <sstream>

namespace
{
class PlaneField : public Geom_ImplicitField
{
public:
  double Value(const gp_Pnt& p) const override { return p.Z(); }

  gp_Vec Gradient(const gp_Pnt&) const override { return gp_Vec(0, 0, 1); }
};

TopoDS_Face sphereFace(const gp_Pnt& low = gp_Pnt(-3, -3, -3), const gp_Pnt& high = gp_Pnt(3, 3, 3))
{
  TopoDS_Face face;
  BRep_Builder().MakeImplicitFace(face,
                                  new Geom_ImplicitSphere(gp_Pnt(0, 0, 0), 2),
                                  low,
                                  high,
                                  TopLoc_Location(),
                                  0.1);
  return face;
}
} // namespace

TEST(BRep_ImplicitFace, NativeRepresentationAndLegacyRejection)
{
  TopoDS_Face f = sphereFace();
  ASSERT_EQ(f.ShapeType(), TopAbs_FACE);
  ASSERT_TRUE(BRep_Tool::IsImplicit(f));
  EXPECT_TRUE(BRep_Tool::IsGeometric(f));
  EXPECT_FALSE(BRep_Tool::HasParametricSurface(f));
  TopLoc_Location l;
  auto            field = BRep_Tool::ImplicitGeometry(f, l);
  EXPECT_TRUE(field->IsKind(STANDARD_TYPE(Geom_ImplicitSphere)));
  EXPECT_DOUBLE_EQ(BRep_Tool::Tolerance(f), 0.1);
  EXPECT_FALSE(BRep_Tool::NaturalRestriction(f));
  EXPECT_TRUE(BRep_Tool::Triangulation(f, l).IsNull());
  EXPECT_THROW(BRep_Tool::Surface(f), Standard_NotImplemented);
  EXPECT_THROW(BRep_Tool::Surface(f, l), Standard_NotImplemented);
  EXPECT_THROW((void)BRepAdaptor_Surface(f), Standard_NotImplemented);
  const auto* data = static_cast<const BRep_TFace*>(f.TShape().get());
  EXPECT_THROW(data->Surface(), Standard_NotImplemented);
  TopoDS_Face nullFace;
  EXPECT_FALSE(BRep_Tool::IsImplicit(nullFace));
}

TEST(BRep_ImplicitFace, ConventionalFacesRemainSupported)
{
  TopoDS_Face  f;
  BRep_Builder builder;
  builder.MakeFace(f, new Geom_Plane(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 0.1);
  EXPECT_FALSE(BRep_Tool::IsImplicit(f));
  EXPECT_TRUE(BRep_Tool::IsGeometric(f));
  EXPECT_TRUE(BRep_Tool::HasParametricSurface(f));
  EXPECT_FALSE(BRep_Tool::Surface(f).IsNull());
  EXPECT_NO_THROW((void)BRepAdaptor_Surface(f));
  TopLoc_Location l;
  EXPECT_THROW(BRep_Tool::ImplicitGeometry(f, l), Standard_TypeMismatch);
  EXPECT_THROW(BRep_ImplicitTool::Value(f, gp_Pnt()), Standard_TypeMismatch);
  std::ostringstream ascii, binary;
  EXPECT_NO_THROW(BRepTools::Write(f, ascii));
  EXPECT_NO_THROW(BinTools::Write(f, binary));
}

TEST(BRep_ImplicitFace, FieldsAndAnalyticalGradients)
{
  Geom_ImplicitSphere sphere(gp_Pnt(), 2);
  EXPECT_DOUBLE_EQ(sphere.Value(gp_Pnt(0, 0, 0)), -2);
  EXPECT_DOUBLE_EQ(sphere.Value(gp_Pnt(3, 0, 0)), 1);
  EXPECT_DOUBLE_EQ(sphere.Gradient(gp_Pnt(3, 0, 0)).X(), 1);
  Geom_ImplicitGyroid gyro(10, 0.17);
  EXPECT_NEAR(gyro.Value(gp_Pnt()), -0.17, 1e-15);
  gp_Pnt p(1.2, -0.7, 0.3);
  gp_Vec g = gyro.Gradient(p);
  for (int axis = 1; axis <= 3; ++axis)
  {
    gp_Pnt lo = p, hi = p;
    lo.SetCoord(axis, p.Coord(axis) - 1e-5);
    hi.SetCoord(axis, p.Coord(axis) + 1e-5);
    EXPECT_NEAR(g.Coord(axis), (gyro.Value(hi) - gyro.Value(lo)) / 2e-5, 1e-9);
  }
  EXPECT_THROW(Geom_ImplicitSphere(gp_Pnt(), 0), Standard_DomainError);
  EXPECT_THROW(Geom_ImplicitGyroid(0), Standard_DomainError);
  EXPECT_THROW(gyro.Value(gp_Pnt(std::numeric_limits<double>::infinity(), 0, 0)),
               Standard_DomainError);
}

TEST(BRep_ImplicitFace, CombinedLocationsAndOrientation)
{
  BRep_Builder builder;
  gp_Trsf      translation;
  translation.SetTranslation(gp_Vec(10, 0, 0));
  TopoDS_Face f;
  builder.MakeImplicitFace(f,
                           new Geom_ImplicitSphere(gp_Pnt(), 2),
                           gp_Pnt(-3, -3, -3),
                           gp_Pnt(3, 3, 3),
                           TopLoc_Location(translation),
                           0.1);
  gp_Trsf rotation;
  rotation.SetRotation(gp_Ax1(gp_Pnt(), gp_Dir(0, 0, 1)), std::acos(-1.0) / 2);
  f.Move(TopLoc_Location(rotation));
  gp_Pnt p(0, 12.5, 0);
  EXPECT_NEAR(BRep_ImplicitTool::Value(f, p), 0.5, 1e-12);
  auto g = BRep_ImplicitTool::Gradient(f, p);
  EXPECT_NEAR(g.X(), 0, 1e-12);
  EXPECT_NEAR(g.Y(), 1, 1e-12);
  EXPECT_NEAR(BRep_ImplicitTool::Normal(f, p).Y(), 1, 1e-12);
  f.Reverse();
  EXPECT_NEAR(BRep_ImplicitTool::Normal(f, p).Y(), -1, 1e-12);
  EXPECT_NEAR(BRep_ImplicitTool::Value(f, p), 0.5, 1e-12);
  auto   box = BRep_ImplicitTool::Bounds(f);
  double x0, y0, z0, x1, y1, z1;
  box.Get(x0, y0, z0, x1, y1, z1);
  EXPECT_NEAR(y0, 7, 1e-9);
  EXPECT_NEAR(y1, 13, 1e-9);
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(f, p));
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(f, gp_Pnt(0, 20, 0)));
  auto hits = BRep_ImplicitTool::IntersectRay(f, gp_Pnt(0, 5, 0.31), gp_Vec(0, 3, 0), 0, 10, 13);
  ASSERT_EQ(hits.size(), 2u);
  EXPECT_NEAR(hits.front().Distance, 5 - std::sqrt(4 - 0.31 * 0.31), 0.05);
}

TEST(BRep_ImplicitFace, GeometricRayTolerancePointOne)
{
  TopoDS_Face  f = sphereFace();
  const gp_Pnt origin(-5, 0.31, 0.27);
  const double offset = std::sqrt(4 - 0.31 * 0.31 - 0.27 * 0.27);
  auto         hits   = BRep_ImplicitTool::IntersectRay(f, origin, gp_Vec(4, 0, 0), 0, 10, 11);
  ASSERT_EQ(hits.size(), 2u);
  for (std::size_t i = 0; i < hits.size(); ++i)
  {
    const double exact = i == 0 ? 5 - offset : 5 + offset;
    EXPECT_LE(hits[i].BracketWidth, 0.1);
    EXPECT_LE(std::abs(hits[i].Distance - exact), 0.05);
    EXPECT_DOUBLE_EQ(hits[i].Residual, std::abs(BRep_ImplicitTool::Value(f, hits[i].Point)));
  }
  auto precise = BRep_ImplicitTool::IntersectRay(f, origin, gp_Vec(1, 0, 0), 0, 10, 11, 1e-6);
  ASSERT_EQ(precise.size(), 2u);
  EXPECT_LE(std::abs(precise.front().Distance - (5 - offset)), 1e-6);
}

TEST(BRep_ImplicitFace, DomainDoesNotCreateCapsAndQueriesFailExplicitly)
{
  auto f = sphereFace(gp_Pnt(-1, -1, -1), gp_Pnt(1, 1, 1));
  EXPECT_TRUE(BRep_ImplicitTool::IntersectRay(f, gp_Pnt(-5, 0, 0), gp_Vec(1, 0, 0), 0, 10).empty());
  EXPECT_THROW(BRep_ImplicitTool::Normal(f, gp_Pnt()), Standard_DomainError);
  EXPECT_THROW(BRep_ImplicitTool::IntersectRay(f, gp_Pnt(), gp_Vec(), 0, 10), Standard_DomainError);
  EXPECT_THROW(BRep_ImplicitTool::IntersectRay(f, gp_Pnt(), gp_Vec(1, 0, 0), 2, 1),
               Standard_DomainError);
  EXPECT_THROW(BRep_ImplicitTool::IntersectRay(f, gp_Pnt(), gp_Vec(1, 0, 0), 0, 10, 0),
               Standard_DomainError);
}

TEST(BRep_ImplicitFace, GyroidRayRootsWithoutUV)
{
  TopoDS_Face f;
  BRep_Builder().MakeImplicitFace(f,
                                  new Geom_ImplicitGyroid(10),
                                  gp_Pnt(-12, -1, -1),
                                  gp_Pnt(12, 1, 1),
                                  TopLoc_Location(),
                                  0.1);
  auto hits =
    BRep_ImplicitTool::IntersectRay(f, gp_Pnt(-12, 0, 0), gp_Vec(1, 0, 0), 0, 24, 127, 0.01);
  ASSERT_EQ(hits.size(), 5u);
  for (std::size_t i = 0; i < hits.size(); ++i)
    EXPECT_NEAR(hits[i].Point.X(), -10 + 5 * int(i), 0.01);
}

TEST(BRep_ImplicitFace, CopyPreservesFieldButDropsTransientMesh)
{
  auto                            f    = sphereFace();
  occ::handle<Poly_Triangulation> mesh = new Poly_Triangulation(3, 1, false);
  mesh->SetNode(1, gp_Pnt(2, 0, 0));
  mesh->SetNode(2, gp_Pnt(0, 2, 0));
  mesh->SetNode(3, gp_Pnt(0, 0, 2));
  mesh->SetTriangle(1, Poly_Triangle(1, 2, 3));
  BRep_Builder().UpdateFace(f, mesh);
  auto copy = TopoDS::Face(f.EmptyCopied());
  ASSERT_TRUE(BRep_Tool::IsImplicit(copy));
  TopLoc_Location a, b;
  EXPECT_EQ(BRep_Tool::ImplicitGeometry(f, a), BRep_Tool::ImplicitGeometry(copy, b));
  EXPECT_TRUE(BRep_Tool::Triangulation(copy, a).IsNull());
  EXPECT_FALSE(BRep_Tool::Triangulation(f, a).IsNull());
  EXPECT_DOUBLE_EQ(BRep_ImplicitTool::Value(copy, gp_Pnt(3, 0, 0)), 1);
  std::ostringstream json;
  EXPECT_NO_THROW(copy.TShape()->DumpJson(json));
  EXPECT_NE(json.str().find("BRep_TImplicitFace"), std::string::npos);
}

TEST(BRep_ImplicitFace, ParametricMutationAndLossyPersistenceRejected)
{
  auto         f = sphereFace();
  BRep_Builder builder;
  EXPECT_THROW(
    builder.UpdateFace(f, new Geom_Plane(gp_Pnt(), gp_Dir(0, 0, 1)), TopLoc_Location(), 0.1),
    Standard_NotImplemented);
  EXPECT_THROW(builder.NaturalRestriction(f, true), Standard_NotImplemented);
  EXPECT_THROW(builder.UpdateFace(f, -1.0), Standard_DomainError);
  EXPECT_DOUBLE_EQ(BRep_Tool::Tolerance(f), 0.1);
  EXPECT_TRUE(BRep_Tool::IsImplicit(f));
  std::ostringstream ascii, binary;
  EXPECT_THROW(BRepTools::Write(f, ascii), Standard_NotImplemented);
  EXPECT_THROW(BinTools::Write(f, binary), Standard_NotImplemented);
}

TEST(BRep_ImplicitFace, ConstructorValidationPreservesExistingFace)
{
  auto         f   = sphereFace();
  auto         old = f.TShape();
  BRep_Builder builder;
  EXPECT_THROW(builder.MakeImplicitFace(f,
                                        occ::handle<Geom_ImplicitField>(),
                                        gp_Pnt(-1, -1, -1),
                                        gp_Pnt(1, 1, 1),
                                        TopLoc_Location(),
                                        0.1),
               Standard_DomainError);
  EXPECT_EQ(old, f.TShape());
  EXPECT_THROW(builder.MakeImplicitFace(f,
                                        new Geom_ImplicitSphere(gp_Pnt(), 2),
                                        gp_Pnt(),
                                        gp_Pnt(),
                                        TopLoc_Location(),
                                        0.1),
               Standard_DomainError);
  EXPECT_EQ(old, f.TShape());
}

TEST(BRep_ImplicitFace, SamplingLimitsAreExplicit)
{
  auto sphere = sphereFace();
  EXPECT_TRUE(
    BRep_ImplicitTool::IntersectRay(sphere, gp_Pnt(-5, 2, 0), gp_Vec(1, 0, 0), 0, 10, 31).empty());
  EXPECT_EQ(
    BRep_ImplicitTool::IntersectRay(sphere, gp_Pnt(-5, 2, 0), gp_Vec(1, 0, 0), 0, 10, 32).size(),
    1u);
  TopoDS_Face plane;
  BRep_Builder().MakeImplicitFace(plane,
                                  new PlaneField(),
                                  gp_Pnt(-1, -1, -1),
                                  gp_Pnt(1, 1, 1),
                                  TopLoc_Location(),
                                  0.1);
  EXPECT_THROW(BRep_ImplicitTool::IntersectRay(plane, gp_Pnt(), gp_Vec(1, 0, 0), 0, 1),
               Standard_DomainError);
}
