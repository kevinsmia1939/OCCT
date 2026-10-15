// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRep_TImplicitFace.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_ImplicitIO.hxx>
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
#include <TopoDS_Wire.hxx>
#include <TopoDS_LockedShape.hxx>
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

TEST(BRep_ImplicitFace, PersistenceRoundTrip)
{
  for (int primitive = 0; primitive < 2; ++primitive)
    for (int orientation = 0; orientation < 4; ++orientation)
    {
      occ::handle<Geom_ImplicitField> field;
      if (primitive == 0)
        field = new Geom_ImplicitSphere(gp_Pnt(0.31, -0.23, 0.57), 2.123456789012345);
      else
        field = new Geom_ImplicitGyroid(3.123456789012345, -0.2718281828459);
      gp_Trsf geometry;
      geometry.SetTranslation(gp_Vec(1.125, 2.25, -3.5));
      gp_Trsf instance;
      instance.SetRotation(gp_Ax1(gp_Pnt(1, 2, 3), gp_Dir(1, 2, 3)), 0.731);
      gp_Trsf scale;
      scale.SetScale(gp_Pnt(0, 0, 0), -1.25);
      instance.Multiply(scale);
      TopoDS_Face source;
      BRep_Builder().MakeImplicitFace(source,
                                      field,
                                      gp_Pnt(-4, -5, -6),
                                      gp_Pnt(7, 8, 9),
                                      TopLoc_Location(geometry),
                                      0.1);
      source.Location(TopLoc_Location(instance), false);
      source.Orientation(TopAbs_Orientation(orientation));
      BRep_Builder().UpdateFace(source, new Poly_Triangulation(3, 1, false));
      std::ostringstream out;
      BRep_ImplicitIO::Write(source, out);
      std::istringstream in(out.str());
      TopoDS_Face        restored;
      BRep_ImplicitIO::Read(in, restored);
      EXPECT_EQ(restored.Orientation(), source.Orientation());
      const auto* data = static_cast<const BRep_TImplicitFace*>(restored.TShape().get());
      EXPECT_DOUBLE_EQ(data->Tolerance(), 0.1);
      EXPECT_DOUBLE_EQ(data->Minimum().Y(), -5);
      EXPECT_DOUBLE_EQ(data->Maximum().Z(), 9);
      TopLoc_Location location;
      EXPECT_TRUE(BRep_Tool::Triangulation(restored, location).IsNull());
      EXPECT_EQ(typeid(*data->Field()), typeid(*field));
      for (int r = 1; r <= 3; ++r)
        for (int c = 1; c <= 4; ++c)
        {
          EXPECT_NEAR(data->Location().Transformation().Value(r, c), geometry.Value(r, c), 1.e-14);
          EXPECT_NEAR(restored.Location().Transformation().Value(r, c),
                      instance.Value(r, c),
                      1.e-14);
        }
      for (const gp_Pnt p : {gp_Pnt(0.123, 0.456, 0.789), gp_Pnt(-1, 2, -3)})
      {
        EXPECT_NEAR(BRep_ImplicitTool::Value(source, p),
                    BRep_ImplicitTool::Value(restored, p),
                    1.e-13);
        EXPECT_NEAR(
          (BRep_ImplicitTool::Gradient(source, p) - BRep_ImplicitTool::Gradient(restored, p))
            .Magnitude(),
          0,
          1.e-12);
      }
      std::ostringstream again;
      BRep_ImplicitIO::Write(restored, again);
      TopoDS_Face        twice;
      std::istringstream second(again.str());
      EXPECT_NO_THROW(BRep_ImplicitIO::Read(second, twice));
    }
}

TEST(BRep_ImplicitFace, PersistenceRejectsWithoutDataLoss)
{
  TopoDS_Face        face = sphereFace();
  std::ostringstream out;
  BRep_ImplicitIO::Write(face, out);
  const std::string valid     = out.str();
  auto              badRecord = [&](const std::string& record) {
    TopoDS_Face        target = face;
    std::istringstream in(record);
    EXPECT_THROW(BRep_ImplicitIO::Read(in, target), Standard_Failure);
    EXPECT_TRUE(target.IsEqual(face));
  };
  for (size_t end = 0; end < valid.find("END") + 3; ++end)
    badRecord(valid.substr(0, end));
  badRecord(valid + "unexpected");
  badRecord(std::string(65537, ' '));
  auto replace = [&](const std::string& from, const std::string& to) {
    std::string changed = valid;
    changed.replace(changed.find(from), from.size(), to);
    badRecord(changed);
  };
  replace("FACE 1", "FACE 2");
  replace("SPHERE", "CUSTOM");
  replace("TOLERANCE 0.1", "TOLERANCE -0.1");
  replace("TOLERANCE 0.1", "TOLERANCE nan");
  replace("ORIENTATION 0", "ORIENTATION 4");
  replace("GEOMETRY_LOCATION\n1 0 0 0", "GEOMETRY_LOCATION\n1 0.2 0 0");
  replace("GEOMETRY_LOCATION\n1 0 0 0", "GEOMETRY_LOCATION\n2 0 0 0");
  replace("GEOMETRY_LOCATION\n1 0 0 0", "GEOMETRY_LOCATION\n0 0 0 0");

  TopoDS_Face locked = sphereFace();
  locked.Locked(true);
  std::istringstream in(valid);
  EXPECT_THROW(BRep_ImplicitIO::Read(in, locked), TopoDS_LockedShape);
  EXPECT_EQ(in.tellg(), std::streampos(0));

  auto unsupported = [&](const TopoDS_Face& source) {
    std::ostringstream stream;
    stream << "unchanged";
    EXPECT_THROW(BRep_ImplicitIO::Write(source, stream), Standard_NotImplemented);
    EXPECT_EQ(stream.str(), "unchanged");
  };
  unsupported(TopoDS_Face());
  TopoDS_Face  custom;
  BRep_Builder builder;
  builder.MakeImplicitFace(custom,
                           new PlaneField,
                           gp_Pnt(-1, -1, -1),
                           gp_Pnt(1, 1, 1),
                           TopLoc_Location(),
                           0.1);
  unsupported(custom);

  // A derived primitive may replace its equation: do not serialize its base.
  class DifferentSphere : public Geom_ImplicitSphere
  {
  public:
    DifferentSphere()
        : Geom_ImplicitSphere(gp_Pnt(), 2)
    {
    }

    double Value(const gp_Pnt&) const override { return 1; }
  };

  builder.MakeImplicitFace(custom,
                           new DifferentSphere,
                           gp_Pnt(-1, -1, -1),
                           gp_Pnt(1, 1, 1),
                           TopLoc_Location(),
                           0.1);
  unsupported(custom);
  TopoDS_Wire wire;
  builder.MakeWire(wire);
  builder.Add(face, wire);
  unsupported(face);
  TopoDS_Face conventional;
  builder.MakeFace(conventional, new Geom_Plane(gp_Pnt(), gp_Dir(0, 0, 1)), 0.1);
  unsupported(conventional);
}

TEST(BRep_ImplicitFace, PersistenceStreamAndTolerancePolicy)
{
  // The wire format must not inherit the caller's decimal locale.
  class CommaDecimal : public std::numpunct<char>
  {
    char do_decimal_point() const override { return ','; }
  };

  TopoDS_Face        source = sphereFace();
  std::ostringstream out;
  out.imbue(std::locale(std::locale::classic(), new CommaDecimal));
  BRep_ImplicitIO::Write(source, out);
  EXPECT_NE(out.str().find("TOLERANCE 0.10000000000000001"), std::string::npos);
  std::istringstream in(out.str());
  in.imbue(out.getloc());
  TopoDS_Face restored;
  BRep_ImplicitIO::Read(in, restored);
  const auto hits =
    BRep_ImplicitTool::IntersectRay(restored, gp_Pnt(-5, .31, .27), gp_Vec(1, 0, 0), 0, 10, 11);
  ASSERT_EQ(hits.size(), 2u);
  const double radiusX = std::sqrt(4 - .31 * .31 - .27 * .27);
  EXPECT_LE(std::abs(hits[0].Distance - (5 - radiusX)), .05);
  EXPECT_LE(std::abs(hits[1].Distance - (5 + radiusX)), .05);
  for (const auto& hit : hits)
    EXPECT_LE(hit.BracketWidth, .1);

  std::ostringstream failedOut;
  failedOut.setstate(std::ios::badbit);
  EXPECT_THROW(BRep_ImplicitIO::Write(source, failedOut), Standard_DomainError);
  std::istringstream failedIn(out.str());
  failedIn.setstate(std::ios::badbit);
  EXPECT_THROW(BRep_ImplicitIO::Read(failedIn, restored), Standard_DomainError);
  EXPECT_TRUE(BRep_Tool::IsImplicit(restored));

  // Persist the stored tolerance, not BRep_Tool's conventional precision floor.
  source     = sphereFace();
  auto* data = static_cast<BRep_TImplicitFace*>(source.TShape().get());
  data->Tolerance(1.e-10);
  std::ostringstream tiny;
  BRep_ImplicitIO::Write(source, tiny);
  std::istringstream tinyIn(tiny.str());
  BRep_ImplicitIO::Read(tinyIn, restored);
  EXPECT_DOUBLE_EQ(static_cast<const BRep_TImplicitFace*>(restored.TShape().get())->Tolerance(),
                   1.e-10);
  data->Tolerance(-1);
  std::ostringstream unchanged;
  unchanged << "sentinel";
  EXPECT_THROW(BRep_ImplicitIO::Write(source, unchanged), Standard_DomainError);
  EXPECT_EQ(unchanged.str(), "sentinel");
}

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
