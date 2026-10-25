// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitCylinder.hxx>
#include <BRep_ImplicitIO.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <Standard_AbortiveTransaction.hxx>
#include <gtest/gtest.h>
#include <cmath>
#include <sstream>

namespace
{
using Ops = BRepClass3d_ImplicitCylinder;
using Settings = BRep_ImplicitCylinderSection::Settings;
TopoDS_Face field(double iso = .21)
{
  TopoDS_Face face;
  BRep_Builder().MakeImplicitFace(face, new Geom_ImplicitGyroid(10, iso),
    gp_Pnt(-10, -10, -10), gp_Pnt(10, 10, 10), TopLoc_Location(), .1);
  return face;
}
Settings settings()
{
  Settings s; s.Radius = .5; s.AngleMinimum = 3;
  s.AngleMaximum = s.AngleMinimum + 2 * std::acos(-1.);
  s.HeightMinimum = 0; s.HeightMaximum = .5;
  return s;
}
double area(const std::vector<TopoDS_Face>& faces)
{
  double value = 0;
  for (const auto& face : faces)
  {
    GProp_GProps props;
    const double error = BRepGProp::SurfaceProperties(face, props, 1.e-10);
    EXPECT_LE(error, 1.e-8);
    value += props.Mass();
  }
  return value;
}
int wires(const TopoDS_Face& face)
{
  int count = 0;
  for (TopExp_Explorer it(face, TopAbs_WIRE); it.More(); it.Next()) ++count;
  return count;
}
}

TEST(BRepClass3dImplicitCylinder, NativeCapsPartitionDiskAndShareSideRootVertices)
{
  const auto source = field(); const auto s = settings();
  std::stringstream before; BRep_ImplicitIO::Write(source, before);
  const auto negative = Ops::BuildCaps(source, s, true);
  const auto positive = Ops::BuildCaps(source, s, false);
  ASSERT_TRUE(negative.TopologyValidated) << negative.Reason;
  ASSERT_TRUE(positive.TopologyValidated) << positive.Reason;
  EXPECT_LE(negative.ContourBound, .1);
  EXPECT_TRUE(negative.Field.IsSame(source));
  EXPECT_EQ(negative.Parameters.WorkCheckpoint, nullptr);
  for (int cap = 0; cap < 2; ++cap)
  {
    ASSERT_EQ(negative.Faces[cap].size(), 1); ASSERT_EQ(positive.Faces[cap].size(), 1);
    EXPECT_NEAR(area(negative.Faces[cap]) + area(positive.Faces[cap]), std::acos(-1.) * .25, 1.e-6);
    for (const auto* result : {&negative, &positive})
    {
      const auto& face = result->Faces[cap].front();
      EXPECT_TRUE(BRepCheck_Analyzer(face).IsValid());
      EXPECT_EQ(wires(face), 1);
      const auto normal = BRepAdaptor_Surface(face).Plane().Axis().Direction();
      EXPECT_NEAR(normal.Dot(s.Frame.Direction()), cap == 0 ? -1 : 1, 1.e-12);
      for (const auto& edge : result->FieldEdges[cap])
      {
        TopoDS_Vertex a, b; TopExp::Vertices(edge, a, b);
        for (const auto& v : {a, b})
        {
          bool shared = false;
          for (const auto& root : result->Source.Disks.Boundary.Roots) shared |= v.IsSame(root.Vertex);
          EXPECT_TRUE(shared);
        }
      }
      int witnesses = 0;
      for (int i = -4; i <= 4; ++i) for (int j = -4; j <= 4; ++j)
      {
        const gp_Pnt p(i * .1, j * .1, cap == 0 ? 0 : .5);
        if (std::hypot(p.X(), p.Y()) >= .48) continue;
        const double value = BRep_ImplicitTool::Value(source, p);
        if (std::abs(value) <= 2.2 * result->ContourBound) continue;
        ++witnesses;
        const auto state = BRepClass_FaceClassifier(face, p, 1.e-7).State();
        EXPECT_EQ(state == TopAbs_IN, result->Negative ? value < 0 : value > 0);
      }
      EXPECT_GT(witnesses, 10);
    }
  }
  std::stringstream after; BRep_ImplicitIO::Write(source, after);
  EXPECT_EQ(before.str(), after.str());
}

TEST(BRepClass3dImplicitCylinder, NativeCapsNestInteriorHolesAndPositiveIslands)
{
  auto s = settings(); s.Frame.SetLocation(gp_Pnt(2.5, 1.25, 0));
  s.Radius = 1.5; s.HeightMaximum = .02;
  const auto source = field(1.3);
  const auto negative = Ops::BuildCaps(source, s, true);
  const auto positive = Ops::BuildCaps(source, s, false);
  ASSERT_TRUE(negative.TopologyValidated) << negative.Reason;
  ASSERT_TRUE(positive.TopologyValidated) << positive.Reason;
  for (int cap = 0; cap < 2; ++cap)
  {
    ASSERT_EQ(negative.Faces[cap].size(), 1); ASSERT_EQ(positive.Faces[cap].size(), 1);
    EXPECT_EQ(wires(negative.Faces[cap][0]), 2); EXPECT_EQ(wires(positive.Faces[cap][0]), 1);
    EXPECT_NEAR(area(negative.Faces[cap]) + area(positive.Faces[cap]), std::acos(-1.) * 2.25, 1.e-5);
    const gp_Pnt center(2.5, 1.25, cap == 0 ? 0 : .02);
    EXPECT_EQ(BRepClass_FaceClassifier(negative.Faces[cap][0], center, 1.e-7).State(), TopAbs_OUT);
    EXPECT_EQ(BRepClass_FaceClassifier(positive.Faces[cap][0], center, 1.e-7).State(), TopAbs_IN);
  }
}

TEST(BRepClass3dImplicitCylinder, NativeCapsFullEmptyFailureAndCancellation)
{
  auto s = settings();
  const auto full = Ops::BuildCaps(field(10), s, true);
  const auto empty = Ops::BuildCaps(field(10), s, false);
  ASSERT_TRUE(full.TopologyValidated) << full.Reason;
  ASSERT_TRUE(empty.TopologyValidated) << empty.Reason;
  for (int cap = 0; cap < 2; ++cap)
  { EXPECT_NEAR(area(full.Faces[cap]), std::acos(-1.) * .25, 1.e-7); EXPECT_TRUE(empty.Faces[cap].empty()); }
  const auto lateFailure = Ops::BuildCaps(field(10), s, true, 1);
  EXPECT_FALSE(lateFailure.TopologyValidated); EXPECT_EQ(lateFailure.Reason, "cap_nesting_budget");
  for (int cap = 0; cap < 2; ++cap)
  { EXPECT_TRUE(lateFailure.Faces[cap].empty()); EXPECT_TRUE(lateFailure.FieldEdges[cap].empty()); }
  s.MaximumCells = 4;
  const auto failed = Ops::BuildCaps(field(), s);
  EXPECT_FALSE(failed.TopologyValidated);
  for (int cap = 0; cap < 2; ++cap)
  { EXPECT_TRUE(failed.Faces[cap].empty()); EXPECT_TRUE(failed.FieldEdges[cap].empty()); }
  s = settings();
  s.WorkCheckpoint = [] { throw Standard_AbortiveTransaction("Cancel native caps"); };
  EXPECT_THROW(Ops::BuildCaps(field(), s), Standard_AbortiveTransaction);
}

TEST(BRepClass3dImplicitCylinder, NativeCapsRigidPlacementAndCriticalContactRefusal)
{
  auto s = settings(); auto source = field();
  gp_Trsf pose; pose.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), .4);
  pose.SetTranslationPart(gp_Vec(3, -4, 2));
  source.Location(TopLoc_Location(pose)); s.Frame.Transform(pose);
  const auto placed = Ops::BuildCaps(source, s);
  ASSERT_TRUE(placed.TopologyValidated) << placed.Reason;
  for (const auto& faces : placed.Faces) for (const auto& face : faces)
    EXPECT_TRUE(BRepCheck_Analyzer(face).IsValid());
  s = settings(); s.Frame.SetLocation(gp_Pnt(2.5, 1.25, 0)); s.Radius = 1.5;
  s.HeightMinimum = -.5; s.HeightMaximum = 0;
  const auto critical = Ops::BuildCaps(field(std::sqrt(2.)), s);
  EXPECT_FALSE(critical.TopologyValidated);
  for (const auto& faces : critical.Faces) EXPECT_TRUE(faces.empty());
}

TEST(BRepClass3dImplicitCylinder, NativeCapsRejectUnrepresentableArea)
{
  TopoDS_Face source;
  BRep_Builder().MakeImplicitFace(source, new Geom_ImplicitGyroid(1.e155, 10),
    gp_Pnt(-2.e155, -2.e155, -2.e155), gp_Pnt(2.e155, 2.e155, 2.e155), TopLoc_Location(), .1);
  auto s = settings(); s.Radius = 1.e155; s.DistanceTolerance = 1.e145;
  const auto result = Ops::BuildCaps(source, s);
  ASSERT_TRUE(result.Source.Verified);
  EXPECT_FALSE(result.TopologyValidated);
  EXPECT_EQ(result.Reason, "cap_cycle_area_unresolved");
  for (const auto& faces : result.Faces) EXPECT_TRUE(faces.empty());
}
