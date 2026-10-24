// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRep_ImplicitCylinderSection.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_Builder.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Standard_AbortiveTransaction.hxx>
#include <Standard_DomainError.hxx>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
using Section = BRep_ImplicitCylinderSection;
TopoDS_Face field(double iso = .21)
{
  TopoDS_Face face;
  BRep_Builder().MakeImplicitFace(face, new Geom_ImplicitGyroid(10, iso),
    gp_Pnt(-10, -10, -10), gp_Pnt(10, 10, 10), TopLoc_Location(), .1);
  return face;
}
Section::Settings settings()
{
  Section::Settings s;
  s.Radius = 2;
  s.AngleMinimum = .27; s.AngleMaximum = 1.1;
  s.HeightMinimum = -4.2; s.HeightMaximum = 3.7;
  return s;
}
double distance(const gp_Pnt& point, const Section::Result& result)
{
  double closest = std::numeric_limits<double>::infinity();
  for (const auto& curve : result.Curves)
    for (size_t i = 1; i < curve.Points.size(); ++i)
    {
      const auto& a = curve.Points[i - 1];
      const gp_Vec ab(a, curve.Points[i]), ap(a, point);
      const double t = std::clamp(ap.Dot(ab) / ab.SquareMagnitude(), 0., 1.);
      closest = std::min(closest, point.Distance(a.Translated(ab * t)));
    }
  return closest;
}
}

TEST(BRepImplicitCylinderSection, AnalyticHeightRootsAndPhysicalCoverage)
{
  const auto face = field();
  const auto s = settings();
  const auto result = Section::Trace(face, s);
  ASSERT_TRUE(result.Resolved) << result.UnresolvedCells << " / " << result.GraphDefects;
  ASSERT_FALSE(result.Curves.empty());
  EXPECT_LE(result.GeometricBound, .1);
  EXPECT_LE(result.MaximumRootBracket, .1 / 32);
  for (const auto& curve : result.Curves)
    for (const auto& p : curve.Points)
    {
      EXPECT_NEAR(std::hypot(p.X(), p.Y()), s.Radius, 1.e-12);
      EXPECT_LE(std::abs(BRep_ImplicitTool::Value(face, p)), .01);
    }
  // Independent analytic height solution of A*cos(k*z)+B*sin(k*z)+C=0.
  const double k = 2 * std::acos(-1.) / 10;
  int witnesses = 0;
  for (int i = 0; i <= 60; ++i)
  {
    const double angle = s.AngleMinimum + (s.AngleMaximum - s.AngleMinimum) * i / 60;
    const double x = s.Radius * std::cos(angle), y = s.Radius * std::sin(angle);
    const double a = std::sin(k * y), b = std::cos(k * x);
    const double c = std::sin(k * x) * std::cos(k * y) - .21;
    const double radius = std::hypot(a, b);
    if (std::abs(c) > radius) continue;
    const double phase = std::atan2(b, a), offset = std::acos(-c / radius);
    for (double sign : {-1., 1.}) for (int turn = -2; turn <= 2; ++turn)
    {
      const double z = (phase + sign * offset + turn * 2 * std::acos(-1.)) / k;
      if (z < s.HeightMinimum || z > s.HeightMaximum) continue;
      ++witnesses;
      EXPECT_LE(distance(gp_Pnt(x, y, z), result), result.GeometricBound + 1.e-10);
    }
  }
  EXPECT_GT(witnesses, 10);
}

TEST(BRepImplicitCylinderSection, EmptyBudgetPlacementAndCancellation)
{
  auto s = settings();
  auto empty = Section::Trace(field(10), s);
  EXPECT_TRUE(empty.Resolved); EXPECT_TRUE(empty.Curves.empty());
  s.MaximumCells = 4;
  EXPECT_FALSE(Section::Trace(field(), s).Resolved);
  s = settings();
  auto face = field();
  gp_Trsf pose; pose.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), .4);
  pose.SetTranslationPart(gp_Vec(3, -4, 2));
  face.Location(TopLoc_Location(pose));
  s.Frame.Transform(pose);
  const auto placed = Section::Trace(face, s);
  ASSERT_TRUE(placed.Resolved);
  for (const auto& curve : placed.Curves) for (const auto& p : curve.Points)
    EXPECT_NEAR(std::hypot(p.Transformed(pose.Inverted()).X(),
                          p.Transformed(pose.Inverted()).Y()), s.Radius, 1.e-11);
  s.WorkCheckpoint = [] { throw Standard_AbortiveTransaction("Cylinder query cancelled"); };
  EXPECT_THROW(Section::Trace(face, s), Standard_AbortiveTransaction);
}

TEST(BRepImplicitCylinderSection, InvalidAndPeriodicChartsReject)
{
  const auto face = field();
  auto s = settings(); s.Radius = -1;
  EXPECT_THROW(Section::Trace(face, s), Standard_DomainError);
  s = settings(); s.AngleMinimum = 0; s.AngleMaximum = 2 * std::acos(-1.);
  EXPECT_THROW(Section::Trace(face, s), Standard_DomainError);
  s = settings(); s.Frame.SetLocation(gp_Pnt(20, 0, 0));
  EXPECT_THROW(Section::Trace(face, s), Standard_DomainError);
  s = settings(); s.DistanceTolerance = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(Section::Trace(face, s), Standard_DomainError);
}

TEST(BRepImplicitCylinderSection, PeriodicSharedSeamClosesOneAnalyticLoop)
{
  const auto face = field();
  auto s = settings();
  s.Radius = .5;
  s.AngleMinimum = .27;
  s.AngleMaximum = s.AngleMinimum + 2 * std::acos(-1.);
  s.HeightMinimum = -1.7; s.HeightMaximum = 1.3;
  const auto result = Section::TracePeriodic(face, s);
  ASSERT_TRUE(result.Resolved) << result.UnresolvedCells << " / " << result.GraphDefects;
  ASSERT_EQ(result.Curves.size(), 1);
  const auto& curve = result.Curves.front();
  EXPECT_TRUE(curve.Closed);
  ASSERT_GT(curve.Points.size(), 20);
  EXPECT_EQ(curve.Points.front().Distance(curve.Points.back()), 0);
  EXPECT_EQ(result.GraphDefects, 0);
  EXPECT_LE(result.GeometricBound, .1);
  const double k = 2 * std::acos(-1.) / 10;
  for (int i = 0; i <= 240; ++i)
  {
    const double angle = s.AngleMinimum + 2 * std::acos(-1.) * i / 240;
    const double x = s.Radius * std::cos(angle), y = s.Radius * std::sin(angle);
    const double a = std::sin(k * y), b = std::cos(k * x);
    const double c = std::sin(k * x) * std::cos(k * y) - .21;
    const double phase = std::atan2(b, a), offset = std::acos(-c / std::hypot(a, b));
    const double z = (phase - offset) / k;
    ASSERT_GE(z, s.HeightMinimum); ASSERT_LE(z, s.HeightMaximum);
    EXPECT_LE(distance(gp_Pnt(x, y, z), result), result.GeometricBound + 1.e-10);
  }
  // Moving the seam changes the chart only, not the physical closed section.
  s.AngleMinimum = .63; s.AngleMaximum = s.AngleMinimum + 2 * std::acos(-1.);
  const auto moved = Section::TracePeriodic(face, s);
  ASSERT_TRUE(moved.Resolved);
  ASSERT_EQ(moved.Curves.size(), 1);
  EXPECT_TRUE(moved.Curves.front().Closed);
  for (const auto& p : curve.Points) EXPECT_LE(distance(p, moved), .1);
  s.MaximumCells = 4;
  EXPECT_FALSE(Section::TracePeriodic(face, s).Resolved);
  s = settings();
  EXPECT_THROW(Section::TracePeriodic(face, s), Standard_DomainError);
}

TEST(BRepImplicitCylinderSection, PeriodicOpenEndsBelongOnlyToHeightBoundary)
{
  auto s = settings();
  s.Radius = .5;
  s.AngleMinimum = 3;
  s.AngleMaximum = s.AngleMinimum + 2 * std::acos(-1.);
  s.HeightMinimum = 0; s.HeightMaximum = 1.3;
  const auto result = Section::TracePeriodic(field(), s);
  ASSERT_TRUE(result.Resolved) << result.UnresolvedCells << " / " << result.GraphDefects;
  ASSERT_FALSE(result.Curves.empty());
  int open = 0;
  for (const auto& curve : result.Curves)
    if (!curve.Closed)
    {
      ++open;
      EXPECT_EQ(curve.Points.front().Z(), 0);
      EXPECT_EQ(curve.Points.back().Z(), 0);
    }
  EXPECT_GT(open, 0);
  // A near/exact zero on a grid/seam endpoint must remain unresolved.
  s.AngleMinimum = 0; s.AngleMaximum = 2 * std::acos(-1.);
  EXPECT_FALSE(Section::TracePeriodic(field(std::sin(std::acos(-1.) / 10)), s).Resolved);
}
