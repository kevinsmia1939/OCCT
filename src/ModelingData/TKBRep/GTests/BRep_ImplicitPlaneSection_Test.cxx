// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRep_ImplicitPlaneSection.hxx>
#include <BRep_ImplicitGyroidOps.hxx>
#include <BRep_TImplicitPlanePatch.hxx>
#include <BRep_ImplicitIO.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Edge.hxx>
#include <BRep_Tool.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_Builder.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_ImplicitSphere.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <gtest/gtest.h>
#include <cmath>
#include <sstream>

namespace
{
TopoDS_Face gyroid(double iso = .21)
{
  TopoDS_Face face;
  BRep_Builder().MakeImplicitFace(face,
                                  new Geom_ImplicitGyroid(10, iso),
                                  gp_Pnt(-10, -10, -10),
                                  gp_Pnt(10, 10, 10),
                                  TopLoc_Location(),
                                  .1);
  return face;
}

BRep_ImplicitPlaneSection::Settings settings()
{
  BRep_ImplicitPlaneSection::Settings s;
  s.Frame    = gp_Ax2(gp_Pnt(0, 0, .37), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
  s.UMinimum = -3.83;
  s.UMaximum = 4.17;
  s.VMinimum = -4.07;
  s.VMaximum = 3.93;
  return s;
}
} // namespace

TEST(BRep_ImplicitPlaneSection, RegularGyroidGeometricBoundsAndConnectivity)
{
  const auto face   = gyroid();
  const auto s      = settings();
  const auto result = BRep_ImplicitPlaneSection::Trace(face, s);
  EXPECT_TRUE(result.Resolved) << result.UnresolvedCells << " unresolved, " << result.AmbiguousCells
                               << " ambiguous";
  ASSERT_FALSE(result.Curves.empty());
  EXPECT_LE(result.GeometricBound, .1);
  EXPECT_LE(result.MaximumRootBracket, .1 / 16);
  EXPECT_EQ(result.UnresolvedCells, 0);
  for (const auto& curve : result.Curves)
  {
    ASSERT_GE(curve.Points.size(), 2u);
    for (const auto& p : curve.Points)
    {
      EXPECT_NEAR(p.Z(), .37, 1.e-14);
      EXPECT_LE(std::abs(BRep_ImplicitTool::Value(face, p)), result.MaximumResidual + 1.e-14);
      EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(face, p));
    }
    if (!curve.Closed)
      for (const auto& p : {curve.Points.front(), curve.Points.back()})
        EXPECT_LT(std::min({std::abs(p.X() - s.UMinimum),
                            std::abs(p.X() - s.UMaximum),
                            std::abs(p.Y() - s.VMinimum),
                            std::abs(p.Y() - s.VMaximum)}),
                  1.e-12);
  }
  const auto again = BRep_ImplicitPlaneSection::Trace(face, s);
  ASSERT_EQ(again.Curves.size(), result.Curves.size());
  for (size_t i = 0; i < result.Curves.size(); ++i)
  {
    ASSERT_EQ(again.Curves[i].Points.size(), result.Curves[i].Points.size());
    for (size_t j = 0; j < result.Curves[i].Points.size(); ++j)
      EXPECT_DOUBLE_EQ(again.Curves[i].Points[j].Distance(result.Curves[i].Points[j]), 0);
  }
}

TEST(BRep_ImplicitPlaneSection, ExactLatticeNodesTraceAllDefaultCubeFacesAndKeepPointContacts)
{
  const auto face = gyroid(0);
  std::ostringstream before, after;
  BRep_ImplicitIO::Write(face, before);
  const double k = 2 * std::acos(-1.) / 10;
  for (int axis = 1; axis <= 3; ++axis)
    for (int high = 0; high < 2; ++high)
    {
      gp_Pnt origin;
      origin.SetCoord(axis, high ? 10 : 0);
      gp_Vec normal, u;
      normal.SetCoord(axis, 1);
      u.SetCoord(axis % 3 + 1, 1);
      BRep_ImplicitPlaneSection::Settings s;
      s.Frame = gp_Ax2(origin, gp_Dir(normal), gp_Dir(u));
      s.UMinimum = s.VMinimum = 0;
      s.UMaximum = s.VMaximum = 10;
      s.AllowExactLatticeNodes = true;
      const auto result = BRep_ImplicitPlaneSection::Trace(face, s);
      ASSERT_TRUE(result.Resolved) << axis << " " << high << ": " << result.UnresolvedCells
                                  << " ambiguous=" << result.AmbiguousCells
                                  << " graph=" << result.GraphDefects;
      ASSERT_EQ(result.Curves.size(), 3u);
      ASSERT_EQ(result.PointContacts.size(), 2u);
      EXPECT_GT(result.ExactLatticeNodeCount, 0);
      EXPECT_LE(result.GeometricBound, .1);
      EXPECT_LE(result.MaximumRootBracket, .1 / 16);
      const auto v = gp_Vec(s.Frame.YDirection());
      EXPECT_EQ(result.PointContacts[0].Distance(origin), 0);
      EXPECT_EQ(result.PointContacts[1].Distance(origin.Translated(u * 10 + v * 10)), 0);
      int cornerEnds = 0;
      for (const auto& curve : result.Curves)
      {
        ASSERT_FALSE(curve.Closed);
        ASSERT_GE(curve.Points.size(), 2u);
        for (const auto& p : {curve.Points.front(), curve.Points.back()})
          cornerEnds += p.Distance(origin.Translated(u * 10)) == 0
                     || p.Distance(origin.Translated(v * 10)) == 0;
        for (const auto& p : curve.Points)
          EXPECT_LE(std::abs(BRep_ImplicitTool::Value(face, p)), result.MaximumResidual + 1.e-13);
      }
      EXPECT_EQ(cornerEnds, 2);
      // Independent complete slice families v=(-atan(sin(k*u))+m*pi)/k.
      // These check geometric coverage, not just residuals at traced samples.
      for (int sample = 0; sample <= 100; ++sample)
        for (int turn = 0; turn <= 2; ++turn)
        {
          const double x = sample * .1;
          const double y = (-std::atan(std::sin(k * x)) + turn * std::acos(-1.)) / k;
          if (y < 0 || y > 10)
            continue;
          const auto p = origin.Translated(u * x + v * y);
          double nearest = 1.e100;
          for (const auto& curve : result.Curves)
            for (size_t i = 1; i < curve.Points.size(); ++i)
            {
              const auto a = curve.Points[i - 1];
              const gp_Vec ab(a, curve.Points[i]);
              const double t = std::clamp(gp_Vec(a, p).Dot(ab) / ab.SquareMagnitude(), 0., 1.);
              nearest = std::min(nearest, p.Distance(a.Translated(ab * t)));
            }
          for (const auto& contact : result.PointContacts)
            nearest = std::min(nearest, p.Distance(contact));
          EXPECT_LE(nearest, result.GeometricBound + 1.e-12);
        }
    }
  BRep_ImplicitPlaneSection::Settings tiny;
  tiny.Frame = gp_Ax2(gp_Pnt(-.125, -.125, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
  tiny.UMinimum = tiny.VMinimum = 0;
  tiny.UMaximum = tiny.VMaximum = .25;
  tiny.AllowExactLatticeNodes = true;
  const auto crossing = BRep_ImplicitPlaneSection::Trace(face, tiny);
  ASSERT_TRUE(crossing.Resolved);
  ASSERT_EQ(crossing.Curves.size(), 1u);
  EXPECT_TRUE(crossing.PointContacts.empty());
  int originUses = 0;
  for (const auto& p : crossing.Curves.front().Points)
    originUses += p.Distance(gp_Pnt()) == 0;
  EXPECT_EQ(originUses, 1); // One shared interior vertex, not four coincident roots.
  tiny.Frame = gp_Ax2(gp_Pnt(), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
  const auto contact = BRep_ImplicitPlaneSection::Trace(face, tiny);
  ASSERT_TRUE(contact.Resolved);
  EXPECT_TRUE(contact.Curves.empty());
  ASSERT_EQ(contact.PointContacts.size(), 1u);
  EXPECT_EQ(contact.PointContacts.front().Distance(gp_Pnt()), 0);
  tiny.Frame = gp_Ax2(gp_Pnt(0, 10, 0), gp_Dir(0, 0, -1), gp_Dir(1, 0, 0));
  tiny.UMaximum = tiny.VMaximum = 10;
  const auto reversed = BRep_ImplicitPlaneSection::Trace(face, tiny);
  ASSERT_TRUE(reversed.Resolved);
  ASSERT_EQ(reversed.Curves.size(), 3u);
  ASSERT_EQ(reversed.PointContacts.size(), 2u);
  EXPECT_EQ(reversed.PointContacts[0].Distance(gp_Pnt(10, 10, 0)), 0);
  EXPECT_EQ(reversed.PointContacts[1].Distance(gp_Pnt()), 0);
  BRep_ImplicitIO::Write(face, after);
  EXPECT_EQ(before.str(), after.str());
}

TEST(BRep_ImplicitPlaneSection, LatticeModeRefusesUnknownZerosRoundedAffineMapsAndPlacements)
{
  auto s = settings();
  s.Frame = gp_Ax2(gp_Pnt(), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
  s.UMinimum = s.VMinimum = 0;
  s.UMaximum = s.VMaximum = 10;
  s.MaximumCells = 64;
  s.AllowExactLatticeNodes = true;
  const auto budget = BRep_ImplicitPlaneSection::Trace(gyroid(0), s);
  EXPECT_FALSE(budget.Resolved);
  EXPECT_GT(budget.UnresolvedCells, 0);
  s.Frame = gp_Ax2(gp_Pnt(.1, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
  s.UMinimum = -.1;
  s.UMaximum = 9.9; // Rounded sum is 10, but the exact requested endpoint is not.
  const auto rounded = BRep_ImplicitPlaneSection::Trace(gyroid(0), s);
  EXPECT_FALSE(rounded.Resolved);
  EXPECT_GT(rounded.AmbiguousCells, 0);
  for (const auto& p : rounded.PointContacts)
    EXPECT_NE(p.X(), 10);
  s.Frame = gp_Ax2(gp_Pnt(), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
  s.UMinimum = 0;
  s.UMaximum = std::nextafter(10., 0.);
  EXPECT_FALSE(BRep_ImplicitPlaneSection::Trace(gyroid(0), s).Resolved);
  s.UMaximum = 10;
  s.Frame = gp_Ax2(gp_Pnt(), gp_Dir(0, 0, 1), gp_Dir(1, 1, 0));
  EXPECT_THROW(BRep_ImplicitPlaneSection::Trace(gyroid(0), s), Standard_NotImplemented);
  s.Frame = gp_Ax2(gp_Pnt(), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
  gp_Trsf shift;
  shift.SetTranslation(gp_Vec(1.e-18, 0, 0));
  EXPECT_THROW(BRep_ImplicitPlaneSection::Trace(TopoDS::Face(gyroid(0).Moved(TopLoc_Location(shift))), s),
               Standard_NotImplemented);
}

TEST(BRep_ImplicitPlaneSection, EmptyVersusUnresolvedAndInvalid)
{
  auto       s     = settings();
  const auto empty = BRep_ImplicitPlaneSection::Trace(gyroid(3), s);
  EXPECT_TRUE(empty.Resolved);
  EXPECT_TRUE(empty.Curves.empty());
  s.MaximumCells    = 4;
  const auto budget = BRep_ImplicitPlaneSection::Trace(gyroid(), s);
  EXPECT_FALSE(budget.Resolved);
  EXPECT_GT(budget.UnresolvedCells, 0);
  s              = settings();
  s.MaximumCells = 32;
  s.Frame        = gp_Ax2(gp_Pnt(), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
  s.UMinimum = s.VMinimum = -5;
  s.UMaximum = s.VMaximum = 5;
  const auto aligned      = BRep_ImplicitPlaneSection::Trace(gyroid(0), s);
  EXPECT_FALSE(aligned.Resolved);
  EXPECT_GT(aligned.AmbiguousCells, 0);
  s                   = settings();
  s.DistanceTolerance = 0;
  EXPECT_THROW(BRep_ImplicitPlaneSection::Trace(gyroid(), s), Standard_DomainError);
  s              = settings();
  s.MaximumCells = 1024;
  EXPECT_THROW(BRep_ImplicitPlaneSection::Trace(gyroid(), s), Standard_DomainError);
  s          = settings();
  s.UMinimum = -11;
  EXPECT_THROW(BRep_ImplicitPlaneSection::Trace(gyroid(), s), Standard_DomainError);
  s = settings();
  EXPECT_THROW(BRep_ImplicitPlaneSection::Trace(TopoDS_Face(), s), Standard_NotImplemented);
}

TEST(BRep_ImplicitPlaneSection, PlacedObliqueRectangle)
{
  const auto source = gyroid();
  auto       placed = source;
  gp_Trsf    t;
  t.SetRotation(gp_Ax1(gp_Pnt(), gp_Dir(1, 2, 3)), .53);
  gp_Trsf scale;
  scale.SetScale(gp_Pnt(), 1.25);
  t.Multiply(scale);
  t.SetTranslationPart(gp_Vec(2, 3, 4));
  placed.Location(TopLoc_Location(t));
  auto s = settings();
  s.Frame.Transform(t);
  s.UMinimum *= 1.25;
  s.UMaximum *= 1.25;
  s.VMinimum *= 1.25;
  s.VMaximum *= 1.25;
  const auto result = BRep_ImplicitPlaneSection::Trace(placed, s);
  EXPECT_TRUE(result.Resolved) << result.UnresolvedCells;
  ASSERT_FALSE(result.Curves.empty());
  for (const auto& curve : result.Curves)
    for (const auto& p : curve.Points)
    {
      const auto local = p.Transformed(t.Inverted());
      EXPECT_NEAR(local.Z(), .37, 1.e-13);
      EXPECT_LE(std::abs(BRep_ImplicitTool::Value(placed, p)), result.MaximumResidual + 1.e-13);
    }
}

TEST(BRep_ImplicitPlaneSection, ClosedGyroidSplitSharesApproximateSeam)
{
  TopoDS_Face face;
  BRep_Builder().MakeImplicitFace(face,
                                  new Geom_ImplicitGyroid(10, 1.45),
                                  gp_Pnt(-.31, -.47, -4),
                                  gp_Pnt(2.93, 3.07, 4),
                                  TopLoc_Location(),
                                  .1);
  const auto split = BRep_ImplicitGyroidOps::SplitZ(face, 1.25);
  ASSERT_TRUE(split.Section.Resolved);
  ASSERT_FALSE(split.Section.Curves.empty());
  for (const auto& curve : split.Section.Curves)
    EXPECT_TRUE(curve.Closed);
  auto wireA = TopoDS_Iterator(split.Negative), wireB = TopoDS_Iterator(split.Positive);
  for (; wireA.More() && wireB.More(); wireA.Next(), wireB.Next())
  {
    const auto a = TopoDS::Edge(TopoDS_Iterator(wireA.Value()).Value());
    const auto b = TopoDS::Edge(TopoDS_Iterator(wireB.Value()).Value());
    EXPECT_TRUE(a.IsSame(b));
    EXPECT_NE(a.Orientation(), b.Orientation());
    double     first, last;
    const auto curve = BRep_Tool::Curve(a, first, last);
    for (int i = 0; i < 100; ++i)
    {
      const auto p = curve->Value(first + (last - first) * i / 100);
      EXPECT_NEAR(p.Z(), 1.25, 1.e-14);
      // Distance control is the section's cell bound, NOT this residual.
      EXPECT_LE(std::abs(BRep_ImplicitTool::Value(face, p)),
                2 * std::sqrt(3.) * (2 * std::acos(-1.) / 10) * split.Section.GeometricBound
                  + 1.e-12);
    }
  }
  EXPECT_FALSE(wireA.More());
  EXPECT_FALSE(wireB.More());
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(split.Negative, gp_Pnt(1, 1, 0)));
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(split.Negative, gp_Pnt(1, 1, 2)));
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(split.Positive, gp_Pnt(1, 1, 2)));
  auto data = occ::down_cast<BRep_TImplicitPlanePatch>(split.Negative.TShape());
  ASSERT_FALSE(data.IsNull());
  EXPECT_FALSE(occ::down_cast<BRep_TImplicitPlanePatch>(data->EmptyCopy()).IsNull());
  std::ostringstream out;
  EXPECT_THROW(BRep_ImplicitIO::Write(split.Negative, out), Standard_NotImplemented);
  EXPECT_TRUE(out.str().empty());
  EXPECT_FALSE(TopoDS_Iterator(face).More());
  EXPECT_THROW(BRep_ImplicitGyroidOps::SplitZ(gyroid(), .37), Standard_NotImplemented);
  EXPECT_THROW(BRep_ImplicitGyroidOps::SplitZ(face, 1.25, .1, 4), Standard_DomainError);
}

TEST(BRep_ImplicitPlaneSection, IndependentAnalyticSliceCoverageAndSaddleRejection)
{
  const auto s      = settings();
  const auto result = BRep_ImplicitPlaneSection::Trace(gyroid(), s);
  ASSERT_TRUE(result.Resolved);
  const double k       = 2 * std::acos(-1.) / 10;
  int          checked = 0;
  // On z=constant, A*cos(k*y)+B*sin(k*y)=C gives independent exact
  // reference points; do not validate a geometric bound using field residuals.
  for (int i = 0; i <= 500; ++i)
  {
    const double x = s.UMinimum + (s.UMaximum - s.UMinimum) * i / 500;
    const double a = std::sin(k * x), b = std::cos(k * .37);
    const double c = .21 - std::sin(k * .37) * std::cos(k * x), radius = std::hypot(a, b);
    if (std::abs(c) > radius)
      continue;
    for (int turn = -2; turn <= 2; ++turn)
      for (double sign : {-1., 1.})
      {
        const double y =
          (std::atan2(b, a) + sign * std::acos(c / radius) + turn * 2 * std::acos(-1.)) / k;
        if (y < s.VMinimum || y > s.VMaximum)
          continue;
        const gp_Pnt p(x, y, .37);
        double       nearest = 1.e100;
        for (const auto& curve : result.Curves)
          for (size_t j = 1; j < curve.Points.size(); ++j)
          {
            const auto&  start = curve.Points[j - 1];
            const gp_Vec edge(start, curve.Points[j]);
            const double t =
              std::clamp(gp_Vec(start, p).Dot(edge) / edge.SquareMagnitude(), 0., 1.);
            nearest = std::min(nearest, p.Distance(start.Translated(edge * t)));
          }
        EXPECT_LE(nearest, result.GeometricBound + 1.e-12);
        ++checked;
      }
  }
  EXPECT_GT(checked, 100);
  auto saddle         = settings();
  saddle.Frame        = gp_Ax2(gp_Pnt(0, 0, 1.25), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
  saddle.MaximumCells = 64;
  const auto critical = BRep_ImplicitPlaneSection::Trace(gyroid(std::sqrt(2.)), saddle);
  EXPECT_FALSE(critical.Resolved);
  EXPECT_GT(critical.PotentialCriticalCells, 0);
}
