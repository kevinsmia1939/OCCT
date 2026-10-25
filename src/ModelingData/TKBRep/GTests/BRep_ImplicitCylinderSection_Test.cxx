// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRep_ImplicitCylinderSection.hxx>
#include <BRep_TImplicitCylinderPatch.hxx>
#include <BRep_ImplicitIO.hxx>
#include <BRep_ImplicitSphereOps.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <Geom_Circle.hxx>
#include <TopExp.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_ImplicitSphere.hxx>
#include <Standard_NotImplemented.hxx>
#include <TopoDS_Wire.hxx>
#include <Standard_AbortiveTransaction.hxx>
#include <Standard_DomainError.hxx>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

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
int analyticDiskWitnesses(const Section::Settings& s, int cap, double iso,
                          const Section::Result& disk)
{
  const double k = 2 * std::acos(-1.) / 10;
  const double z = s.Frame.Location().Z() + (cap == 0 ? s.HeightMinimum : s.HeightMaximum);
  int witnesses = 0;
  for (int i = 0; i <= 100; ++i)
  {
    const double localY = s.Radius * (2. * i / 100 - 1);
    const double y = s.Frame.Location().Y() + localY;
    const double a = std::cos(k * y), b = std::sin(k * z);
    const double c = std::sin(k * y) * std::cos(k * z) - iso;
    const double radius = std::hypot(a, b);
    if (std::abs(c) > radius || radius == 0) continue;
    const double phase = std::atan2(b, a), angle = std::asin(-c / radius);
    for (double branch : {angle, std::acos(-1.) - angle})
      for (int turn = -2; turn <= 2; ++turn)
      {
        const double x = (branch - phase + turn * 2 * std::acos(-1.)) / k;
        if (std::hypot(x - s.Frame.Location().X(), localY) > s.Radius) continue;
        ++witnesses;
        EXPECT_LE(distance(gp_Pnt(x, y, z), disk), disk.GeometricBound + 1.e-10);
      }
  }
  return witnesses;
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

TEST(BRepImplicitCylinderSection, EndCapRootsRetainOwnedVertexAndCurveIdentity)
{
  const auto face = field();
  auto s = settings();
  s.Radius = .5; s.AngleMinimum = 3;
  s.AngleMaximum = s.AngleMinimum + 2 * std::acos(-1.);
  s.HeightMinimum = 0; s.HeightMaximum = .5;
  const auto result = Section::TracePeriodicBoundary(face, s);
  ASSERT_TRUE(result.Section.Resolved);
  ASSERT_TRUE(result.BoundaryIncidenceVerified);
  ASSERT_EQ(result.Roots.size(), 4);
  ASSERT_EQ(result.CurveRoots.size(), result.Section.Curves.size());
  std::vector<int> uses(result.Roots.size(), 0);
  std::array<int, 2> caps {{0, 0}};
  for (size_t index = 0; index < result.Roots.size(); ++index)
  {
    const auto& root = result.Roots[index];
    ++caps[root.HeightBoundary];
    EXPECT_FALSE(root.Vertex.IsNull());
    const auto p = BRep_Tool::Pnt(root.Vertex);
    const double tolerance = BRep_Tool::Tolerance(root.Vertex);
    EXPECT_NEAR(std::hypot(p.X(), p.Y()), s.Radius, 1.e-13);
    EXPECT_EQ(p.Z(), root.HeightBoundary == 0 ? s.HeightMinimum : s.HeightMaximum);
    EXPECT_LE(p.Distance(root.First), tolerance);
    EXPECT_LE(p.Distance(root.Last), tolerance);
    EXPECT_GE(root.AngleLast, root.AngleFirst);
    EXPECT_LE(s.Radius * (root.AngleLast - root.AngleFirst), .1 / 32 + 1.e-14);
    EXPECT_LE(BRep_ImplicitTool::Value(face, root.First) * BRep_ImplicitTool::Value(face, root.Last), 0);
    if (index && root.HeightBoundary == result.Roots[index - 1].HeightBoundary)
      EXPECT_GT(root.AngleFirst, result.Roots[index - 1].AngleLast);
  }
  EXPECT_EQ(caps[0], 2); EXPECT_EQ(caps[1], 2);
  for (size_t index = 0; index < result.CurveRoots.size(); ++index)
  {
    const auto& curve = result.Section.Curves[index];
    ASSERT_FALSE(curve.Closed);
    for (int end = 0; end < 2; ++end)
    {
      const size_t id = result.CurveRoots[index][end];
      ASSERT_GT(id, 0); ASSERT_LE(id, result.Roots.size());
      ++uses[id - 1];
      EXPECT_EQ(BRep_Tool::Pnt(result.Roots[id - 1].Vertex).Distance(
        end == 0 ? curve.Points.front() : curve.Points.back()), 0);
    }
  }
  for (auto count : uses) EXPECT_EQ(count, 1);
  auto copy = result;
  EXPECT_TRUE(copy.Roots.front().Vertex.IsSame(result.Roots.front().Vertex));
  s.MaximumCells = 4;
  const auto failed = Section::TracePeriodicBoundary(face, s);
  EXPECT_FALSE(failed.Section.Resolved); EXPECT_FALSE(failed.BoundaryIncidenceVerified);
  EXPECT_TRUE(failed.Roots.empty()); EXPECT_TRUE(failed.CurveRoots.empty());
}

TEST(BRepImplicitCylinderSection, AnalyticRingArcsShareRootsAndPartitionSigns)
{
  const auto face = field();
  auto s = settings();
  s.Radius = .5; s.AngleMinimum = 3;
  s.AngleMaximum = s.AngleMinimum + 2 * std::acos(-1.);
  s.HeightMinimum = 0; s.HeightMaximum = .5;
  const auto result = Section::BuildBoundaryRings(face, s);
  ASSERT_TRUE(result.Verified);
  for (int cap = 0; cap < 2; ++cap)
  {
    ASSERT_EQ(result.Rings[cap].size(), 2);
    double total = 0;
    for (const auto& arc : result.Rings[cap])
    {
      EXPECT_TRUE(BRepCheck_Analyzer(arc.Edge).IsValid());
      EXPECT_EQ(arc.Edge.Orientation(), cap == 0 ? TopAbs_REVERSED : TopAbs_FORWARD);
      double first, last;
      const auto curve = BRep_Tool::Curve(arc.Edge, first, last);
      ASSERT_FALSE(occ::down_cast<Geom_Circle>(curve).IsNull());
      EXPECT_GT(last, first); total += last - first;
      TopoDS_Vertex a, b; TopExp::Vertices(arc.Edge, a, b);
      EXPECT_TRUE(a.IsSame(result.Boundary.Roots[arc.LowerRoot - 1].Vertex));
      EXPECT_TRUE(b.IsSame(result.Boundary.Roots[arc.UpperRoot - 1].Vertex));
      EXPECT_LE(curve->Value(first).Distance(BRep_Tool::Pnt(a)), BRep_Tool::Tolerance(a));
      EXPECT_LE(curve->Value(last).Distance(BRep_Tool::Pnt(b)), BRep_Tool::Tolerance(b));
      // Independent dense circle samples away from root uncertainty must agree.
      const auto& lower = result.Boundary.Roots[arc.LowerRoot - 1];
      const auto& upper = result.Boundary.Roots[arc.UpperRoot - 1];
      const double low = lower.AngleLast, high = upper.AngleFirst + (arc.UpperRoot < arc.LowerRoot ? 2 * std::acos(-1.) : 0);
      for (int i = 1; i < 30; ++i)
      {
        const double angle = low + (high - low) * i / 30;
        const auto point = curve->Value(angle);
        const double k = 2 * std::acos(-1.) / 10;
        const double value = std::sin(k * point.X()) * std::cos(k * point.Y())
          + std::sin(k * point.Y()) * std::cos(k * point.Z())
          + std::sin(k * point.Z()) * std::cos(k * point.X()) - .21;
        EXPECT_EQ(value < 0, arc.Negative);
      }
    }
    EXPECT_NEAR(total, 2 * std::acos(-1.), 1.e-13);
    EXPECT_NE(result.Rings[cap][0].Negative, result.Rings[cap][1].Negative);
  }
  gp_Trsf pose;
  pose.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), .4);
  pose.SetTranslationPart(gp_Vec(3, -4, 2));
  auto placedFace = face;
  placedFace.Location(TopLoc_Location(pose));
  auto placedSettings = s;
  placedSettings.Frame.Transform(pose);
  const auto placed = Section::BuildBoundaryRings(placedFace, placedSettings);
  ASSERT_TRUE(placed.Verified);
  for (int cap = 0; cap < 2; ++cap)
    for (size_t i = 0; i < result.Rings[cap].size(); ++i)
    {
      EXPECT_EQ(placed.Rings[cap][i].Negative, result.Rings[cap][i].Negative);
      EXPECT_TRUE(BRepCheck_Analyzer(placed.Rings[cap][i].Edge).IsValid());
    }
  s.MaximumCells = 4;
  const auto failed = Section::BuildBoundaryRings(face, s);
  EXPECT_FALSE(failed.Verified); EXPECT_TRUE(failed.Rings[0].empty()); EXPECT_TRUE(failed.Rings[1].empty());
}

TEST(BRepImplicitCylinderSection, RootFreeRingsAreActualClosedCirclePartitions)
{
  auto s = settings();
  s.Radius = .5; s.AngleMaximum = s.AngleMinimum + 2 * std::acos(-1.);
  const auto result = Section::BuildBoundaryRings(field(10), s);
  ASSERT_TRUE(result.Verified);
  for (const auto& ring : result.Rings)
  {
    ASSERT_EQ(ring.size(), 2);
    TopoDS_Vertex a, b, c, d;
    TopExp::Vertices(ring[0].Edge, a, b); TopExp::Vertices(ring[1].Edge, c, d);
    EXPECT_TRUE(a.IsSame(d)); EXPECT_TRUE(b.IsSame(c));
    for (const auto& arc : ring)
    {
      EXPECT_TRUE(arc.Negative); EXPECT_EQ(arc.LowerRoot, 0); EXPECT_EQ(arc.UpperRoot, 0);
      EXPECT_TRUE(BRepCheck_Analyzer(arc.Edge).IsValid());
    }
  }
  const auto positive = Section::BuildBoundaryRings(field(-10), s);
  ASSERT_TRUE(positive.Verified);
  for (const auto& ring : positive.Rings) for (const auto& arc : ring) EXPECT_FALSE(arc.Negative);
}

TEST(BRepImplicitCylinderSection, RingPrecisionCannotExceedRequestedTolerance)
{
  auto s = settings();
  s.Radius = .5; s.AngleMinimum = 3;
  s.AngleMaximum = s.AngleMinimum + 2 * std::acos(-1.);
  s.Frame.SetLocation(gp_Pnt(1000, 0, 0));
  s.DistanceTolerance = 2.e-9;
  auto face = field(10);
  gp_Trsf pose; pose.SetTranslation(gp_Vec(1000, 0, 0));
  face.Location(TopLoc_Location(pose));
  const auto result = Section::BuildBoundaryRings(face, s);
  ASSERT_TRUE(result.Boundary.Section.Resolved);
  EXPECT_FALSE(result.Verified);
  EXPECT_TRUE(result.Rings[0].empty()); EXPECT_TRUE(result.Rings[1].empty());
}

TEST(BRepImplicitCylinderSection, CylinderRestrictionQueriesPlacementCopyAndPersistence)
{
  const auto source = field();
  const gp_Ax2 frame(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
  const auto face = BRep_TImplicitCylinderPatch::Restrict(source, frame, 2, -1, 1);
  EXPECT_TRUE(BRep_Tool::IsImplicit(face));
  EXPECT_FALSE(BRep_Tool::HasParametricSurface(face));
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(face, gp_Pnt(2, 0, 0)));
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(face, gp_Pnt(0, 0, 1)));
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(face, gp_Pnt(1.5, 1.5, 0)));
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(face, gp_Pnt(0, 0, 1.01)));
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(source, gp_Pnt(1.5, 1.5, 0)));
  EXPECT_EQ(BRep_ImplicitTool::Value(face, gp_Pnt(.3, .2, .1)),
            BRep_ImplicitTool::Value(source, gp_Pnt(.3, .2, .1)));
  const auto data = occ::down_cast<BRep_TImplicitCylinderPatch>(face.TShape());
  const auto copy = occ::down_cast<BRep_TImplicitCylinderPatch>(data->EmptyCopy());
  ASSERT_FALSE(copy.IsNull()); EXPECT_EQ(copy->Radius(), 2);
  gp_Trsf pose; pose.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), .4);
  pose.SetTranslationPart(gp_Vec(3, -4, 2));
  auto placed = face; placed.Location(TopLoc_Location(pose));
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(placed, gp_Pnt(.4, .3, .2).Transformed(pose)));
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(placed, gp_Pnt(1.5, 1.5, 0).Transformed(pose)));
  std::stringstream stream;
  BRep_ImplicitIO::Write(placed, stream);
  EXPECT_EQ(stream.str().find("OCCT_IMPLICIT_FACE 3\n"), 0u);
  TopoDS_Face restored; BRep_ImplicitIO::Read(stream, restored);
  ASSERT_FALSE(occ::down_cast<BRep_TImplicitCylinderPatch>(restored.TShape()).IsNull());
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(restored, gp_Pnt(.4, .3, .2).Transformed(pose)));
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(restored, gp_Pnt(1.5, 1.5, 0).Transformed(pose)));
  std::stringstream primitive; BRep_ImplicitIO::Write(source, primitive);
  EXPECT_EQ(primitive.str().find("OCCT_IMPLICIT_FACE 1\n"), 0u);
  EXPECT_THROW(BRep_TImplicitCylinderPatch::Restrict(face, frame, 1, -1, 1), Standard_NotImplemented);
  EXPECT_THROW(BRep_TImplicitCylinderPatch::Restrict(source, frame, -1, -1, 1), Standard_DomainError);
  EXPECT_THROW(BRep_TImplicitCylinderPatch::Restrict(source, frame, 20, -1, 1), Standard_DomainError);
  // Invalid cylinder records leave a pre-existing destination unchanged.
  auto invalid = stream.str();
  const auto at = invalid.find("2 -1 1\n"); ASSERT_NE(at, std::string::npos);
  invalid.replace(at, 7, "-2 -1 1\n");
  std::istringstream bad(invalid);
  const auto before = restored.TShape();
  EXPECT_THROW(BRep_ImplicitIO::Read(bad, restored), Standard_DomainError);
  EXPECT_EQ(restored.TShape(), before);
  TopoDS_Face withChild; withChild.TShape(data->EmptyCopy());
  TopoDS_Wire wire; BRep_Builder().MakeWire(wire); BRep_Builder().Add(withChild, wire);
  std::stringstream rejected; rejected << "prefix";
  EXPECT_THROW(BRep_ImplicitIO::Write(withChild, rejected), Standard_NotImplemented);
  EXPECT_EQ(rejected.str(), "prefix");
  const gp_Ax2 oblique(gp_Pnt(.3, .2, .1), gp_Dir(1, 2, 3), gp_Dir(1, 0, 0));
  const auto tilted = BRep_TImplicitCylinderPatch::Restrict(source, oblique, 1, -.7, .9);
  const auto inside = oblique.Location().Translated(gp_Vec(oblique.Direction()) * .2
                                                    + gp_Vec(oblique.XDirection()) * .3);
  const auto outside = oblique.Location().Translated(gp_Vec(oblique.XDirection()) * 1.2);
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(tilted, inside));
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(tilted, outside));
  std::stringstream obliqueRecord; BRep_ImplicitIO::Write(tilted, obliqueRecord);
  TopoDS_Face tiltedRestored; BRep_ImplicitIO::Read(obliqueRecord, tiltedRestored);
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(tiltedRestored, inside));
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(tiltedRestored, outside));
}

TEST(BRepImplicitCylinderSection, CylinderRestrictionGuardsProjectionRaysAndGeometryLocation)
{
  TopoDS_Face sphere;
  BRep_Builder().MakeImplicitFace(sphere, new Geom_ImplicitSphere(gp_Pnt(0, 0, 0), 1),
    gp_Pnt(-3, -3, -3), gp_Pnt(3, 3, 3), TopLoc_Location(), .1);
  const gp_Ax2 frame(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
  const auto clipped = BRep_TImplicitCylinderPatch::Restrict(sphere, frame, .5, -2, 2);
  EXPECT_THROW(BRep_ImplicitSphereOps::Tessellate(clipped, .1), Standard_NotImplemented);
  EXPECT_TRUE(BRep_ImplicitTool::Project(clipped, gp_Pnt(0, 0, 1.04)).Converged);
  EXPECT_FALSE(BRep_ImplicitTool::Project(clipped, gp_Pnt(1.04, 0, 0)).Converged);
  EXPECT_TRUE(BRep_ImplicitTool::IntersectRay(clipped, gp_Pnt(-3, 0, 0), gp_Vec(1, 0, 0), 0, 6, 256, .001).empty());
  EXPECT_EQ(BRep_ImplicitTool::IntersectRay(clipped, gp_Pnt(0, 0, -3), gp_Vec(0, 0, 1), 0, 6, 256, .001).size(), 2);
  gp_Trsf pose; pose.SetTranslation(gp_Vec(3, -4, 2));
  TopoDS_Face geometryPlaced;
  BRep_Builder().MakeImplicitFace(geometryPlaced, new Geom_ImplicitGyroid(10, .21),
    gp_Pnt(-10, -10, -10), gp_Pnt(10, 10, 10), TopLoc_Location(pose), .1);
  const auto restricted = BRep_TImplicitCylinderPatch::Restrict(geometryPlaced, frame, 2, -1, 1);
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(restricted, gp_Pnt(.4, .3, .2).Transformed(pose)));
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(restricted, gp_Pnt(1.5, 1.5, 0).Transformed(pose)));
  std::stringstream record; BRep_ImplicitIO::Write(restricted, record);
  TopoDS_Face restored; BRep_ImplicitIO::Read(record, restored);
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(restored, gp_Pnt(.4, .3, .2).Transformed(pose)));
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(restored, gp_Pnt(1.5, 1.5, 0).Transformed(pose)));
}

TEST(BRepImplicitCylinderSection, AcceptedPlaneGraphRetainsCellAndSharedRootIdentities)
{
  BRep_ImplicitPlaneSection::Settings s;
  s.UMinimum = s.VMinimum = -.5; s.UMaximum = s.VMaximum = .5;
  s.DistanceTolerance = .05;
  const auto graph = BRep_ImplicitPlaneSection::TraceGraph(field(), s);
  ASSERT_TRUE(graph.Section.Resolved);
  ASSERT_FALSE(graph.Arcs.empty());
  EXPECT_LE(graph.Section.GeometricBound, .05);
  const double step = 1. / graph.Section.Resolution;
  for (const auto& arc : graph.Arcs)
    for (auto vertex : {arc.First, arc.Last})
    {
      ASSERT_LT(vertex, graph.Vertices.size());
      const auto& p = graph.Vertices[vertex];
      EXPECT_GE(p.X(), -.5 + step * arc.U - 1.e-12);
      EXPECT_LE(p.X(), -.5 + step * (arc.U + 1) + 1.e-12);
      EXPECT_GE(p.Y(), -.5 + step * arc.V - 1.e-12);
      EXPECT_LE(p.Y(), -.5 + step * (arc.V + 1) + 1.e-12);
    }
  s.MaximumCells = 4;
  const auto failed = BRep_ImplicitPlaneSection::TraceGraph(field(), s);
  EXPECT_FALSE(failed.Section.Resolved); EXPECT_TRUE(failed.Vertices.empty()); EXPECT_TRUE(failed.Arcs.empty());
}

TEST(BRepImplicitCylinderSection, EndCapDisksAttachAllCircleRootsByCellIdentity)
{
  auto s = settings();
  s.Radius = .5; s.AngleMinimum = 3; s.AngleMaximum = s.AngleMinimum + 2 * std::acos(-1.);
  s.HeightMinimum = 0; s.HeightMaximum = .5;
  const auto result = Section::TraceEndCapDisks(field(), s);
  ASSERT_TRUE(result.Verified) << result.Disks[0].AmbiguousCells << "/" << result.Disks[0].GraphDefects
    << "; " << result.Disks[1].AmbiguousCells << "/" << result.Disks[1].GraphDefects;
  ASSERT_EQ(result.Boundary.Roots.size(), 4);
  std::vector<int> uses(4, 0);
  for (int cap = 0; cap < 2; ++cap)
  {
    const auto& disk = result.Disks[cap];
    EXPECT_TRUE(disk.Resolved); EXPECT_LE(disk.GeometricBound, .1);
    EXPECT_GT(analyticDiskWitnesses(s, cap, .21, disk), 10);
    ASSERT_FALSE(disk.Curves.empty());
    ASSERT_EQ(disk.Curves.size(), result.CurveRoots[cap].size());
    for (size_t i = 0; i < disk.Curves.size(); ++i)
    {
      const auto& curve = disk.Curves[i];
      EXPECT_FALSE(curve.Closed);
      for (int end = 0; end < 2; ++end)
      {
        const auto id = result.CurveRoots[cap][i][end];
        ASSERT_GT(id, 0); ASSERT_LE(id, 4); ++uses[id - 1];
        EXPECT_EQ(result.Boundary.Roots[id - 1].HeightBoundary, cap);
        EXPECT_EQ(BRep_Tool::Pnt(result.Boundary.Roots[id - 1].Vertex).Distance(
          end == 0 ? curve.Points.front() : curve.Points.back()), 0);
      }
      for (const auto& p : curve.Points)
      { EXPECT_LE(std::hypot(p.X(), p.Y()), .5 + 1.e-12); EXPECT_LE(std::abs(BRep_ImplicitTool::Value(field(), p)), .01); }
    }
  }
  for (auto count : uses) EXPECT_EQ(count, 1);
  s.MaximumCells = 4;
  const auto failed = Section::TraceEndCapDisks(field(), s);
  EXPECT_FALSE(failed.Verified);
  for (int cap = 0; cap < 2; ++cap)
  { EXPECT_FALSE(failed.Disks[cap].Resolved); EXPECT_TRUE(failed.Disks[cap].Curves.empty()); EXPECT_TRUE(failed.CurveRoots[cap].empty()); }
}

TEST(BRepImplicitCylinderSection, EndCapDisksRetainClosedInteriorIslandsWithoutCircleRoots)
{
  auto s = settings();
  s.Frame.SetLocation(gp_Pnt(2.5, 1.25, 0)); s.Radius = 1.5;
  s.AngleMinimum = 3; s.AngleMaximum = s.AngleMinimum + 2 * std::acos(-1.);
  s.HeightMinimum = 0; s.HeightMaximum = .02;
  const auto result = Section::TraceEndCapDisks(field(1.3), s);
  ASSERT_TRUE(result.Verified) << result.Boundary.Section.UnresolvedCells << "; "
    << result.Disks[0].UnresolvedCells << "/" << result.Disks[0].GraphDefects;
  EXPECT_TRUE(result.Boundary.Roots.empty());
  for (int cap = 0; cap < 2; ++cap)
  {
    ASSERT_EQ(result.Disks[cap].Curves.size(), 1);
    EXPECT_TRUE(result.Disks[cap].Curves[0].Closed);
    EXPECT_EQ(result.CurveRoots[cap][0][0], 0); EXPECT_EQ(result.CurveRoots[cap][0][1], 0);
    EXPECT_GT(result.Disks[cap].Curves[0].Points.size(), 10);
    EXPECT_GT(analyticDiskWitnesses(s, cap, 1.3, result.Disks[cap]), 10);
  }
}

TEST(BRepImplicitCylinderSection, EndCapDisksHandleInteriorZeroEmptyPlacementAndCancellation)
{
  auto s = settings();
  s.Radius = .5; s.AngleMinimum = 3; s.AngleMaximum = s.AngleMinimum + 2 * std::acos(-1.);
  s.HeightMinimum = 0; s.HeightMaximum = .5;
  const auto zero = Section::TraceEndCapDisks(field(0), s);
  ASSERT_TRUE(zero.Verified) << zero.Boundary.Section.UnresolvedCells << "; "
    << zero.Disks[0].UnresolvedCells << "/" << zero.Disks[0].GraphDefects;
  EXPECT_GT(zero.Disks[0].ExactLatticeNodeCount, 0);
  EXPECT_GT(analyticDiskWitnesses(s, 0, 0, zero.Disks[0]), 10);
  const auto empty = Section::TraceEndCapDisks(field(10), s);
  EXPECT_TRUE(empty.Verified);
  for (const auto& disk : empty.Disks) EXPECT_TRUE(disk.Curves.empty());
  gp_Trsf pose; pose.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), .4);
  pose.SetTranslationPart(gp_Vec(3, -4, 2));
  auto placedFace = field(); placedFace.Location(TopLoc_Location(pose));
  auto placedSettings = s; placedSettings.Frame.Transform(pose);
  const auto placed = Section::TraceEndCapDisks(placedFace, placedSettings);
  EXPECT_TRUE(placed.Verified);
  s.HeightMaximum = .37;
  const auto nonlattice = Section::TraceEndCapDisks(field(std::sin(2 * std::acos(-1.) / 10 * .37)), s);
  EXPECT_TRUE(nonlattice.Verified);
  auto criticalSettings = s;
  criticalSettings.Frame.SetLocation(gp_Pnt(2.5, 1.25, 0));
  criticalSettings.Radius = 1.5;
  criticalSettings.HeightMinimum = -.5; criticalSettings.HeightMaximum = 0;
  const auto lateFailure = Section::TraceEndCapDisks(field(std::sqrt(2.)), criticalSettings);
  EXPECT_TRUE(lateFailure.Boundary.Section.Resolved);
  EXPECT_FALSE(lateFailure.Verified);
  for (const auto& disk : lateFailure.Disks)
  { EXPECT_FALSE(disk.Resolved); EXPECT_TRUE(disk.Curves.empty()); }
  s.WorkCheckpoint = [] { throw Standard_AbortiveTransaction("Disk query cancelled"); };
  EXPECT_THROW(Section::TraceEndCapDisks(field(), s), Standard_AbortiveTransaction);
}
