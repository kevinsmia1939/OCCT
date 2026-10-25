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
#include <Geom2d_Curve.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Wire.hxx>
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

TEST(BRepClass3dImplicitCylinder, SeamedCapsReuseSideChartEdgesAndRootHandles)
{
  const auto s = settings();
  for (double iso : {.21, 10.})
    for (bool negative : {true, false})
    {
      const auto result = Ops::BuildChartCaps(field(iso), s, negative);
      ASSERT_TRUE(result.Chart.Verified);
      ASSERT_TRUE(result.Regions.TopologyValidated) << result.Regions.Reason;
      const auto& graphRoots = result.Chart.Graph.Boundary.Roots;
      const auto& capRoots = result.Regions.Source.Disks.Boundary.Roots;
      ASSERT_EQ(graphRoots.size(), capRoots.size());
      for (size_t i = 0; i < graphRoots.size(); ++i)
        EXPECT_TRUE(graphRoots[i].Vertex.IsSame(capRoots[i].Vertex));
      for (int cap = 0; cap < 2; ++cap)
      {
        for (const auto& arc : result.Chart.Rings.Rings[cap])
        {
          int uses = 0;
          for (const auto& face : result.Regions.Faces[cap])
          {
            EXPECT_TRUE(BRepCheck_Analyzer(face).IsValid());
            for (TopExp_Explorer edge(face, TopAbs_EDGE); edge.More(); edge.Next())
              if (edge.Current().IsSame(arc.Edge)) ++uses;
          }
          EXPECT_EQ(uses, arc.Negative == negative ? 1 : 0);
        }
      }
    }
  const auto a = Ops::BuildChartCaps(field(), s, true);
  const auto b = Ops::BuildChartCaps(field(), s, false);
  for (int cap = 0; cap < 2; ++cap)
    EXPECT_NEAR(area(a.Regions.Faces[cap]) + area(b.Regions.Faces[cap]),
                std::acos(-1.) * s.Radius * s.Radius, 1.e-6);
  const auto failed = Ops::BuildChartCaps(field(), s, true, 1);
  EXPECT_FALSE(failed.Regions.TopologyValidated);
  EXPECT_TRUE(failed.Regions.Faces[0].empty()); EXPECT_TRUE(failed.Regions.Faces[1].empty());
}

TEST(BRepClass3dImplicitCylinder, ChartSeamBuildsValidOutwardCylinderFace)
{
  const auto s = settings();
  const auto result = Ops::BuildChartCaps(field(10), s);
  ASSERT_TRUE(result.Regions.TopologyValidated);
  const auto& chart = result.Chart;
  ASSERT_EQ(chart.SeamEdges.size(), 1);
  BRep_Builder builder;
  TopoDS_Face face; builder.MakeFace(face, chart.Surface, 1.e-7);
  TopoDS_Wire wire; builder.MakeWire(wire);
  for (const auto& arc : chart.Rings.Rings[0])
  { auto edge = arc.Edge; edge.Reverse(); builder.Add(wire, edge); }
  builder.Add(wire, chart.SeamEdges.front());
  for (auto arc = chart.Rings.Rings[1].rbegin(); arc != chart.Rings.Rings[1].rend(); ++arc)
  { auto edge = arc->Edge; edge.Reverse(); builder.Add(wire, edge); }
  auto seam = chart.SeamEdges.front(); seam.Reverse(); builder.Add(wire, seam);
  wire.Closed(true); builder.Add(face, wire);
  EXPECT_TRUE(BRepCheck_Analyzer(face).IsValid());
  EXPECT_NEAR(area({face}), 2 * std::acos(-1.) * s.Radius
                          * (s.HeightMaximum - s.HeightMinimum), 1.e-8);
  EXPECT_EQ(BRepClass_FaceClassifier(face, gp_Pnt(s.Radius, 0, .25), 1.e-7).State(), TopAbs_IN);
  EXPECT_EQ(BRepClass_FaceClassifier(face, gp_Pnt(s.Radius, 0, 1), 1.e-7).State(), TopAbs_OUT);
}

TEST(BRepClass3dImplicitCylinder, SideFieldEdgesShareChartVerticesAndBoundCylinderDeviation)
{
  for (bool closed : {false, true})
    for (bool negative : {false, true})
    {
      auto s = settings();
      if (closed)
      { s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
        s.HeightMinimum = -1.7; s.HeightMaximum = 1.3; }
      const auto result = Ops::BuildSideBoundary(field(), s, negative);
      ASSERT_TRUE(result.Verified) << result.Reason;
      ASSERT_FALSE(result.FieldEdges.empty());
      EXPECT_LE(result.SurfaceDeviationBound, s.DistanceTolerance);
      const auto& chart = result.Source.Chart;
      for (const auto& edge : result.FieldEdges)
      {
        double first, last, pf, pl;
        const auto curve = BRep_Tool::Curve(edge, first, last);
        const auto pc = BRep_Tool::CurveOnSurface(edge, chart.Surface, TopLoc_Location(), pf, pl);
        ASSERT_FALSE(pc.IsNull()); EXPECT_EQ(first, pf); EXPECT_EQ(last, pl);
        for (int i = 0; i <= 100; ++i)
        {
          const double t = first + (last - first) * i / 100;
          const auto uv = pc->Value(t);
          EXPECT_GE(uv.X(), s.AngleMinimum - 1.e-12);
          EXPECT_LE(uv.X(), s.AngleMaximum + 1.e-12);
          EXPECT_LE(curve->Value(t).Distance(chart.Surface->Value(uv.X(), uv.Y())),
                    result.SurfaceDeviationBound);
        }
        TopoDS_Vertex a, b; TopExp::Vertices(edge, a, b);
        for (const auto& v : {a, b})
        {
          bool shared = false;
          for (const auto& node : chart.Graph.Vertices)
            if (!node.Vertex.IsNull() && node.Vertex.IsSame(v)) shared = true;
          EXPECT_TRUE(shared);
        }
      }
    }
  auto s = settings(); s.MaximumCells = 4;
  const auto failed = Ops::BuildSideBoundary(field(), s);
  EXPECT_FALSE(failed.Verified); EXPECT_TRUE(failed.FieldEdges.empty());
}

TEST(BRepClass3dImplicitCylinder, NativeSideRegionsPartitionCylinderAndClassifyMaterial)
{
  for (bool winding : {false, true})
  {
    auto s = settings();
    if (winding)
    { s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
      s.HeightMinimum = -1.7; s.HeightMaximum = 1.3; }
    const auto source = field();
    const auto negative = Ops::BuildSides(source, s, true);
    const auto positive = Ops::BuildSides(source, s, false);
    ASSERT_TRUE(negative.TopologyValidated) << negative.Reason;
    ASSERT_TRUE(positive.TopologyValidated) << positive.Reason;
    EXPECT_NEAR(area(negative.Faces) + area(positive.Faces),
                2 * std::acos(-1.) * s.Radius * (s.HeightMaximum - s.HeightMinimum), 1.e-6);
    for (int i = 0; i < 31; ++i) for (int j = 0; j < 23; ++j)
    {
      const double u = s.AngleMinimum + (s.AngleMaximum - s.AngleMinimum) * (i + .37) / 31;
      const double v = s.HeightMinimum + (s.HeightMaximum - s.HeightMinimum) * (j + .43) / 23;
      const auto p = negative.Boundary.Source.Chart.Surface->Value(u, v);
      const double value = BRep_ImplicitTool::Value(source, p);
      if (std::abs(value) < .25) continue;
      for (const auto* result : {&negative, &positive})
      {
        int inside = 0;
        for (const auto& face : result->Faces)
        {
          EXPECT_TRUE(BRepCheck_Analyzer(face).IsValid());
          const auto state = BRepClass_FaceClassifier(face, gp_Pnt2d(u, v), 1.e-7).State();
          EXPECT_NE(state, TopAbs_ON); inside += state == TopAbs_IN;
        }
        EXPECT_EQ(inside, (value < 0) == result->Boundary.Source.Regions.Negative ? 1 : 0);
      }
    }
  }
  const auto full = Ops::BuildSides(field(10), settings());
  ASSERT_TRUE(full.TopologyValidated) << full.Reason;
  ASSERT_EQ(full.Faces.size(), 1);
  const auto empty = Ops::BuildSides(field(10), settings(), false);
  EXPECT_TRUE(empty.TopologyValidated); EXPECT_TRUE(empty.Faces.empty());
}

TEST(BRepClass3dImplicitCylinder, NativeSidesPlacementCancellationAndFailureCleanup)
{
  auto s = settings();
  const auto source = field();
  std::stringstream before; BRep_ImplicitIO::Write(source, before);
  const auto original = Ops::BuildSides(source, s);
  ASSERT_TRUE(original.TopologyValidated) << original.Reason;
  gp_Trsf move; move.SetTranslation(gp_Vec(1.2, -.4, 2.1));
  auto placed = source; placed.Location(TopLoc_Location(move));
  s.Frame.Transform(move);
  const auto moved = Ops::BuildSides(placed, s);
  ASSERT_TRUE(moved.TopologyValidated) << moved.Reason;
  EXPECT_NEAR(area(original.Faces), area(moved.Faces), 1.e-6);
  std::stringstream after; BRep_ImplicitIO::Write(source, after);
  EXPECT_EQ(before.str(), after.str());
  s = settings(); s.MaximumCells = 4;
  const auto budget = Ops::BuildSides(source, s);
  EXPECT_FALSE(budget.TopologyValidated); EXPECT_TRUE(budget.Faces.empty());
  s = settings();
  static size_t calls = 0;
  calls = 0; s.WorkCheckpoint = [] { ++calls; };
  ASSERT_TRUE(Ops::BuildSides(source, s).TopologyValidated);
  ASSERT_GT(calls, 2);
  --calls;
  s.WorkCheckpoint = [] { if (--calls == 0) throw Standard_AbortiveTransaction("late side cancellation"); };
  EXPECT_THROW(Ops::BuildSides(source, s), Standard_AbortiveTransaction);
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
