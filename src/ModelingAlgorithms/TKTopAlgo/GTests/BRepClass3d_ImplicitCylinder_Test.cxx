// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitCylinder.hxx>
#include <BRepClass3d_ImplicitCylinderIO.hxx>
#include <BRep_ImplicitIO.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepGProp.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <GProp_GProps.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Plane.hxx>
#include <Geom2d_Curve.hxx>
#include <BRep_TImplicitCylinderPatch.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS.hxx>
#include <Standard_AbortiveTransaction.hxx>
#include <Standard_DomainError.hxx>
#include <gtest/gtest.h>
#include <cmath>
#include <sstream>
#include <limits>

namespace
{
using Ops = BRepClass3d_ImplicitCylinder;
using Settings = BRep_ImplicitCylinderSection::Settings;
TopoDS_Face field(double iso = .21, double tolerance = .1)
{
  TopoDS_Face face;
  BRep_Builder().MakeImplicitFace(face, new Geom_ImplicitGyroid(10, iso),
    gp_Pnt(-10, -10, -10), gp_Pnt(10, 10, 10), TopLoc_Location(), tolerance);
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

TEST(BRepClass3dImplicitCylinder, TrimWiresCloseBySharedIdentityWithOppositeFaceIncidence)
{
  for (int fixture = 0; fixture < 4; ++fixture)
    for (bool negative : {false, true})
    {
      auto s = settings();
      double iso = .21;
      if (fixture == 1)
      { s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
        s.HeightMinimum = -1.7; s.HeightMaximum = 1.3; }
      if (fixture == 2) iso = 10;
      if (fixture == 3)
      { iso = 1.3; s.Frame.SetLocation(gp_Pnt(2.5, 1.25, 0)); s.Radius = 1.5;
        s.HeightMinimum = 0; s.HeightMaximum = .02; }
      const auto result = Ops::BuildTrimBoundary(field(iso), s, negative);
      ASSERT_TRUE(result.IncidenceVerified) << fixture << ": " << result.Reason;
      std::vector<TopoDS_Edge> edges = result.Source.Boundary.FieldEdges;
      for (const auto& cap : result.Source.Boundary.Source.Regions.FieldEdges)
        edges.insert(edges.end(), cap.begin(), cap.end());
      ASSERT_EQ(result.EdgeWires.size(), edges.size());
      if (fixture == 2) EXPECT_TRUE(result.Wires.empty());
      else EXPECT_FALSE(result.Wires.empty());
      if (fixture == 3) EXPECT_EQ(result.Wires.size(), 2);
      for (size_t i = 0; i < result.Wires.size(); ++i)
      {
        const auto& wire = result.Wires[i];
        EXPECT_TRUE(wire.Closed()); EXPECT_TRUE(BRepCheck_Analyzer(wire).IsValid());
        size_t count = 0;
        for (TopExp_Explorer it(wire, TopAbs_EDGE); it.More(); it.Next())
        {
          int match = -1;
          for (size_t e = 0; e < edges.size(); ++e) if (it.Current().IsSame(edges[e])) match = int(e);
          ASSERT_GE(match, 0);
          EXPECT_EQ(result.EdgeWires[size_t(match)], int(i));
          EXPECT_NE(it.Current().Orientation(), edges[size_t(match)].Orientation()); ++count;
        }
        EXPECT_GT(count, 0);
      }
    }
  auto s = settings(); s.MaximumCells = 4;
  const auto failed = Ops::BuildTrimBoundary(field(), s);
  EXPECT_FALSE(failed.IncidenceVerified); EXPECT_TRUE(failed.Wires.empty());
  EXPECT_TRUE(failed.EdgeWires.empty());
  const auto source = field();
  std::stringstream before; BRep_ImplicitIO::Write(source, before);
  s = settings();
  static size_t calls = 0;
  calls = 0; s.WorkCheckpoint = [] { ++calls; };
  ASSERT_TRUE(Ops::BuildTrimBoundary(source, s).IncidenceVerified);
  ASSERT_GT(calls, 2); --calls;
  s.WorkCheckpoint = [] { if (--calls == 0) throw Standard_AbortiveTransaction("late trim cancellation"); };
  EXPECT_THROW(Ops::BuildTrimBoundary(source, s), Standard_AbortiveTransaction);
  std::stringstream after; BRep_ImplicitIO::Write(source, after);
  EXPECT_EQ(before.str(), after.str());
}

TEST(BRepClass3dImplicitCylinder, AxialSheetConesBoundEveryFibreAndRefuseUnprovedCases)
{
  const auto source = field();
  for (bool reverse : {false, true})
  {
    auto s = settings();
    if (reverse)
    { s.Frame = gp_Ax2(gp_Pnt(), gp_Dir(0, 0, -1), gp_Dir(1, 0, 0));
      s.HeightMinimum = -.5; s.HeightMaximum = 0; }
    const auto result = Ops::CheckAxialSheet(source, s);
    ASSERT_TRUE(result.MonotonicityVerified) << result.Reason;
    EXPECT_EQ(result.Direction, reverse ? -1 : 1);
    EXPECT_GT(result.DerivativeLowerBound, 0); EXPECT_GT(result.Cells, 1);
    for (int i = 0; i <= 12; ++i) for (int j = 0; j <= 12; ++j) for (int h = 0; h <= 8; ++h)
    {
      const double x = s.Radius * (2. * i / 12 - 1), y = s.Radius * (2. * j / 12 - 1);
      if (std::hypot(x, y) > s.Radius) continue;
      const double z = s.HeightMinimum + (s.HeightMaximum - s.HeightMinimum) * h / 8;
      const auto p = s.Frame.Location().Translated(gp_Vec(s.Frame.XDirection()) * x
        + gp_Vec(s.Frame.YDirection()) * y + gp_Vec(s.Frame.Direction()) * z);
      EXPECT_GE(result.Direction * BRep_ImplicitTool::Gradient(source, p).Dot(gp_Vec(s.Frame.Direction())),
                result.DerivativeLowerBound);
    }
  }
  auto s = settings();
  const auto budget = Ops::CheckAxialSheet(source, s, true, 2000000, 24, 1);
  EXPECT_FALSE(budget.MonotonicityVerified); EXPECT_EQ(budget.Cells, 1);
  EXPECT_EQ(budget.Direction, 0); EXPECT_EQ(budget.DerivativeLowerBound, 0);
  EXPECT_FALSE(Ops::CheckAxialSheet(source, s, true, 2000000, 0).MonotonicityVerified);
  EXPECT_THROW(Ops::CheckAxialSheet(source, s, true, 2000000, 24, 0), Standard_DomainError);
  gp_Trsf move; move.SetTranslation(gp_Vec(1.2, -.4, 2.1));
  auto placed = source; placed.Location(TopLoc_Location(move)); s.Frame.Transform(move);
  const auto moved = Ops::CheckAxialSheet(placed, s);
  ASSERT_TRUE(moved.MonotonicityVerified) << moved.Reason;
  EXPECT_EQ(moved.Direction, 1);
  const auto empty = Ops::CheckAxialSheet(field(10), settings());
  ASSERT_TRUE(empty.MonotonicityVerified);
  EXPECT_TRUE(empty.Boundary.Wires.empty()); // Uniqueness is not an existence claim.
  s = settings();
  s.Frame.SetLocation(gp_Pnt(2.5, 1.25, 0)); s.Radius = 1.5;
  s.HeightMinimum = 0; s.HeightMaximum = .02;
  const auto folded = Ops::CheckAxialSheet(field(1.3), s, true, 2000000, 12, 10000);
  ASSERT_TRUE(folded.Boundary.IncidenceVerified) << folded.Boundary.Reason;
  EXPECT_FALSE(folded.MonotonicityVerified);
  EXPECT_EQ(folded.Direction, 0); EXPECT_EQ(folded.DerivativeLowerBound, 0);
  s = settings();
  static size_t calls = 0;
  calls = 0; s.WorkCheckpoint = [] { ++calls; };
  ASSERT_TRUE(Ops::CheckAxialSheet(source, s).MonotonicityVerified);
  ASSERT_GT(calls, 2); --calls;
  s.WorkCheckpoint = [] { if (--calls == 0) throw Standard_AbortiveTransaction("late axial cancellation"); };
  EXPECT_THROW(Ops::CheckAxialSheet(source, s), Standard_AbortiveTransaction);
}

TEST(BRepClass3dImplicitCylinder, ProjectedRegionsRetainNativeWireMappingWithoutPromotingFieldProof)
{
  for (bool winding : {false, true})
  {
    auto s = settings();
    if (winding)
    { s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
      s.HeightMinimum = -1.7; s.HeightMaximum = 1.3; }
    const auto source = field();
    const auto a = Ops::BuildProjectedRegions(source, s, true);
    const auto b = Ops::BuildProjectedRegions(source, s, false);
    ASSERT_TRUE(a.ProxyTopologyValidated) << a.Reason;
    ASSERT_TRUE(b.ProxyTopologyValidated) << b.Reason;
    EXPECT_FALSE(a.MathematicalCorrespondenceVerified);
    EXPECT_FALSE(b.MathematicalCorrespondenceVerified);
    ASSERT_FALSE(a.Faces.empty());
    EXPECT_NEAR(area(a.Faces), area(b.Faces), 1.e-6);
    ASSERT_EQ(a.WireRegions.size(), a.Source.Boundary.Wires.size());
    for (int id : a.WireRegions) { EXPECT_GE(id, 0); EXPECT_LT(size_t(id), a.Faces.size()); }
    for (const auto& face : a.Faces) EXPECT_TRUE(BRepCheck_Analyzer(face).IsValid());
    if (winding) EXPECT_NEAR(area(a.Faces), std::acos(-1.) * .25, 1.e-8);
    for (int i = 0; i <= 20; ++i) for (int j = 0; j <= 20; ++j)
    {
      const double x = s.Radius * (2. * i / 20 - 1), y = s.Radius * (2. * j / 20 - 1);
      if (std::hypot(x, y) >= s.Radius - .05) continue;
      const double lo = BRep_ImplicitTool::Value(source, gp_Pnt(x, y, s.HeightMinimum));
      const double hi = BRep_ImplicitTool::Value(source, gp_Pnt(x, y, s.HeightMaximum));
      if (std::min(std::abs(lo), std::abs(hi)) < .2) continue;
      int inside = 0;
      for (const auto& face : a.Faces)
        inside += BRepClass_FaceClassifier(face, gp_Pnt(x, y, 0), 1.e-7).State() == TopAbs_IN;
      EXPECT_EQ(inside, lo < 0 && hi > 0 ? 1 : 0);
    }
  }
  const auto empty = Ops::BuildProjectedRegions(field(10), settings());
  ASSERT_TRUE(empty.ProxyTopologyValidated) << empty.Reason;
  EXPECT_TRUE(empty.Faces.empty()); EXPECT_TRUE(empty.WireRegions.empty());
  auto s = settings();
  const auto failed = Ops::BuildProjectedRegions(field(), s, true, 2000000, 24, 1);
  EXPECT_FALSE(failed.ProxyTopologyValidated); EXPECT_TRUE(failed.Faces.empty());
  EXPECT_TRUE(failed.WireRegions.empty()); EXPECT_FALSE(failed.MathematicalCorrespondenceVerified);
}

TEST(BRepClass3dImplicitCylinder, FullDiskCorrespondenceAttachesActualCylinderRestrictedImplicitFace)
{
  auto s = settings(); s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
  s.HeightMinimum = -1.7; s.HeightMaximum = 1.3;
  const auto source = field();
  std::stringstream before; BRep_ImplicitIO::Write(source, before);
  for (bool negative : {false, true})
  {
    const auto result = Ops::BuildImplicitFaces(source, s, negative);
    ASSERT_TRUE(result.FacesAttached) << result.Reason;
    ASSERT_TRUE(result.Reference.MathematicalCorrespondenceVerified);
    EXPECT_GT(result.Reference.CorrespondenceCells, 0);
    ASSERT_EQ(result.Faces.size(), 1);
    const auto& face = result.Faces.front();
    EXPECT_FALSE(occ::down_cast<BRep_TImplicitCylinderPatch>(face.TShape()).IsNull());
    EXPECT_EQ(face.Orientation(), negative ? TopAbs_FORWARD : TopAbs_REVERSED);
    EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(face, gp_Pnt(0, 0, 0)));
    EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(face, gp_Pnt(.6, 0, 0)));
    EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(face, gp_Pnt(0, 0, 2)));
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
    {
      bool shared = false;
      for (TopExp_Explorer expected(result.Reference.Source.Boundary.Wires.front(), TopAbs_EDGE);
           expected.More(); expected.Next()) if (it.Current().IsEqual(expected.Current())) shared = true;
      EXPECT_TRUE(shared);
    }
  }
  gp_Trsf move; move.SetTranslation(gp_Vec(1.2, -.4, 2.1));
  auto placed = source; placed.Location(TopLoc_Location(move)); s.Frame.Transform(move);
  const auto translated = Ops::BuildImplicitFaces(placed, s);
  ASSERT_TRUE(translated.FacesAttached) << translated.Reason;
  EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(translated.Faces[0], gp_Pnt(1.2, -.4, 2.1)));
  EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(translated.Faces[0], gp_Pnt(1.8, -.4, 2.1)));
  const auto partial = Ops::BuildImplicitFaces(source, settings());
  EXPECT_FALSE(partial.FacesAttached); EXPECT_TRUE(partial.Faces.empty());
  EXPECT_TRUE(partial.Reference.ProxyTopologyValidated);
  EXPECT_FALSE(partial.Reference.MathematicalCorrespondenceVerified);
  auto reversed = settings();
  reversed.Frame = gp_Ax2(gp_Pnt(), gp_Dir(0, 0, -1), gp_Dir(1, 0, 0));
  reversed.HeightMinimum = -1.3; reversed.HeightMaximum = 1.7;
  reversed.AngleMinimum = .27; reversed.AngleMaximum = .27 + 2 * std::acos(-1.);
  const auto reverse = Ops::BuildImplicitFaces(source, reversed);
  ASSERT_TRUE(reverse.FacesAttached) << reverse.Reason;
  EXPECT_EQ(reverse.Reference.Source.Direction, -1);
  auto cancellable = settings(); cancellable.AngleMinimum = .27;
  cancellable.AngleMaximum = .27 + 2 * std::acos(-1.);
  cancellable.HeightMinimum = -1.7; cancellable.HeightMaximum = 1.3;
  static size_t calls = 0;
  calls = 0; cancellable.WorkCheckpoint = [] { ++calls; };
  ASSERT_TRUE(Ops::BuildImplicitFaces(source, cancellable).FacesAttached);
  ASSERT_GT(calls, 2); --calls;
  cancellable.WorkCheckpoint = [] { if (--calls == 0) throw Standard_AbortiveTransaction("late face cancellation"); };
  EXPECT_THROW(Ops::BuildImplicitFaces(source, cancellable), Standard_AbortiveTransaction);
  std::stringstream after; BRep_ImplicitIO::Write(source, after);
  EXPECT_EQ(before.str(), after.str());
}

TEST(BRepClass3dImplicitCylinder, NativeCylinderShellPairsSeamsAndHasSphericalVertexLinks)
{
  const auto source = field();
  for (bool reverse : {false, true}) for (bool negative : {false, true})
  {
    auto s = settings(); s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
    s.HeightMinimum = -1.7; s.HeightMaximum = 1.3;
    if (reverse)
    { s.Frame = gp_Ax2(gp_Pnt(), gp_Dir(0, 0, -1), gp_Dir(1, 0, 0));
      s.HeightMinimum = -1.3; s.HeightMaximum = 1.7; }
    const auto result = Ops::BuildNativeShell(source, s, negative);
    ASSERT_TRUE(result.TopologyValidated) << reverse << "/" << negative << ": " << result.Reason;
    EXPECT_FALSE(result.Shell.IsNull()); EXPECT_TRUE(result.Shell.Closed());
    EXPECT_EQ(result.Faces.size(), 3); EXPECT_EQ(result.Euler, 2);
    EXPECT_EQ(result.Vertices, 3); EXPECT_EQ(result.Edges, 4);
    int implicit = 0;
    for (TopExp_Explorer it(result.Shell, TopAbs_FACE); it.More(); it.Next())
      implicit += BRep_Tool::IsImplicit(TopoDS::Face(it.Current()));
    EXPECT_EQ(implicit, 1);
    const int retained = negative == !reverse ? 0 : 1;
    const auto& caps = result.Source.Reference.Source.Boundary.Source.Boundary.Source.Regions.Faces;
    EXPECT_EQ(caps[retained].size(), 1); EXPECT_TRUE(caps[1 - retained].empty());
  }
  auto s = settings(); s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
  s.HeightMinimum = -1.7; s.HeightMaximum = 1.3;
  gp_Trsf move; move.SetTranslation(gp_Vec(1.2, -.4, 2.1));
  auto placed = source; placed.Location(TopLoc_Location(move)); s.Frame.Transform(move);
  const auto translated = Ops::BuildNativeShell(placed, s);
  ASSERT_TRUE(translated.TopologyValidated) << translated.Reason;
  const auto partial = Ops::BuildNativeShell(source, settings());
  EXPECT_FALSE(partial.TopologyValidated); EXPECT_TRUE(partial.Shell.IsNull()); EXPECT_TRUE(partial.Faces.empty());
  static size_t calls = 0;
  calls = 0; s.WorkCheckpoint = [] { ++calls; };
  ASSERT_TRUE(Ops::BuildNativeShell(placed, s).TopologyValidated);
  ASSERT_GT(calls, 2); --calls;
  s.WorkCheckpoint = [] { if (--calls == 0) throw Standard_AbortiveTransaction("late shell cancellation"); };
  EXPECT_THROW(Ops::BuildNativeShell(placed, s), Standard_AbortiveTransaction);
}

TEST(BRepClass3dImplicitCylinder, ShellToleranceSeparatesSourcePolicyFromCombinedTrimBound)
{
  auto s = settings(); s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
  s.HeightMinimum = -1.7; s.HeightMaximum = 1.3; s.DistanceTolerance = .05;
  for (double tolerance : {1.e-6, .5}) for (bool negative : {false, true})
  {
    const auto source = field(.21, tolerance);
    std::stringstream before; BRep_ImplicitIO::Write(source, before);
    const auto result = Ops::CheckShellTolerance(source, s, negative);
    ASSERT_TRUE(result.ToleranceVerified) << result.Reason;
    EXPECT_GT(result.ContourBound, 0); EXPECT_GT(result.SurfaceDeviationBound, 0);
    EXPECT_EQ(result.FieldWireBound, result.ContourBound + result.SurfaceDeviationBound);
    EXPECT_LE(result.FieldWireBound, s.DistanceTolerance);
    EXPECT_EQ(result.MaximumFaceTolerance, s.DistanceTolerance);
    EXPECT_LE(result.MaximumEdgeTolerance, s.DistanceTolerance);
    EXPECT_LE(result.MaximumVertexTolerance, s.DistanceTolerance);
    EXPECT_EQ(BRep_Tool::Tolerance(result.Source.Source.Faces[0]), s.DistanceTolerance);
    EXPECT_EQ(BRep_Tool::Tolerance(source), tolerance);
    std::stringstream after; BRep_ImplicitIO::Write(source, after);
    EXPECT_EQ(before.str(), after.str());
  }
  auto small = s; small.MaximumCells = 4;
  const auto failed = Ops::CheckShellTolerance(field(), small);
  EXPECT_FALSE(failed.ToleranceVerified); EXPECT_FALSE(failed.Source.TopologyValidated);
  static size_t calls = 0;
  calls = 0; s.WorkCheckpoint = [] { ++calls; };
  ASSERT_TRUE(Ops::CheckShellTolerance(field(), s).ToleranceVerified);
  ASSERT_GT(calls, 2); --calls;
  s.WorkCheckpoint = [] { if (--calls == 0) throw Standard_AbortiveTransaction("late tolerance cancellation"); };
  EXPECT_THROW(Ops::CheckShellTolerance(field(), s), Standard_AbortiveTransaction);
}

TEST(BRepClass3dImplicitCylinder, ShellGeometryBoundsBothDirectionsOverTheCylinder)
{
  const auto source = field();
  for (double tolerance : {.1, .05}) for (bool negative : {false, true})
  {
    auto s = settings(); s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
    s.HeightMinimum = -1.7; s.HeightMaximum = 1.3; s.DistanceTolerance = tolerance;
    const auto result = Ops::CheckShellGeometry(source, s, negative);
    ASSERT_TRUE(result.BoundaryVerified) << result.Reason << " cells " << result.Cells << " tests " << result.SegmentTests;
    EXPECT_GT(result.Cells, 1); EXPECT_GT(result.SegmentTests, 1);
    EXPECT_LE(result.MathematicalToNativeBound, tolerance);
    EXPECT_LE(result.NativeToMathematicalBound, tolerance);
    EXPECT_LE(result.BoundaryBound, tolerance);
    EXPECT_GE(result.BoundaryBound, result.Source.FieldWireBound);
    const auto& sides = result.Source.Source.Source.Reference.Source.Boundary.Source;
    const double k = 2 * std::acos(-1.) / 10;
    for (int i = 0; i < 24; ++i)
    {
      const double angle = s.AngleMinimum + 2 * std::acos(-1.) * (i + .37) / 24;
      const double x = s.Radius * std::cos(angle), y = s.Radius * std::sin(angle);
      const double a = std::sin(k * y), b = std::cos(k * x), c = std::sin(k * x) * std::cos(k * y) - .21;
      const double root = (std::atan2(b, a) - std::acos(-c / std::hypot(a, b))) / k;
      TopoDS_Vertex point; BRep_Builder().MakeVertex(point, gp_Pnt(x, y, root), 1.e-10);
      BRepExtrema_DistShapeShape distance(point, sides.Faces[0]);
      ASSERT_TRUE(distance.IsDone());
      EXPECT_LE(distance.Value(), result.MathematicalToNativeBound + 1.e-8);
      for (int j = 0; j < 12; ++j)
      {
        const double z = s.HeightMinimum + (s.HeightMaximum - s.HeightMinimum) * (j + .41) / 12;
        const auto state = BRepClass_FaceClassifier(sides.Faces[0], gp_Pnt2d(angle, z), 1.e-9).State();
        if (state == TopAbs_IN && ((negative && z > root) || (!negative && z < root)))
          EXPECT_LE(std::abs(z - root), result.NativeToMathematicalBound + 1.e-8);
      }
    }
  }
}

TEST(BRepClass3dImplicitCylinder, ShellGeometryRefusesBudgetsAndSupportsPlacementCancellation)
{
  auto s = settings(); s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
  s.HeightMinimum = -1.7; s.HeightMaximum = 1.3;
  const auto source = field();
  const auto cells = Ops::CheckShellGeometry(source, s, true, 2000000, 24, 200000, 24, 1);
  EXPECT_FALSE(cells.BoundaryVerified); EXPECT_EQ(cells.Cells, 1); EXPECT_EQ(cells.BoundaryBound, 0);
  const auto segments = Ops::CheckShellGeometry(source, s, true, 2000000, 24, 200000, 24, 200000, 1);
  EXPECT_FALSE(segments.BoundaryVerified); EXPECT_EQ(segments.SegmentTests, 1);
  EXPECT_FALSE(Ops::CheckShellGeometry(source, s, true, 2000000, 24, 200000, 0).BoundaryVerified);
  EXPECT_THROW(Ops::CheckShellGeometry(source, s, true, 2000000, 24, 200000, 24, 0), Standard_DomainError);
  gp_Trsf move; move.SetTranslation(gp_Vec(1.2, -.4, 2.1));
  auto placed = source; placed.Location(TopLoc_Location(move)); s.Frame.Transform(move);
  const auto translated = Ops::CheckShellGeometry(placed, s);
  ASSERT_TRUE(translated.BoundaryVerified) << translated.Reason;
  auto reverse = settings(); reverse.AngleMinimum = .27; reverse.AngleMaximum = .27 + 2 * std::acos(-1.);
  reverse.Frame = gp_Ax2(gp_Pnt(), gp_Dir(0, 0, -1), gp_Dir(1, 0, 0));
  reverse.HeightMinimum = -1.3; reverse.HeightMaximum = 1.7;
  const auto reversed = Ops::CheckShellGeometry(source, reverse);
  ASSERT_TRUE(reversed.BoundaryVerified) << reversed.Reason;
  const auto partial = Ops::CheckShellGeometry(source, settings());
  EXPECT_FALSE(partial.BoundaryVerified); EXPECT_EQ(partial.BoundaryBound, 0);
  static size_t calls = 0;
  calls = 0; s.WorkCheckpoint = [] { ++calls; };
  ASSERT_TRUE(Ops::CheckShellGeometry(placed, s).BoundaryVerified);
  ASSERT_GT(calls, 2); --calls;
  s.WorkCheckpoint = [] { if (--calls == 0) throw Standard_AbortiveTransaction("late geometry cancellation"); };
  EXPECT_THROW(Ops::CheckShellGeometry(placed, s), Standard_AbortiveTransaction);
}

TEST(BRepClass3dImplicitCylinder, MaterialClassificationUsesSpatialBoundsAndKeepsAmbiguity)
{
  auto s = settings(); s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
  s.HeightMinimum = -1.7; s.HeightMaximum = 1.3;
  const double root = std::asin(.21) * 10 / (2 * std::acos(-1.));
  const auto source = field();
  std::stringstream before; BRep_ImplicitIO::Write(source, before);
  for (bool negative : {false, true})
  {
    const auto geometry = Ops::CheckShellGeometry(source, s, negative);
    ASSERT_TRUE(geometry.BoundaryVerified) << geometry.Reason;
    const auto classify = [&](const gp_Pnt& p) { return Ops::ClassifyMaterial(geometry, p); };
    EXPECT_EQ(classify(gp_Pnt(0, 0, -.5)).State, negative ? TopAbs_IN : TopAbs_OUT);
    EXPECT_EQ(classify(gp_Pnt(0, 0, .8)).State, negative ? TopAbs_OUT : TopAbs_IN);
    EXPECT_EQ(classify(gp_Pnt(0, 0, root)).State, TopAbs_ON);
    EXPECT_EQ(classify(gp_Pnt(0, 0, -1.7)).State, negative ? TopAbs_ON : TopAbs_OUT);
    EXPECT_EQ(classify(gp_Pnt(0, 0, 1.3)).State, negative ? TopAbs_OUT : TopAbs_ON);
    EXPECT_EQ(classify(gp_Pnt(.5, 0, -.8)).State, negative ? TopAbs_ON : TopAbs_OUT);
    EXPECT_EQ(classify(gp_Pnt(.5, 0, 1)).State, negative ? TopAbs_OUT : TopAbs_ON);
    EXPECT_EQ(classify(gp_Pnt(0, 0, root - 1.e-5)).State, negative ? TopAbs_IN : TopAbs_OUT);
    EXPECT_EQ(classify(gp_Pnt(0, 0, root + 1.e-5)).State, negative ? TopAbs_OUT : TopAbs_IN);
    const auto far = classify(gp_Pnt(2, 0, 0));
    EXPECT_EQ(far.State, TopAbs_OUT); EXPECT_TRUE(far.BoundaryAllowanceSeparated);
    const auto near = classify(gp_Pnt(0, 0, root));
    EXPECT_FALSE(near.BoundaryAllowanceSeparated); EXPECT_LE(near.DistanceUpperBound, 1.e-7);
    const auto margin = classify(gp_Pnt(0, 0, root + .05));
    ASSERT_GT(margin.DistanceUpperBound, margin.DistanceLowerBound);
    const double tolerance = (margin.DistanceLowerBound + margin.DistanceUpperBound) * .5;
    EXPECT_EQ(Ops::ClassifyMaterial(geometry, gp_Pnt(0, 0, root + .05), tolerance).State, TopAbs_UNKNOWN);
    EXPECT_THROW(Ops::ClassifyMaterial(geometry, gp_Pnt(), 0), Standard_DomainError);
    EXPECT_THROW(Ops::ClassifyMaterial(geometry, gp_Pnt(std::numeric_limits<double>::quiet_NaN(), 0, 0)), Standard_DomainError);
  }
  EXPECT_THROW(Ops::ClassifyMaterial(Ops::ShellGeometry(), gp_Pnt()), Standard_DomainError);
  gp_Trsf move; move.SetTranslation(gp_Vec(1.2, -.4, 2.1));
  auto placed = source; placed.Location(TopLoc_Location(move));
  auto reversed = s;
  reversed.Frame = gp_Ax2(gp_Pnt(), gp_Dir(0, 0, -1), gp_Dir(1, 0, 0));
  reversed.Frame.Transform(move); reversed.HeightMinimum = -1.3; reversed.HeightMaximum = 1.7;
  const auto geometry = Ops::CheckShellGeometry(placed, reversed);
  ASSERT_TRUE(geometry.BoundaryVerified) << geometry.Reason;
  EXPECT_EQ(Ops::ClassifyMaterial(geometry, gp_Pnt(0, 0, -.5).Transformed(move)).State, TopAbs_IN);
  EXPECT_EQ(Ops::ClassifyMaterial(geometry, gp_Pnt(0, 0, root).Transformed(move)).State, TopAbs_ON);
  EXPECT_EQ(Ops::ClassifyMaterial(geometry, gp_Pnt(.5, 0, -.8).Transformed(move)).State, TopAbs_ON);
  std::stringstream after; BRep_ImplicitIO::Write(source, after);
  EXPECT_EQ(before.str(), after.str());
}

TEST(BRepClass3dImplicitCylinder, MaterialVolumeBoundsAnalyticReferenceAndBuildsActualSolid)
{
  auto s = settings(); s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
  s.HeightMinimum = -1.7; s.HeightMaximum = 1.3;
  const auto source = field();
  const double k = 2 * std::acos(-1.) / 10;
  const auto reference = [&](int n) {
    double volume = 0;
    for (int i = 0; i < n; ++i)
    {
      const double r0 = s.Radius * i / n, r1 = s.Radius * (i + 1) / n;
      const double radius = (r0 + r1) * .5, weight = std::acos(-1.) * (r1 * r1 - r0 * r0) / n;
      for (int j = 0; j < n; ++j)
      {
        const double theta = 2 * std::acos(-1.) * (j + .5) / n;
        const double x = radius * std::cos(theta), y = radius * std::sin(theta);
        const double a = std::sin(k * y), b = std::cos(k * x), c = std::sin(k * x) * std::cos(k * y) - .21;
        const double root = (std::atan2(b, a) - std::acos(-c / std::hypot(a, b))) / k;
        volume += weight * (root - s.HeightMinimum);
      }
    }
    return volume;
  };
  const double coarse = reference(128), exact = reference(256);
  ASSERT_NEAR(coarse, exact, 1.e-5);
  Ops::VolumeSettings settings;
  double lows = 0, highs = 0;
  for (bool negative : {false, true})
  {
    const auto geometry = Ops::CheckShellGeometry(source, s, negative);
    ASSERT_TRUE(geometry.BoundaryVerified) << geometry.Reason;
    const auto result = Ops::BuildNativeSolid(geometry, settings);
    ASSERT_TRUE(result.SolidValidated) << result.Reason << " cells " << result.Volume.Cells << " roots " << result.Volume.RootTests;
    ASSERT_TRUE(result.Volume.Verified);
    EXPECT_EQ(result.Solid.ShapeType(), TopAbs_SOLID); EXPECT_TRUE(result.Solid.Closed());
    EXPECT_GT(result.Volume.Lower, 0); EXPECT_LT(result.Volume.Upper, result.Volume.CylinderVolume);
    EXPECT_LE(result.Volume.Upper - result.Volume.Lower, settings.AbsoluteTolerance);
    const double expected = negative ? exact : result.Volume.CylinderVolume - exact;
    EXPECT_LE(result.Volume.Lower, expected + 1.e-5); EXPECT_GE(result.Volume.Upper, expected - 1.e-5);
    lows += result.Volume.Lower; highs += result.Volume.Upper;
    int shells = 0, implicit = 0;
    for (TopExp_Explorer it(result.Solid, TopAbs_SHELL); it.More(); it.Next()) ++shells;
    for (TopExp_Explorer it(result.Solid, TopAbs_FACE); it.More(); it.Next())
      implicit += BRep_Tool::IsImplicit(TopoDS::Face(it.Current()));
    EXPECT_EQ(shells, 1); EXPECT_EQ(implicit, 1);
    EXPECT_EQ(Ops::ClassifyMaterial(result.Source, gp_Pnt(0, 0, -.5)).State, negative ? TopAbs_IN : TopAbs_OUT);
  }
  const double full = std::acos(-1.) * s.Radius * s.Radius * (s.HeightMaximum - s.HeightMinimum);
  EXPECT_LE(lows, full); EXPECT_GE(highs, full);
}

TEST(BRepClass3dImplicitCylinder, SolidVolumeRefusesBudgetsPreservesSourceAndHandlesPlacementCancellation)
{
  auto s = settings(); s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
  s.Frame = gp_Ax2(gp_Pnt(), gp_Dir(0, 0, -1), gp_Dir(1, 0, 0));
  s.HeightMinimum = -1.3; s.HeightMaximum = 1.7;
  const auto source = field();
  std::stringstream before; BRep_ImplicitIO::Write(source, before);
  gp_Trsf move; move.SetTranslation(gp_Vec(1.2, -.4, 2.1));
  auto placed = source; placed.Location(TopLoc_Location(move)); s.Frame.Transform(move);
  const auto geometry = Ops::CheckShellGeometry(placed, s);
  ASSERT_TRUE(geometry.BoundaryVerified) << geometry.Reason;
  Ops::VolumeSettings limits;
  const auto solid = Ops::BuildNativeSolid(geometry, limits);
  ASSERT_TRUE(solid.SolidValidated) << solid.Reason;
  EXPECT_EQ(Ops::ClassifyMaterial(solid.Source, gp_Pnt(0, 0, -.5).Transformed(move)).State, TopAbs_IN);
  limits.MaximumCells = 1;
  const auto cells = Ops::BoundMaterialVolume(geometry, limits);
  EXPECT_FALSE(cells.Verified); EXPECT_EQ(cells.Cells, 1); EXPECT_EQ(cells.Lower, 0); EXPECT_EQ(cells.Upper, 0);
  EXPECT_TRUE(Ops::BuildNativeSolid(geometry, limits).Solid.IsNull());
  limits = Ops::VolumeSettings(); limits.MaximumRootTests = 1;
  const auto roots = Ops::BoundMaterialVolume(geometry, limits);
  EXPECT_FALSE(roots.Verified); EXPECT_EQ(roots.RootTests, 1);
  limits = Ops::VolumeSettings(); limits.MaximumDepth = 0;
  EXPECT_FALSE(Ops::BoundMaterialVolume(geometry, limits).Verified);
  limits = Ops::VolumeSettings(); limits.AbsoluteTolerance = 100;
  const auto coarse = Ops::BuildNativeSolid(geometry, limits);
  EXPECT_TRUE(coarse.Volume.Verified); EXPECT_FALSE(coarse.SolidValidated); EXPECT_TRUE(coarse.Solid.IsNull());
  limits.AbsoluteTolerance = 0;
  EXPECT_THROW(Ops::BoundMaterialVolume(geometry, limits), Standard_DomainError);
  limits = Ops::VolumeSettings();
  EXPECT_THROW(Ops::BoundMaterialVolume(Ops::ShellGeometry(), limits), Standard_DomainError);
  static size_t calls = 0;
  calls = 0; limits.WorkCheckpoint = [] { ++calls; };
  ASSERT_TRUE(Ops::BoundMaterialVolume(geometry, limits).Verified);
  ASSERT_GT(calls, 2);
  const size_t full = calls;
  calls = full - 1; limits.WorkCheckpoint = [] { if (--calls == 0) throw Standard_AbortiveTransaction("late volume cancellation"); };
  EXPECT_THROW(Ops::BoundMaterialVolume(geometry, limits), Standard_AbortiveTransaction);
  calls = full + 1;
  EXPECT_THROW(Ops::BuildNativeSolid(geometry, limits), Standard_AbortiveTransaction);
  std::stringstream after; BRep_ImplicitIO::Write(source, after);
  EXPECT_EQ(before.str(), after.str());
}

TEST(BRepClass3dImplicitCylinder, StoredCylinderSolidPreservesSharedGraphAndDropsCertificates)
{
  auto s = settings(); s.AngleMinimum = .27; s.AngleMaximum = .27 + 2 * std::acos(-1.);
  s.HeightMinimum = -1.7; s.HeightMaximum = 1.3;
  for (bool negative : {false, true})
  {
    auto source = field();
    if (negative)
    { gp_Trsf move; move.SetTranslation(gp_Vec(1.2, -.4, 2.1)); source.Location(TopLoc_Location(move)); s.Frame.Transform(move); }
    const auto geometry = Ops::CheckShellGeometry(source, s, negative);
    ASSERT_TRUE(geometry.BoundaryVerified) << geometry.Reason;
    const auto solid = Ops::BuildNativeSolid(geometry, Ops::VolumeSettings());
    ASSERT_TRUE(solid.SolidValidated) << solid.Reason;
    std::ostringstream record; BRepClass3d_ImplicitCylinderIO::Write(solid, record);
    BRepClass3d_ImplicitCylinderIO::StoredSolid restored;
    std::istringstream input(record.str()); BRepClass3d_ImplicitCylinderIO::Read(input, restored);
    ASSERT_TRUE(restored.TopologyValidated); EXPECT_EQ(restored.Shape.ShapeType(), TopAbs_SOLID);
    EXPECT_TRUE(restored.Shape.Closed()); EXPECT_EQ(restored.Faces.size(), 3);
    EXPECT_FALSE(restored.GeometryVerified); EXPECT_FALSE(restored.VolumeVerified); EXPECT_FALSE(restored.SolidValidated);
    EXPECT_EQ(restored.Negative, negative);
    EXPECT_EQ(restored.StoredBoundaryBound, geometry.BoundaryBound);
    EXPECT_EQ(restored.StoredVolumeLower, solid.Volume.Lower); EXPECT_EQ(restored.StoredVolumeUpper, solid.Volume.Upper);
    std::ostringstream a,b; BRep_ImplicitIO::Write(source,a); BRep_ImplicitIO::Write(restored.Primitive,b);
    EXPECT_EQ(a.str(), b.str());
    int implicit = 0;
    for (TopExp_Explorer it(restored.Shape, TopAbs_FACE); it.More(); it.Next())
      implicit += BRep_Tool::IsImplicit(TopoDS::Face(it.Current()));
    EXPECT_EQ(implicit, 1);
    const auto coherence = BRepClass3d_ImplicitCylinderIO::CheckCurveCoherence(restored);
    EXPECT_TRUE(coherence.Verified) << coherence.Reason;
    EXPECT_GT(coherence.Intervals, 0);
    EXPECT_LE(coherence.Bound, s.DistanceTolerance);
    const auto exhausted = BRepClass3d_ImplicitCylinderIO::CheckCurveCoherence(restored, 1);
    EXPECT_FALSE(exhausted.Verified);
    EXPECT_EQ(exhausted.Reason, "cylinder_io_coherence_budget");
    EXPECT_FALSE(restored.GeometryVerified);
    const auto seam = BRepClass3d_ImplicitCylinderIO::CheckPrimitiveSeam(restored);
    EXPECT_TRUE(seam.Verified) << seam.Reason;
    EXPECT_LE(seam.Bound, s.DistanceTolerance);
    const auto maskedSeam = BRepClass3d_ImplicitCylinderIO::CheckMaskedSeam(restored);
    EXPECT_TRUE(maskedSeam.Verified) << maskedSeam.Reason;
    EXPECT_LE(maskedSeam.Bound, s.DistanceTolerance);
    EXPECT_FALSE(BRepClass3d_ImplicitCylinderIO::CheckMaskedSeam(restored,1).Verified);
    const auto coverage = BRepClass3d_ImplicitCylinderIO::CheckSeamCoverage(restored);
    EXPECT_TRUE(coverage.Verified) << coverage.Reason;
    EXPECT_LE(coverage.Bound,s.DistanceTolerance);
    EXPECT_FALSE(BRepClass3d_ImplicitCylinderIO::CheckSeamCoverage(restored,1).Verified);
    const auto axial = Ops::CheckAxialField(restored.Primitive,restored.Parameters);
    EXPECT_TRUE(axial.MonotonicityVerified) << axial.Reason;
    EXPECT_TRUE(axial.Boundary.Wires.empty());
    EXPECT_FALSE(axial.Boundary.IncidenceVerified);
    const auto retainedCap = BRepClass3d_ImplicitCylinderIO::CheckRetainedCap(restored);
    EXPECT_TRUE(retainedCap.Verified) << retainedCap.Reason;
    EXPECT_LE(retainedCap.Bound,s.DistanceTolerance);
    EXPECT_FALSE(BRepClass3d_ImplicitCylinderIO::CheckRetainedCap(restored,1).Verified);
    auto provenance = restored;
    provenance.StoredAxialDirection *= -1;
    EXPECT_TRUE(BRepClass3d_ImplicitCylinderIO::CheckRetainedCap(provenance).Verified);
    const auto capPlane = occ::down_cast<Geom_Plane>(BRep_Tool::Surface(restored.Faces[0]));
    const auto capLocation = capPlane->Location();
    capPlane->SetLocation(capLocation.Translated(gp_Vec(s.Frame.Direction())*.2));
    EXPECT_FALSE(BRepClass3d_ImplicitCylinderIO::CheckRetainedCap(restored).Verified);
    capPlane->SetLocation(capLocation);
    EXPECT_FALSE(BRepClass3d_ImplicitCylinderIO::CheckPrimitiveSeam(restored,1).Verified);
    const auto previous = restored.Shape;
    std::vector<std::string> corruptions {record.str().substr(0, record.str().size()/2), record.str()+"trailing"};
    const auto replace = [&](const std::string& from, const std::string& to) {
      auto bad = record.str(); const auto at = bad.find(from);
      EXPECT_NE(at, std::string::npos); if (at != std::string::npos) bad.replace(at, from.size(), to);
      corruptions.push_back(std::move(bad));
    };
    replace("OCCT_IMPLICIT_CYLINDER_SOLID 1 ", "OCCT_IMPLICIT_CYLINDER_SOLID 9 ");
    replace("VERTICES 3\n", "VERTICES 4\n");
    replace("EDGES 4\n", "EDGES 1000000\n");
    auto poles = record.str(); const auto at = poles.find(" B ");
    ASSERT_NE(at, std::string::npos); const auto end = poles.find(' ', at + 3);
    ASSERT_NE(end, std::string::npos); poles.replace(at + 3, end - at - 3, "65537");
    corruptions.push_back(std::move(poles));
    for (const auto& bad : corruptions)
    {
      std::istringstream corrupt(bad);
      EXPECT_THROW(BRepClass3d_ImplicitCylinderIO::Read(corrupt, restored), Standard_Failure);
      EXPECT_TRUE(restored.Shape.IsSame(previous));
    }
    if (negative)
    {
      std::istringstream oversized(std::string(16 * 1024 * 1024 + 1, 'x'));
      EXPECT_THROW(BRepClass3d_ImplicitCylinderIO::Read(oversized, restored), Standard_Failure);
      EXPECT_TRUE(restored.Shape.IsSame(previous));
    }
    auto invalid = solid; invalid.SolidValidated = false;
    std::ostringstream untouched; untouched << "prefix";
    EXPECT_THROW(BRepClass3d_ImplicitCylinderIO::Write(invalid, untouched), Standard_Failure);
    EXPECT_EQ(untouched.str(), "prefix");
    bool changed = false;
    for (TopExp_Explorer it(restored.Shape, TopAbs_EDGE); it.More() && !changed; it.Next())
    {
      double first, last;
      const auto spline = occ::down_cast<Geom_BSplineCurve>(
        BRep_Tool::Curve(TopoDS::Edge(it.Current()), first, last));
      if (spline.IsNull() || spline->NbPoles() < 3) continue;
      const int pole = 1 + spline->NbPoles()/2;
      spline->SetPole(pole, spline->Pole(pole).Translated(gp_Vec(0,0,1)));
      changed = true;
    }
    ASSERT_TRUE(changed);
    EXPECT_FALSE(BRepClass3d_ImplicitCylinderIO::CheckCurveCoherence(restored).Verified);
    EXPECT_FALSE(BRepClass3d_ImplicitCylinderIO::CheckPrimitiveSeam(restored).Verified);
    EXPECT_FALSE(BRepClass3d_ImplicitCylinderIO::CheckMaskedSeam(restored).Verified);
    EXPECT_FALSE(BRepClass3d_ImplicitCylinderIO::CheckSeamCoverage(restored).Verified);
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
