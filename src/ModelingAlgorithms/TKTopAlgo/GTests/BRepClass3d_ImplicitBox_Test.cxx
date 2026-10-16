// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitBox.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_ImplicitIO.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <Standard_NotImplemented.hxx>
#include <Standard_DomainError.hxx>
#include <gtest/gtest.h>
#include <sstream>
#include <cmath>
#include <array>
#include <limits>

namespace
{
TopoDS_Face field(double iso = 1.45)
{
  TopoDS_Face face;
  BRep_Builder().MakeImplicitFace(face,
                                  new Geom_ImplicitGyroid(10, iso),
                                  gp_Pnt(-5, -5, -5),
                                  gp_Pnt(16, 16, 6),
                                  TopLoc_Location(),
                                  .1);
  return face;
}

TopoDS_Shape box()
{
  return BRepPrimAPI_MakeBox(gp_Pnt(-.31, -.37, 1.25), gp_Pnt(2.93, 2.91, 4)).Shape();
}
} // namespace

TEST(BRepClass3d_ImplicitBox, CommonAndDirectionalCutShareOppositeSeamsAndCapPCurves)
{
  const auto source  = field();
  const auto operand = box();
  const bool wasFree = operand.Free();
  for (auto operation : {BRepClass3d_ImplicitBox::Operation::Common,
                         BRepClass3d_ImplicitBox::Operation::BoxMinusField})
  {
    const auto result = BRepClass3d_ImplicitBox::Build(source, operand, operation);
    EXPECT_EQ(result.Shape.ShapeType(), TopAbs_COMPOUND);
    ASSERT_EQ(result.SeamCount, 1);
    ASSERT_EQ(result.Sections.size(), 6u);
    EXPECT_LE(result.GeometricBound, .1);
    EXPECT_EQ(result.Faces.size(),
              operation == BRepClass3d_ImplicitBox::Operation::Common ? 7u : 2u);
    const auto implicit = result.Faces.back();
    EXPECT_TRUE(BRep_Tool::IsImplicit(implicit));
    EXPECT_FALSE(BRep_Tool::HasParametricSurface(implicit));
    EXPECT_EQ(implicit.Orientation(),
              operation == BRepClass3d_ImplicitBox::Operation::Common ? TopAbs_FORWARD
                                                                      : TopAbs_REVERSED);
    EXPECT_FALSE(BRep_ImplicitTool::IsInDomain(implicit, gp_Pnt(0, 0, 0)));
    int count = 0;
    for (TopExp_Explorer seam(implicit, TopAbs_EDGE); seam.More(); seam.Next())
    {
      int uses = 0;
      for (size_t i = 0; i + 1 < result.Faces.size(); ++i)
      {
        EXPECT_TRUE(BRepCheck_Analyzer(result.Faces[i]).IsValid());
        for (TopExp_Explorer edge(result.Faces[i], TopAbs_EDGE); edge.More(); edge.Next())
          if (edge.Current().IsSame(seam.Current()))
          {
            EXPECT_NE(edge.Current().Orientation(), seam.Current().Orientation());
            double     first, last;
            const auto pc =
              BRep_Tool::CurveOnSurface(TopoDS::Edge(edge.Current()), result.Faces[i], first, last);
            EXPECT_FALSE(pc.IsNull());
            EXPECT_GT(last, first);
            ++uses;
          }
      }
      EXPECT_EQ(uses, 1);
      ++count;
    }
    EXPECT_EQ(count, 1);
    std::ostringstream record;
    EXPECT_THROW(BRep_ImplicitIO::Write(implicit, record), Standard_NotImplemented);
  }
  EXPECT_FALSE(TopoDS_Iterator(source).More());
  EXPECT_EQ(operand.Free(), wasFree);
  EXPECT_TRUE(BRepCheck_Analyzer(operand).IsValid());
}

TEST(BRepClass3d_ImplicitBox, UnsupportedInputsRejectRatherThanReturnPartialBoundary)
{
  using Ops = BRepClass3d_ImplicitBox;
  EXPECT_THROW(Ops::Build(field(0), BRepPrimAPI_MakeBox(3, 3, 3).Shape(), Ops::Operation::Common),
               Standard_NotImplemented); // A corner is exactly on F=0.
  EXPECT_THROW(Ops::Build(field(), BRepPrimAPI_MakeCylinder(1, 2).Shape(), Ops::Operation::Common),
               Standard_NotImplemented);
  EXPECT_THROW(Ops::Build(field(), box(), Ops::Operation::Common, .1, 4), Standard_NotImplemented);
  EXPECT_THROW(Ops::Build(field(), box().Reversed(), Ops::Operation::Common),
               Standard_NotImplemented);
  EXPECT_THROW(Ops::Build(field(), box(), Ops::Operation::Common, 0), Standard_DomainError);
  gp_Trsf translation;
  translation.SetTranslation(gp_Vec(20, 0, 0));
  EXPECT_THROW(
    Ops::Build(field(), box().Moved(TopLoc_Location(translation)), Ops::Operation::Common),
    Standard_NotImplemented);
}

TEST(BRepClass3d_ImplicitBox, OpenContoursJoinAcrossFacesWithAnalyticEdgeReference)
{
  using Ops = BRepClass3d_ImplicitBox;
  const gp_Pnt low(-3.83, -4.07, -2.67), high(4.17, 3.93, 3.37);
  const auto   operand = BRepPrimAPI_MakeBox(low, high).Shape();
  // Independent harmonic solutions along each box edge: A sin(k*t) +
  // B cos(k*t) + C = 0. This does not use the kernel root solver/gradients.
  const double        pi = std::acos(-1.), k = 2 * pi / 10;
  std::vector<gp_Pnt> reference;
  for (int axis = 0; axis < 3; ++axis)
    for (int bits = 0; bits < 8; ++bits)
    {
      if (bits & (1 << axis))
        continue;
      gp_Pnt       p(bits & 1 ? high.X() : low.X(),
                     bits & 2 ? high.Y() : low.Y(),
                     bits & 4 ? high.Z() : low.Z());
      const int    next = (axis + 1) % 3, prev = (axis + 2) % 3;
      const double a = std::cos(k * p.Coord(next + 1));
      const double b = std::sin(k * p.Coord(prev + 1));
      const double c = std::sin(k * p.Coord(next + 1)) * std::cos(k * p.Coord(prev + 1));
      const double r = std::hypot(a, b);
      if (std::abs(c) >= r)
        continue;
      const double phase = std::atan2(b, a), angle = std::asin(-c / r);
      for (double base : {angle, pi - angle})
        for (int n = -3; n <= 3; ++n)
        {
          const double coordinate = (base - phase + 2 * pi * n) / k;
          if (coordinate <= low.Coord(axis + 1) || coordinate >= high.Coord(axis + 1))
            continue;
          p.SetCoord(axis + 1, coordinate);
          reference.push_back(p);
        }
    }
  ASSERT_GT(reference.size(), 4u);
  for (auto operation : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
  {
    const auto result = Ops::Build(field(0), operand, operation);
    EXPECT_EQ(result.BoxEdgeRootCount, int(reference.size()));
    EXPECT_LE(result.GeometricBound, .1);
    std::vector<int> matches(reference.size(), 0);
    int              open = 0;
    for (const auto& section : result.Sections)
    {
      EXPECT_TRUE(section.Resolved);
      for (const auto& curve : section.Curves)
      {
        if (curve.Closed)
          continue;
        ++open;
        for (const auto& point : {curve.Points.front(), curve.Points.back()})
        {
          size_t nearest = 0;
          for (size_t i = 1; i < reference.size(); ++i)
            if (point.Distance(reference[i]) < point.Distance(reference[nearest]))
              nearest = i;
          EXPECT_LE(point.Distance(reference[nearest]), .1 / 256);
          ++matches[nearest];
        }
      }
    }
    EXPECT_GT(open, 0);
    for (int count : matches)
      EXPECT_EQ(count, 2);
    for (size_t i = 0; i + 1 < result.Faces.size(); ++i)
      EXPECT_TRUE(BRepCheck_Analyzer(result.Faces[i]).IsValid());
    const auto implicit = result.Faces.back();
    for (TopoDS_Iterator it(implicit); it.More(); it.Next())
      EXPECT_TRUE(it.Value().Closed());
  }
}

TEST(BRepClass3d_ImplicitBox, TranslatedBoxAndNoSeamCases)
{
  using Ops = BRepClass3d_ImplicitBox;
  gp_Trsf translation;
  translation.SetTranslation(gp_Vec(.01, .02, .03));
  const auto moved    = box().Moved(TopLoc_Location(translation));
  const auto boundary = Ops::Build(field(), moved, Ops::Operation::Common);
  EXPECT_EQ(boundary.SeamCount, 1);
  const auto emptyField = Ops::Build(field(-3), box(), Ops::Operation::Common);
  EXPECT_EQ(emptyField.SeamCount, 0);
  EXPECT_EQ(emptyField.Faces.size(), 1u); // A bounded implicit face can be empty.
  const auto fullBox = Ops::Build(field(3), box(), Ops::Operation::Common);
  EXPECT_EQ(fullBox.SeamCount, 0);
  EXPECT_EQ(fullBox.Faces.size(), 7u); // Six caps plus possibly empty native face.
}

TEST(BRepClass3d_ImplicitBox, MultipleClosedLoopsCreateHolesOrSeparateRetainedCaps)
{
  using Ops = BRepClass3d_ImplicitBox;
  const auto operand =
    BRepPrimAPI_MakeBox(gp_Pnt(-.31, -.37, 1.25), gp_Pnt(12.93, 2.91, 4)).Shape();
  const auto common = Ops::Build(field(), operand, Ops::Operation::Common);
  const auto cut    = Ops::Build(field(), operand, Ops::Operation::BoxMinusField);
  EXPECT_EQ(common.BoxEdgeRootCount, 0);
  EXPECT_EQ(common.SeamCount, 2);
  EXPECT_EQ(cut.SeamCount, 2);
  EXPECT_EQ(common.Faces.size(), 7u);
  EXPECT_EQ(cut.Faces.size(), 3u);
  bool holes = false;
  for (size_t i = 0; i + 1 < common.Faces.size(); ++i)
  {
    int wires = 0;
    for (TopoDS_Iterator it(common.Faces[i]); it.More(); it.Next())
      ++wires;
    holes |= wires == 3;
  }
  EXPECT_TRUE(holes);
}

TEST(BRepClass3d_ImplicitBox, LevelChangesAndPlacedSourceRemainExplicit)
{
  using Ops = BRepClass3d_ImplicitBox;
  const auto operand =
    BRepPrimAPI_MakeBox(gp_Pnt(-3.83, -4.07, -2.67), gp_Pnt(4.17, 3.93, 3.37)).Shape();
  for (double iso : {-.21, .21})
    for (auto operation : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
    {
      const auto result = Ops::Build(field(iso), operand, operation);
      EXPECT_GT(result.BoxEdgeRootCount, 0);
      EXPECT_LE(result.GeometricBound, .1);
    }
  gp_Trsf moved;
  moved.SetTranslation(gp_Vec(.1, 0, 0));
  EXPECT_THROW(Ops::Build(TopoDS::Face(field().Moved(TopLoc_Location(moved))),
                          operand,
                          Ops::Operation::Common),
               Standard_NotImplemented);
  // A plane through an isolated maximum does not become an empty-success
  // section merely because the new Taylor/harmonic exclusions are available.
  EXPECT_THROW(Ops::Build(field(1.5), box(), Ops::Operation::Common), Standard_NotImplemented);
}

TEST(BRepClass3d_ImplicitBox, ClassifierUsesRegionAndPhysicalDistanceNotFieldResidual)
{
  using Ops         = BRepClass3d_ImplicitBox;
  const auto common = Ops::Build(field(), box(), Ops::Operation::Common);
  const auto cut    = Ops::Build(field(), box(), Ops::Operation::BoxMinusField);
  for (const auto& point : {gp_Pnt(0, 0, 2), gp_Pnt(1.25, 1.25, 1.5)})
  {
    const bool negative = point.X() == 0;
    const auto a = Ops::Classify(common, point), b = Ops::Classify(cut, point);
    EXPECT_EQ(a.State, negative ? Ops::PointState::Inside : Ops::PointState::Outside);
    EXPECT_EQ(b.State, negative ? Ops::PointState::Outside : Ops::PointState::Inside);
    EXPECT_GT(a.DistanceLowerBound, .1);
    EXPECT_GT(b.DistanceLowerBound, .1);
  }
  // F at the stationary maximum is only .05, dimensionless, yet the nearest
  // field boundary is more than .1 mm away. Taylor exclusion proves Outside.
  const auto peak = Ops::Classify(common, gp_Pnt(1.25, 1.25, 1.25));
  EXPECT_EQ(peak.State, Ops::PointState::Outside);
  EXPECT_GT(peak.DistanceLowerBound, .1);
  for (double x : {-.30, -.33})
  {
    const gp_Pnt point(x, 0, 2);
    const auto   a = Ops::Classify(common, point);
    EXPECT_EQ(a.State, Ops::PointState::NearBoundary);
    EXPECT_TRUE(a.HasBoundaryWitness);
    EXPECT_NEAR(a.WitnessFirst.X(), -.31, 1.e-14);
    EXPECT_LE(a.WitnessDistanceUpperBound, .1);
    EXPECT_EQ(Ops::Classify(cut, point).State, Ops::PointState::Outside);
  }
  EXPECT_EQ(Ops::Classify(common, gp_Pnt(30, 20, 10)).State, Ops::PointState::Outside);
  EXPECT_THROW(Ops::Classify(common, gp_Pnt(0, 0, 0), 0), Standard_DomainError);
  EXPECT_THROW(Ops::Classify(common, gp_Pnt(std::numeric_limits<double>::infinity(), 0, 0)),
               Standard_DomainError);
  EXPECT_THROW(Ops::Classify(Ops::Boundary(), gp_Pnt(0, 0, 0)), Standard_NotImplemented);
}

TEST(BRepClass3d_ImplicitBox, NearFieldBoundaryHasClippedSignBracketAndCriticalCaseStaysUnresolved)
{
  using Ops = BRepClass3d_ImplicitBox;
  const auto operand =
    BRepPrimAPI_MakeBox(gp_Pnt(-3.83, -4.07, -2.67), gp_Pnt(4.17, 3.93, 3.37)).Shape();
  const auto zero  = Ops::Build(field(0), operand, Ops::Operation::Common);
  const auto query = Ops::Classify(zero, gp_Pnt(0, 0, 0));
  EXPECT_EQ(query.State, Ops::PointState::NearBoundary);
  EXPECT_TRUE(query.HasBoundaryWitness);
  EXPECT_LE(query.WitnessBracketWidth, .1 / 256);
  EXPECT_LE(query.WitnessDistanceUpperBound, .1);
  const auto& face = zero.Faces.back();
  EXPECT_LT(BRep_ImplicitTool::Value(face, query.WitnessFirst)
              * BRep_ImplicitTool::Value(face, query.WitnessLast),
            0);
  for (const auto& point : {query.WitnessFirst, query.WitnessLast})
    EXPECT_TRUE(BRep_ImplicitTool::IsInDomain(face, point));
  // A coarse tolerance ball intersects a boundary but the stationary gradient
  // offers no root bracket and no cap witness: fail explicitly, never guess ON.
  const auto common   = Ops::Build(field(), box(), Ops::Operation::Common);
  const auto critical = Ops::Classify(common, gp_Pnt(1.25, 1.25, 1.25), .4);
  EXPECT_EQ(critical.State, Ops::PointState::Unresolved);
  EXPECT_FALSE(critical.HasBoundaryWitness);
}

TEST(BRepClass3d_ImplicitBox, ClassifierLengthToleranceScalesIndependentlyOfFieldUnits)
{
  using Ops = BRepClass3d_ImplicitBox;
  TopoDS_Face small;
  BRep_Builder().MakeImplicitFace(small,
                                  new Geom_ImplicitGyroid(1, 1.45),
                                  gp_Pnt(-.5, -.5, -.5),
                                  gp_Pnt(.6, .6, .6),
                                  TopLoc_Location(),
                                  .01);
  const auto smallBox =
    BRepPrimAPI_MakeBox(gp_Pnt(-.031, -.037, .125), gp_Pnt(.293, .291, .4)).Shape();
  const auto   a = Ops::Build(field(), box(), Ops::Operation::BoxMinusField);
  const auto   b = Ops::Build(small, smallBox, Ops::Operation::BoxMinusField, .01);
  const gp_Pnt p(1.25, 1.25, 1.5), q(.125, .125, .15);
  EXPECT_NEAR(BRep_ImplicitTool::Value(a.Faces.back(), p),
              BRep_ImplicitTool::Value(b.Faces.back(), q),
              1.e-14);
  EXPECT_EQ(Ops::Classify(a, p, .1).State, Ops::PointState::Inside);
  EXPECT_EQ(Ops::Classify(b, q, .1).State, Ops::PointState::NearBoundary);
  EXPECT_EQ(Ops::Classify(b, q, .01).State, Ops::PointState::Inside);
}
