// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitBox.hxx>
#include <BRepClass3d_ImplicitBoxIO.hxx>
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
#include <Poly_Triangulation.hxx>

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

TEST(BRepClass3d_ImplicitBox,
     TopologyRoundTripPreservesSharedSeamsAndClassificationWithoutRetracing)
{
  using Ops = BRepClass3d_ImplicitBox;
  const auto openBox =
    BRepPrimAPI_MakeBox(gp_Pnt(-3.83, -4.07, -2.67), gp_Pnt(4.17, 3.93, 3.37)).Shape();
  for (bool open : {false, true})
    for (auto operation : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
    {
      const auto source = Ops::Build(field(open ? 0 : 1.45), open ? openBox : box(), operation);
      ASSERT_TRUE(source.SeamBoundVerified);
      std::ostringstream record;
      BRepClass3d_ImplicitBoxIO::Write(source, record);
      EXPECT_EQ(record.str().find("OCCT_IMPLICIT_BOX_TOPOLOGY 1 "), 0u);
      Ops::Boundary      restored;
      std::istringstream input(record.str());
      BRepClass3d_ImplicitBoxIO::Read(input, restored);
      EXPECT_FALSE(restored.SeamBoundVerified);
      EXPECT_EQ(restored.Faces.size(), source.Faces.size());
      EXPECT_EQ(restored.Sources, source.Sources);
      EXPECT_EQ(restored.SeamCount, source.SeamCount);
      EXPECT_EQ(restored.BoxEdgeRootCount, source.BoxEdgeRootCount);
      EXPECT_DOUBLE_EQ(restored.GeometricBound, source.GeometricBound);
      EXPECT_EQ(restored.Faces.back().Orientation(), source.Faces.back().Orientation());
      EXPECT_TRUE(BRep_Tool::IsImplicit(restored.Faces.back()));
      EXPECT_FALSE(BRep_Tool::HasParametricSurface(restored.Faces.back()));
      EXPECT_FALSE(restored.Sections.front().Resolved); // No tracing certificate invented by I/O.
      for (const auto& p : {gp_Pnt(0, 0, 2), gp_Pnt(1.25, 1.25, 1.5), gp_Pnt(-.33, 0, 2)})
        EXPECT_EQ(Ops::Classify(restored, p).State, Ops::Classify(source, p).State);
      for (TopExp_Explorer edge(restored.Faces.back(), TopAbs_EDGE); edge.More(); edge.Next())
      {
        int shared = 0;
        for (size_t i = 0; i + 1 < restored.Faces.size(); ++i)
          for (TopExp_Explorer cap(restored.Faces[i], TopAbs_EDGE); cap.More(); cap.Next())
            if (cap.Current().IsSame(edge.Current()))
            {
              EXPECT_NE(cap.Current().Orientation(), edge.Current().Orientation());
              double first, last;
              EXPECT_FALSE(BRep_Tool::CurveOnSurface(TopoDS::Edge(cap.Current()),
                                                     restored.Faces[i],
                                                     first,
                                                     last)
                             .IsNull());
              ++shared;
            }
        EXPECT_EQ(shared, 1);
      }
      std::ostringstream again;
      BRepClass3d_ImplicitBoxIO::Write(restored, again);
      Ops::Boundary      twice;
      std::istringstream second(again.str());
      BRepClass3d_ImplicitBoxIO::Read(second, twice);
      EXPECT_EQ(twice.BoxEdgeRootCount, source.BoxEdgeRootCount);
    }
}

TEST(BRepClass3d_ImplicitBox, TopologyIORejectsMalformedRecordsAtomicallyAndOmitsMeshes)
{
  using Ops                        = BRepClass3d_ImplicitBox;
  auto               source        = Ops::Build(field(), box(), Ops::Operation::Common);
  const auto         originalShape = source.Shape;
  std::ostringstream output;
  auto               broken = source;
  ++broken.SeamCount;
  EXPECT_THROW(BRepClass3d_ImplicitBoxIO::Write(broken, output), Standard_DomainError);
  EXPECT_TRUE(output.str().empty());
  const occ::handle<Poly_Triangulation> mesh = new Poly_Triangulation(3, 1, false);
  mesh->SetNode(1, gp_Pnt(-.31, -.37, 1.25));
  mesh->SetNode(2, gp_Pnt(-.31, 2.91, 1.25));
  mesh->SetNode(3, gp_Pnt(-.31, -.37, 4));
  mesh->SetTriangle(1, Poly_Triangle(1, 2, 3));
  BRep_Builder().UpdateFace(source.Faces.front(), mesh);
  BRepClass3d_ImplicitBoxIO::Write(source, output);
  const auto record = output.str();
  for (const auto& bad : {std::string("garbage"),
                          record + "trailing",
                          record.substr(0, record.size() / 2),
                          std::string("OCCT_IMPLICIT_BOX_TOPOLOGY 99 1 1 1 0 0 0\n-1 \nxx")})
  {
    std::istringstream input(bad);
    EXPECT_THROW(BRepClass3d_ImplicitBoxIO::Read(input, source), Standard_DomainError);
    EXPECT_TRUE(source.Shape.IsSame(originalShape));
  }
  Ops::Boundary      restored;
  std::istringstream input(record);
  BRepClass3d_ImplicitBoxIO::Read(input, restored);
  TopLoc_Location location;
  EXPECT_FALSE(BRep_Tool::Triangulation(source.Faces.front(), location).IsNull());
  EXPECT_TRUE(BRep_Tool::Triangulation(restored.Faces.front(), location).IsNull());
  source.Shape.Locked(true);
  std::istringstream locked(record);
  EXPECT_THROW(BRepClass3d_ImplicitBoxIO::Read(locked, source), Standard_DomainError);
  source.Shape.Locked(false);
}

TEST(BRepClass3d_ImplicitBox, WholeBoxStratifiedRegularityBothOperationsAndStoredTopology)
{
  using Ops = BRepClass3d_ImplicitBox;
  const auto open =
    BRepPrimAPI_MakeBox(gp_Pnt(-3.83, -4.07, -2.67), gp_Pnt(4.17, 3.93, 3.37)).Shape();
  for (bool crossing : {false, true})
    for (auto op : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
    {
      const auto source = Ops::Build(field(crossing ? 0 : 1.45), crossing ? open : box(), op);
      std::ostringstream record;
      BRepClass3d_ImplicitBoxIO::Write(source, record);
      const auto report = Ops::CheckRegularity(source.Faces.back());
      std::cout << "regularity " << crossing << " " << int(op) << " resolved=" << report.Resolved
                << " dimension=" << report.UnresolvedDimension << " visited=" << report.Visited[3]
                << " bounds=" << report.MinimumDerivative[3] << "," << report.MinimumDerivative[2]
                << "," << report.MinimumDerivative[1] << std::endl;
      EXPECT_TRUE(report.Resolved);
      EXPECT_EQ(report.UnresolvedDimension, -1);
      for (bool complete : report.Complete)
        EXPECT_TRUE(complete);
      EXPECT_EQ(report.Visited[0], 8u);
      EXPECT_EQ(report.Excluded[0], 8u);
      EXPECT_GT(report.MinimumDerivative[3], 0);
      EXPECT_GT(report.MinimumDerivative[2], 0);
      if (crossing)
        EXPECT_GT(report.MinimumDerivative[1], 0);
      EXPECT_FALSE(Ops::CheckRegularity(source.Faces.back(), 14, 1).Resolved);
      Ops::Boundary      loaded;
      std::istringstream input(record.str());
      BRepClass3d_ImplicitBoxIO::Read(input, loaded);
      EXPECT_FALSE(loaded.SeamBoundVerified);
      EXPECT_TRUE(Ops::CheckRegularity(loaded.Faces.back()).Resolved);
      std::ostringstream after;
      BRepClass3d_ImplicitBoxIO::Write(source, after);
      EXPECT_EQ(record.str(), after.str());
    }
}

TEST(BRepClass3d_ImplicitBox, RegularityRefusesVolumeFaceEdgeAndCornerCriticalCases)
{
  using Ops       = BRepClass3d_ImplicitBox;
  const auto make = [](double iso, const gp_Pnt& low, const gp_Pnt& high) {
    TopoDS_Face face;
    BRep_Builder()
      .MakeImplicitFace(face, new Geom_ImplicitGyroid(10, iso), low, high, TopLoc_Location(), .001);
    return face;
  };
  const auto critical = make(1.5, gp_Pnt(.7, .8, .9), gp_Pnt(1.8, 1.9, 2.0));
  const auto volume   = Ops::CheckRegularity(critical);
  EXPECT_FALSE(volume.Resolved);
  EXPECT_EQ(volume.UnresolvedDimension, 3);
  const auto tangent = make(1, gp_Pnt(-.5, 2, 1.8), gp_Pnt(.5, 3, 2.5));
  const auto face    = Ops::CheckRegularity(tangent);
  EXPECT_FALSE(face.Resolved);
  EXPECT_TRUE(face.Complete[3]);
  EXPECT_EQ(face.UnresolvedDimension, 2);
  const double k = 2 * std::acos(-1.) / 10, x = .4 / k, y = .6 / k;
  const double z         = std::atan2(std::cos(.4), std::sin(.6)) / k;
  const double iso       = std::sin(.4) * std::cos(.6) + std::hypot(std::sin(.6), std::cos(.4));
  const auto   edgeField = make(iso, gp_Pnt(x, y, z - .05), gp_Pnt(x + .05, y + .05, z + .05));
  const auto   edge      = Ops::CheckRegularity(edgeField);
  EXPECT_FALSE(edge.Resolved);
  EXPECT_TRUE(edge.Complete[3]);
  EXPECT_TRUE(edge.Complete[2]);
  EXPECT_EQ(edge.UnresolvedDimension, 1);
  const auto corner = Ops::CheckRegularity(make(0, gp_Pnt(0, 0, 0), gp_Pnt(.2, .3, .4)));
  EXPECT_FALSE(corner.Resolved);
  EXPECT_TRUE(corner.Complete[3]);
  EXPECT_TRUE(corner.Complete[2]);
  EXPECT_TRUE(corner.Complete[1]);
  EXPECT_EQ(corner.UnresolvedDimension, 0);
  EXPECT_THROW(Ops::CheckRegularity(critical, 0), Standard_DomainError);
  EXPECT_THROW(Ops::CheckRegularity(critical, 14, 0), Standard_DomainError);
}

TEST(BRepClass3d_ImplicitBox, RegularityDerivativeBoundsScaleWithModelLength)
{
  using Ops       = BRepClass3d_ImplicitBox;
  const auto make = [](double scale, bool translated) {
    const gp_Vec shift = translated ? gp_Vec(20 * scale, -10 * scale, 30 * scale) : gp_Vec(0, 0, 0);
    TopoDS_Face  face;
    BRep_Builder().MakeImplicitFace(
      face,
      new Geom_ImplicitGyroid(10 * scale, 1.45),
      gp_Pnt(-.31 * scale, -.37 * scale, 1.25 * scale).Translated(shift),
      gp_Pnt(2.93 * scale, 2.91 * scale, 4 * scale).Translated(shift),
      TopLoc_Location(),
      .001);
    return face;
  };
  const auto original = Ops::CheckRegularity(make(1, false));
  ASSERT_TRUE(original.Resolved);
  for (double scale : {.5, 2.})
  {
    const auto report = Ops::CheckRegularity(make(scale, false));
    EXPECT_TRUE(report.Resolved);
    for (int dimension = 1; dimension <= 3; ++dimension)
      EXPECT_NEAR(report.MinimumDerivative[size_t(dimension)] * scale,
                  original.MinimumDerivative[size_t(dimension)],
                  1.e-10);
    EXPECT_NEAR(report.HessianBound * scale * scale, original.HessianBound, 1.e-12);
    EXPECT_TRUE(Ops::CheckRegularity(make(scale, true)).Resolved);
  }
  const auto conventional = TopoDS::Face(TopExp_Explorer(box(), TopAbs_FACE).Current());
  EXPECT_THROW(Ops::CheckRegularity(conventional), Standard_NotImplemented);
}
