// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitBox.hxx>
#include <BRepClass3d_ImplicitBoxIO.hxx>
#include <BRepClass3d_ImplicitSolidIO.hxx>
#include <BRepMesh_ImplicitBoundary.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRep_ImplicitTool.hxx>
#include <BRep_ImplicitIO.hxx>
#include <BRep_TImplicitFace.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <ElSLib.hxx>
#include <TopExp.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <set>
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
#include <algorithm>
#include <map>
#include <limits>
#include <Poly_Triangulation.hxx>
#include <locale>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>

namespace
{
TopoDS_Face topologyField(const gp_Pnt& low, const gp_Pnt& high, double iso, double scale = 1)
{
  TopoDS_Face face;
  BRep_Builder().MakeImplicitFace(face,
                                  new Geom_ImplicitGyroid(10 * scale, iso),
                                  low,
                                  high,
                                  TopLoc_Location(),
                                  .1 * scale);
  return face;
}

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

TEST(BRepClass3d_ImplicitBox, WholeBoxHomotopyReferenceAndBoundaryTopology)
{
  for (double iso : {1.45, 0.0})
  {
    const auto face   = topologyField(gp_Pnt(-.31, -.37, 1.25), gp_Pnt(2.93, 2.91, 4), iso);
    const auto result = BRepClass3d_ImplicitBox::CheckTopology(face, 32, 16);
    ASSERT_TRUE(result.HomotopyResolved)
      << iso << " " << result.Reason << " dim " << result.UnresolvedDimension << " at "
      << result.UnresolvedMinimum.X() << "," << result.UnresolvedMinimum.Y() << ","
      << result.UnresolvedMinimum.Z();
    ASSERT_TRUE(result.ReferenceManifold) << result.Reason;
    EXPECT_FALSE(result.Nodes.empty());
    EXPECT_FALSE(result.Triangles.empty());
    EXPECT_FALSE(result.BoundaryLoops.empty());
    EXPECT_EQ(result.ComponentEuler, std::vector<int>({1}));
    EXPECT_EQ(result.BoundaryLoops.size(), 1u);
    EXPECT_GT(result.MinimumHomotopyDerivative, 0);
    for (size_t dimension = 0; dimension < 4; ++dimension)
      EXPECT_GT(result.Visited[dimension], 0u);
    auto reversed = face;
    reversed.Reverse();
    const auto other = BRepClass3d_ImplicitBox::CheckTopology(reversed, 32, 16);
    ASSERT_TRUE(other.HomotopyResolved);
    ASSERT_TRUE(other.ReferenceManifold);
    EXPECT_EQ(other.ComponentEuler, result.ComponentEuler);
    EXPECT_EQ(other.Nodes.size(), result.Nodes.size());
    EXPECT_EQ(other.Triangles.size(), result.Triangles.size());
    for (size_t i = 0; i < result.Triangles.size(); ++i)
    {
      EXPECT_EQ(other.Triangles[i][0], result.Triangles[i][0]);
      EXPECT_EQ(other.Triangles[i][1], result.Triangles[i][2]);
      EXPECT_EQ(other.Triangles[i][2], result.Triangles[i][1]);
    }
  }
}

TEST(BRepClass3d_ImplicitBox, WholeBoxHomotopyDetectsHiddenComponentsAndLimits)
{
  const auto hidden = topologyField(gp_Pnt(.3, .3, .3), gp_Pnt(2.3, 2.3, 2.3), 1.4999);
  // All resolution-four nodal samples are negative, but the exact central
  // maximum is positive: accepting signs/incidence alone would lose a sphere.
  for (int z = 0; z <= 4; ++z)
    for (int y = 0; y <= 4; ++y)
      for (int x = 0; x <= 4; ++x)
        EXPECT_LT(
          Geom_ImplicitGyroid(10, 1.4999).Value(gp_Pnt(.3 + .5 * x, .3 + .5 * y, .3 + .5 * z)),
          0);
  EXPECT_GT(Geom_ImplicitGyroid(10, 1.4999).Value(gp_Pnt(1.25, 1.25, 1.25)), 0);
  const auto missing = BRepClass3d_ImplicitBox::CheckTopology(hidden, 4);
  EXPECT_FALSE(missing.HomotopyResolved);
  EXPECT_FALSE(missing.ReferenceManifold);
  const auto face    = topologyField(gp_Pnt(-.31, -.37, 1.25), gp_Pnt(2.93, 2.91, 4), 1.45);
  const auto limited = BRepClass3d_ImplicitBox::CheckTopology(face, 16, 10, 1);
  EXPECT_FALSE(limited.HomotopyResolved);
  EXPECT_EQ(limited.Reason, "homotopy_cell_budget_exhausted");
  const auto nodal = topologyField(gp_Pnt(0, 0, 0), gp_Pnt(2, 2, 2), 0);
  const auto zero  = BRepClass3d_ImplicitBox::CheckTopology(nodal);
  EXPECT_FALSE(zero.HomotopyResolved);
  EXPECT_EQ(zero.Reason, "reference_vertex_sign_unresolved");
  EXPECT_THROW(BRepClass3d_ImplicitBox::CheckTopology(face, 3), Standard_DomainError);
  EXPECT_THROW(BRepClass3d_ImplicitBox::CheckTopology(face, 16, -1), Standard_DomainError);
  EXPECT_THROW(BRepClass3d_ImplicitBox::CheckTopology(face, 16, 10, 0), Standard_DomainError);
}

TEST(BRepClass3d_ImplicitBox, WholeBoxHomotopyClosedSphereAndEmptyReference)
{
  for (int resolution : {16, 32})
  {
    const auto sphere = topologyField(gp_Pnt(.3, .3, .3), gp_Pnt(2.3, 2.3, 2.3), 1.45);
    const auto result = BRepClass3d_ImplicitBox::CheckTopology(sphere, resolution);
    ASSERT_TRUE(result.HomotopyResolved) << result.Reason << " " << result.UnresolvedDimension;
    ASSERT_TRUE(result.ReferenceManifold);
    EXPECT_EQ(result.ComponentEuler, std::vector<int>({2}));
    EXPECT_TRUE(result.BoundaryLoops.empty());
  }
  const auto empty  = topologyField(gp_Pnt(.3, .3, .3), gp_Pnt(2.3, 2.3, 2.3), 2);
  const auto result = BRepClass3d_ImplicitBox::CheckTopology(empty, 4);
  ASSERT_TRUE(result.HomotopyResolved) << result.Reason;
  EXPECT_TRUE(result.ReferenceManifold);
  EXPECT_TRUE(result.Nodes.empty());
  EXPECT_TRUE(result.ComponentEuler.empty());
}

TEST(BRepClass3d_ImplicitBox, WholeBoxHomotopyScalingAndSnapshotImmutability)
{
  const auto source =
    BRepClass3d_ImplicitBox::Build(field(), box(), BRepClass3d_ImplicitBox::Operation::Common);
  std::ostringstream before;
  BRepClass3d_ImplicitBoxIO::Write(source, before);
  const auto baseline = BRepClass3d_ImplicitBox::CheckTopology(source.Faces.back());
  ASSERT_TRUE(baseline.HomotopyResolved) << baseline.Reason;
  ASSERT_TRUE(baseline.ReferenceManifold);
  std::ostringstream after;
  BRepClass3d_ImplicitBoxIO::Write(source, after);
  EXPECT_EQ(before.str(), after.str());
  for (double scale : {.1, 10.0})
  {
    const gp_Pnt low(-.31 * scale, -.37 * scale, 1.25 * scale),
      high(2.93 * scale, 2.91 * scale, 4 * scale);
    const auto face  = topologyField(low, high, 1.45, scale);
    const auto other = BRepClass3d_ImplicitBox::CheckTopology(face);
    ASSERT_TRUE(other.HomotopyResolved) << scale << " " << other.Reason;
    ASSERT_TRUE(other.ReferenceManifold);
    EXPECT_EQ(other.ComponentEuler, baseline.ComponentEuler);
    EXPECT_EQ(other.BoundaryLoops.size(), baseline.BoundaryLoops.size());
    EXPECT_NEAR(other.MaximumSimplexDiameter, baseline.MaximumSimplexDiameter * scale, 1.e-12);
  }
  auto    moved = source.Faces.back();
  gp_Trsf transform;
  transform.SetTranslation(gp_Vec(1, 0, 0));
  moved.Location(TopLoc_Location(transform));
  EXPECT_THROW(BRepClass3d_ImplicitBox::CheckTopology(moved), Standard_NotImplemented);
  EXPECT_THROW(BRepClass3d_ImplicitBox::CheckTopology(TopoDS_Face()), Standard_NotImplemented);
}

TEST(BRepClass3d_ImplicitBox, WholeBoxHomotopyRefusesCriticalBoxStrata)
{
  const double k = 2 * std::acos(-1.) / 10, x = .4 / k, y = .6 / k;
  const double z   = std::atan2(std::cos(.4), std::sin(.6)) / k;
  const double iso = std::sin(.4) * std::cos(.6) + std::hypot(std::sin(.6), std::cos(.4));
  for (const auto& face :
       {topologyField(gp_Pnt(.7, .8, .9), gp_Pnt(1.8, 1.9, 2), 1.5),
        topologyField(gp_Pnt(-.5, 2, 1.8), gp_Pnt(.5, 3, 2.5), 1),
        topologyField(gp_Pnt(x, y, z - .05), gp_Pnt(x + .05, y + .05, z + .05), iso),
        topologyField(gp_Pnt(0, 0, 0), gp_Pnt(.2, .3, .4), 0)})
  {
    const auto result = BRepClass3d_ImplicitBox::CheckTopology(face);
    EXPECT_FALSE(result.HomotopyResolved);
    EXPECT_FALSE(result.ReferenceManifold);
    EXPECT_FALSE(result.Reason.empty());
  }
}

namespace
{
double regionVolume(const BRepClass3d_ImplicitBox::RegionTopologyWitness& region)
{
  double total = 0;
  for (double volume : region.ComponentSignedVolumes)
    total += volume;
  return total;
}

void checkClosedReference(const BRepClass3d_ImplicitBox::RegionTopologyWitness& region)
{
  ASSERT_TRUE(region.Surface.HomotopyResolved) << region.Reason;
  ASSERT_TRUE(region.Surface.ReferenceManifold) << region.Reason;
  ASSERT_TRUE(region.ClosedReferenceManifold) << region.Reason;
  ASSERT_EQ(region.Triangles.size(), region.TriangleSides.size());
  ASSERT_EQ(region.Triangles.size(), region.TriangleComponents.size());
  ASSERT_EQ(region.ComponentEuler.size(), region.ComponentSignedVolumes.size());
  ASSERT_EQ(region.ComponentEuler.size(), region.ComponentRegions.size());
  ASSERT_EQ(size_t(region.RegionCount), region.RegionVolumes.size());
  for (size_t i = 0; i < region.ComponentRegions.size(); ++i)
  {
    ASSERT_GE(region.ComponentRegions[i], 0);
    ASSERT_LT(region.ComponentRegions[i], region.RegionCount);
  }
  for (int r = 0; r < region.RegionCount; ++r)
  {
    double sum = 0;
    for (size_t c = 0; c < region.ComponentRegions.size(); ++c)
      if (region.ComponentRegions[c] == r)
        sum += region.ComponentSignedVolumes[c];
    EXPECT_NEAR(sum, region.RegionVolumes[size_t(r)], 1.e-10);
    EXPECT_GT(region.RegionVolumes[size_t(r)], 0);
  }
  ASSERT_EQ(region.Surface.Nodes.size(), region.Surface.GridEdges.size());
  for (size_t i = 0; i < region.Surface.Nodes.size(); ++i)
    EXPECT_EQ(region.Nodes[i].Distance(region.Surface.Nodes[i]), 0);
  std::map<std::pair<int, int>, std::pair<int, int>> edges;
  for (size_t i = 0; i < region.Triangles.size(); ++i)
  {
    const auto& t = region.Triangles[i];
    ASSERT_GE(region.TriangleComponents[i], 0);
    ASSERT_LT(size_t(region.TriangleComponents[i]), region.ComponentEuler.size());
    for (int j = 0; j < 3; ++j)
    {
      const int a = t[size_t(j)], b = t[size_t((j + 1) % 3)];
      ASSERT_GE(a, 0);
      ASSERT_LT(size_t(a), region.Nodes.size());
      auto& edge = edges[std::minmax(a, b)];
      ++edge.first;
      edge.second += a < b ? 1 : -1;
    }
    if (i < region.Surface.Triangles.size())
    {
      EXPECT_EQ(t, region.Surface.Triangles[i]);
      EXPECT_EQ(region.TriangleSides[i], -1);
    }
    else
    {
      const int side = region.TriangleSides[i];
      ASSERT_GE(side, 0);
      ASSERT_LT(side, 6);
      const gp_Vec area =
        gp_Vec(region.Nodes[size_t(t[0])], region.Nodes[size_t(t[1])])
          .Crossed(gp_Vec(region.Nodes[size_t(t[0])], region.Nodes[size_t(t[2])]));
      EXPECT_GT(area.Coord(side / 2 + 1) * (side % 2 ? 1 : -1), 0);
    }
  }
  for (const auto& edge : edges)
  {
    EXPECT_EQ(edge.second.first, 2);
    EXPECT_EQ(edge.second.second, 0);
  }
  int euler = 0;
  for (int value : region.ComponentEuler)
    euler += value;
  EXPECT_EQ(int(region.Nodes.size()) - int(edges.size()) + int(region.Triangles.size()), euler);
}
} // namespace

TEST(BRepClass3d_ImplicitBox, ClosedRegionReferenceBothOperationsAndComplementaryVolume)
{
  for (bool open : {false, true})
  {
    const gp_Pnt low     = open ? gp_Pnt(-3.83, -4.07, -2.67) : gp_Pnt(-.31, -.37, 1.25);
    const gp_Pnt high    = open ? gp_Pnt(4.17, 3.93, 3.37) : gp_Pnt(2.93, 2.91, 4);
    const auto   operand = BRepPrimAPI_MakeBox(low, high).Shape();
    double       total   = 0;
    for (auto operation : {BRepClass3d_ImplicitBox::Operation::Common,
                           BRepClass3d_ImplicitBox::Operation::BoxMinusField})
    {
      const auto source =
        BRepClass3d_ImplicitBox::Build(field(open ? 0 : 1.45), operand, operation);
      std::ostringstream before;
      BRepClass3d_ImplicitBoxIO::Write(source, before);
      const auto region = BRepClass3d_ImplicitBox::CheckRegionTopology(source.Faces.back());
      checkClosedReference(region);
      ASSERT_TRUE(region.ClosedReferenceManifold) << region.Reason;
      EXPECT_GT(region.Nodes.size(), region.Surface.Nodes.size());
      EXPECT_GT(region.Triangles.size(), region.Surface.Triangles.size());
      EXPECT_GT(regionVolume(region), 0);
      total += regionVolume(region);
      std::ostringstream after;
      BRepClass3d_ImplicitBoxIO::Write(source, after);
      EXPECT_EQ(before.str(), after.str());
    }
    EXPECT_NEAR(total, (high.X() - low.X()) * (high.Y() - low.Y()) * (high.Z() - low.Z()), 1.e-9);
  }
}

TEST(BRepClass3d_ImplicitBox, ClosedRegionReferenceDistinguishesCavityAndEmptyFullRegions)
{
  const gp_Pnt low(.3, .3, .3), high(2.3, 2.3, 2.3);
  auto         source = topologyField(low, high, 1.45);
  const auto   common = BRepClass3d_ImplicitBox::CheckRegionTopology(source, 16);
  checkClosedReference(common);
  ASSERT_EQ(common.ComponentSignedVolumes.size(), 2u);
  EXPECT_EQ(common.ComponentEuler, std::vector<int>({2, 2}));
  EXPECT_EQ(common.RegionCount, 1);
  EXPECT_EQ(common.ComponentRegions, std::vector<int>({0, 0}));
  int positive = 0, negative = 0;
  for (double volume : common.ComponentSignedVolumes)
  {
    positive += volume > 0;
    negative += volume < 0;
  }
  EXPECT_EQ(positive, 1);
  EXPECT_EQ(negative, 1);
  BRepMesh_ImplicitBoundary::Result auditInput;
  auditInput.Nodes             = common.Nodes;
  auditInput.ImplicitTriangles = common.Triangles;
  EXPECT_EQ(BRepMesh_ImplicitBoundary::Audit(auditInput).Status,
            BRepMesh_ImplicitBoundary::MeshAudit::State::Clear);
  source.Reverse();
  const auto cut = BRepClass3d_ImplicitBox::CheckRegionTopology(source, 16);
  checkClosedReference(cut);
  EXPECT_EQ(cut.ComponentEuler, std::vector<int>({2}));
  EXPECT_EQ(cut.RegionCount, 1);
  EXPECT_GT(regionVolume(cut), 0);
  EXPECT_EQ(cut.Nodes.size(), cut.Surface.Nodes.size());
  EXPECT_EQ(cut.Triangles.size(), cut.Surface.Triangles.size());
  EXPECT_NEAR(regionVolume(common) + regionVolume(cut), 8, 1.e-12);
  for (double iso : {-2.0, 2.0})
    for (bool reversed : {false, true})
    {
      auto face = topologyField(low, high, iso);
      if (reversed)
        face.Reverse();
      const auto result = BRepClass3d_ImplicitBox::CheckRegionTopology(face, 4);
      checkClosedReference(result);
      const bool full = (iso > 0) != reversed;
      EXPECT_EQ(result.RegionCount, full ? 1 : 0);
      EXPECT_EQ(result.ComponentEuler, full ? std::vector<int>({2}) : std::vector<int>());
      EXPECT_NEAR(regionVolume(result), full ? 8 : 0, 1.e-12);
      EXPECT_TRUE(result.Surface.Triangles.empty());
    }
}

TEST(BRepClass3d_ImplicitBox, ClosedRegionReferenceMultipleComponentsAndFailureLimits)
{
  auto face = topologyField(gp_Pnt(.3, .3, .3), gp_Pnt(12.3, 2.3, 2.3), 1.45);
  face.Reverse();
  const auto islands = BRepClass3d_ImplicitBox::CheckRegionTopology(face, 64);
  checkClosedReference(islands);
  ASSERT_TRUE(islands.ClosedReferenceManifold) << islands.Reason;
  EXPECT_EQ(islands.ComponentEuler, std::vector<int>({2, 2}));
  EXPECT_EQ(islands.RegionCount, 2);
  ASSERT_EQ(islands.ComponentRegions.size(), 2u);
  EXPECT_NE(islands.ComponentRegions[0], islands.ComponentRegions[1]);
  for (double volume : islands.ComponentSignedVolumes)
    EXPECT_GT(volume, 0);
  const auto limited = BRepClass3d_ImplicitBox::CheckRegionTopology(face, 32, 16, 1);
  EXPECT_FALSE(limited.Surface.HomotopyResolved);
  EXPECT_FALSE(limited.ClosedReferenceManifold);
  EXPECT_TRUE(limited.Nodes.empty());
  EXPECT_FALSE(limited.Reason.empty());
  const auto critical = BRepClass3d_ImplicitBox::CheckRegionTopology(
    topologyField(gp_Pnt(.7, .8, .9), gp_Pnt(1.8, 1.9, 2), 1.5));
  EXPECT_FALSE(critical.ClosedReferenceManifold);
  EXPECT_THROW(BRepClass3d_ImplicitBox::CheckRegionTopology(face, 3), Standard_DomainError);
}

TEST(BRepClass3d_ImplicitBox, StoredCapsAttachToReferenceShellsAndRegionsBothOperations)
{
  for (bool open : {false, true})
    for (auto operation : {BRepClass3d_ImplicitBox::Operation::Common,
                           BRepClass3d_ImplicitBox::Operation::BoxMinusField})
    {
      const auto operand =
        open ? BRepPrimAPI_MakeBox(gp_Pnt(-3.83, -4.07, -2.67), gp_Pnt(4.17, 3.93, 3.37)).Shape()
             : box();
      const auto source =
        BRepClass3d_ImplicitBox::Build(field(open ? 0 : 1.45), operand, operation);
      std::ostringstream before;
      BRepClass3d_ImplicitBoxIO::Write(source, before);
      const auto report = BRepClass3d_ImplicitBox::CheckCapAttachment(source);
      ASSERT_TRUE(report.Resolved) << open << " " << int(operation) << " " << report.Reason;
      EXPECT_TRUE(report.StoredGraphResolved);
      EXPECT_TRUE(report.RegionReferenceResolved);
      EXPECT_EQ(report.ReferenceCaps.size(), source.Faces.size() - 1);
      EXPECT_EQ(report.ReferenceCapCount, int(source.Faces.size()) - 1);
      std::set<int> unique;
      const auto&   gyroid = *static_cast<const Geom_ImplicitGyroid*>(
        occ::down_cast<BRep_TImplicitFace>(source.Faces.back().TShape())->Field().get());
      for (size_t cap = 0; cap < report.ReferenceCaps.size(); ++cap)
      {
        EXPECT_TRUE(unique.insert(report.ReferenceCaps[cap]).second);
        EXPECT_GT(report.Clearances[cap], report.RequiredClearance);
        EXPECT_GE(report.ReferenceShells[cap], 0);
        EXPECT_LT(report.ReferenceShells[cap], report.ReferenceShellCount);
        EXPECT_GE(report.ReferenceRegions[cap], 0);
        EXPECT_LT(report.ReferenceRegions[cap], report.ReferenceRegionCount);
        EXPECT_EQ(BRepClass_FaceClassifier(source.Faces[cap], report.Witnesses[cap], 1.e-9).State(),
                  TopAbs_IN);
        EXPECT_GT(gyroid.Value(report.Witnesses[cap])
                    * (operation == BRepClass3d_ImplicitBox::Operation::Common ? -1 : 1),
                  0);
      }
      std::ostringstream after;
      BRepClass3d_ImplicitBoxIO::Write(source, after);
      EXPECT_EQ(before.str(), after.str());
      std::istringstream                input(before.str());
      BRepClass3d_ImplicitBox::Boundary loaded;
      BRepClass3d_ImplicitBoxIO::Read(input, loaded);
      EXPECT_FALSE(loaded.SeamBoundVerified);
      EXPECT_TRUE(BRepClass3d_ImplicitBox::CheckCapAttachment(loaded).Resolved);
      EXPECT_FALSE(loaded.SeamBoundVerified);
    }
}

TEST(BRepClass3d_ImplicitBox, CapAttachmentMultipleRegionsReorderingAndLimits)
{
  auto source = BRepClass3d_ImplicitBox::Build(
    field(),
    BRepPrimAPI_MakeBox(gp_Pnt(-.31, -.37, 1.25), gp_Pnt(12.93, 2.91, 4)).Shape(),
    BRepClass3d_ImplicitBox::Operation::BoxMinusField);
  const auto report = BRepClass3d_ImplicitBox::CheckCapAttachment(source, 64);
  ASSERT_TRUE(report.Resolved) << report.Reason;
  EXPECT_EQ(report.ReferenceRegionCount, 2);
  ASSERT_EQ(report.ReferenceCaps.size(), 2u);
  EXPECT_NE(report.ReferenceRegions[0], report.ReferenceRegions[1]);
  std::swap(source.Faces[0], source.Faces[1]);
  std::swap(source.Sources[0], source.Sources[1]);
  BRep_Builder builder;
  builder.MakeCompound(source.Shape);
  for (const auto& face : source.Faces)
    builder.Add(source.Shape, face);
  const auto reordered = BRepClass3d_ImplicitBox::CheckCapAttachment(source, 64);
  ASSERT_TRUE(reordered.Resolved) << reordered.Reason;
  EXPECT_EQ(reordered.ReferenceCaps[0], report.ReferenceCaps[1]);
  EXPECT_EQ(reordered.ReferenceCaps[1], report.ReferenceCaps[0]);
  const auto limited = BRepClass3d_ImplicitBox::CheckCapAttachment(source, 64, 16, 4000000, 1);
  EXPECT_FALSE(limited.Resolved);
  EXPECT_FALSE(limited.Reason.empty());
  EXPECT_THROW(BRepClass3d_ImplicitBox::CheckCapAttachment(source, 3), Standard_DomainError);
}

TEST(BRepClass3d_ImplicitBox, DualCapTrimAttachmentBothOperationsAndDirectedWires)
{
  for (bool open : {false, true})
    for (auto operation : {BRepClass3d_ImplicitBox::Operation::Common,
                           BRepClass3d_ImplicitBox::Operation::BoxMinusField})
    {
      const auto operand =
        open ? BRepPrimAPI_MakeBox(gp_Pnt(-3.83, -4.07, -2.67), gp_Pnt(4.17, 3.93, 3.37)).Shape()
             : box();
      const auto source =
        BRepClass3d_ImplicitBox::Build(field(open ? 0 : 1.45), operand, operation);
      std::ostringstream before;
      BRepClass3d_ImplicitBoxIO::Write(source, before);
      const auto report = BRepClass3d_ImplicitBox::CheckTrimAttachment(source);
      ASSERT_TRUE(report.Resolved) << open << " " << int(operation) << " " << report.Reason;
      EXPECT_TRUE(report.CapRegionsAttached);
      EXPECT_TRUE(report.ComplementCapsAttached);
      EXPECT_EQ(report.ReferenceArcCount, source.SeamCount);
      EXPECT_EQ(report.SeamReferenceArcs.size(), size_t(source.SeamCount));
      EXPECT_EQ(report.ReferenceWireCount, open ? 3 : 1);
      EXPECT_EQ(report.ReferenceSurfaceComponentCount, 1);
      EXPECT_EQ(report.WireReferenceLoops.size(), size_t(report.ReferenceWireCount));
      std::set<int> arcs(report.SeamReferenceArcs.begin(), report.SeamReferenceArcs.end());
      EXPECT_EQ(arcs.size(), report.SeamReferenceArcs.size());
      std::set<int> loops(report.WireReferenceLoops.begin(), report.WireReferenceLoops.end());
      EXPECT_EQ(loops.size(), report.WireReferenceLoops.size());
      for (int component : report.WireSurfaceComponents)
        EXPECT_EQ(component, 0);
      std::ostringstream after;
      BRepClass3d_ImplicitBoxIO::Write(source, after);
      EXPECT_EQ(before.str(), after.str());
    }
}

TEST(BRepClass3d_ImplicitBox, DualCapTrimAttachmentDistinguishesHolesAndSurfaceComponents)
{
  const auto operand =
    BRepPrimAPI_MakeBox(gp_Pnt(-.31, -.37, 1.25), gp_Pnt(12.93, 2.91, 4)).Shape();
  for (auto operation : {BRepClass3d_ImplicitBox::Operation::Common,
                         BRepClass3d_ImplicitBox::Operation::BoxMinusField})
  {
    const auto source = BRepClass3d_ImplicitBox::Build(field(), operand, operation);
    const auto report = BRepClass3d_ImplicitBox::CheckTrimAttachment(source, 64);
    ASSERT_TRUE(report.Resolved) << report.Reason;
    EXPECT_EQ(report.ReferenceArcCount, 2);
    EXPECT_EQ(report.ReferenceWireCount, 2);
    EXPECT_EQ(report.ReferenceSurfaceComponentCount, 2);
    ASSERT_EQ(report.SeamReferenceArcs.size(), 2u);
    EXPECT_NE(report.SeamReferenceArcs[0], report.SeamReferenceArcs[1]);
    ASSERT_EQ(report.WireSurfaceComponents.size(), 2u);
    EXPECT_NE(report.WireSurfaceComponents[0], report.WireSurfaceComponents[1]);
  }
}

TEST(BRepClass3d_ImplicitBox, DualCapTrimAttachmentRestoreAndBudgetFailures)
{
  const auto source =
    BRepClass3d_ImplicitBox::Build(field(), box(), BRepClass3d_ImplicitBox::Operation::Common);
  std::ostringstream output;
  BRepClass3d_ImplicitBoxIO::Write(source, output);
  std::istringstream                input(output.str());
  BRepClass3d_ImplicitBox::Boundary loaded;
  BRepClass3d_ImplicitBoxIO::Read(input, loaded);
  ASSERT_TRUE(BRepClass3d_ImplicitBox::CheckTrimAttachment(loaded).Resolved);
  EXPECT_FALSE(loaded.SeamBoundVerified);
  const auto limited = BRepClass3d_ImplicitBox::CheckTrimAttachment(loaded, 32, 16, 4000000, 1);
  EXPECT_FALSE(limited.Resolved);
  EXPECT_FALSE(limited.Reason.empty());
  const auto low = BRepClass3d_ImplicitBox::CheckTrimAttachment(loaded, 32, 16, 1);
  EXPECT_FALSE(low.Resolved);
  EXPECT_FALSE(low.CapRegionsAttached);
  EXPECT_THROW(BRepClass3d_ImplicitBox::CheckTrimAttachment(loaded, 3), Standard_DomainError);
}

TEST(BRepClass3d_ImplicitBox, GenuineNativeShellsBothOperationsAndFieldFaces)
{
  for (bool open : {false, true})
    for (auto operation : {BRepClass3d_ImplicitBox::Operation::Common,
                           BRepClass3d_ImplicitBox::Operation::BoxMinusField})
    {
      const auto operand =
        open ? BRepPrimAPI_MakeBox(gp_Pnt(-3.83, -4.07, -2.67), gp_Pnt(4.17, 3.93, 3.37)).Shape()
             : box();
      const auto source =
        BRepClass3d_ImplicitBox::Build(field(open ? 0 : 1.45), operand, operation);
      std::ostringstream before;
      BRepClass3d_ImplicitBoxIO::Write(source, before);
      const auto result = BRepClass3d_ImplicitBox::BuildNativeShells(source);
      ASSERT_TRUE(result.FacesAttached) << result.Reason;
      ASSERT_TRUE(result.ShellsConstructed)
        << open << " " << int(operation) << " " << result.Reason;
      EXPECT_TRUE(result.TopologyValidated);
      EXPECT_EQ(result.Faces.size(), source.Faces.size());
      EXPECT_EQ(result.ComponentDomains.size(), 1u);
      EXPECT_EQ(result.Separators, 0);
      int implicit = 0;
      for (size_t i = 0; i < result.Faces.size(); ++i)
        if (BRep_Tool::IsImplicit(result.Faces[i]))
        {
          ++implicit;
          EXPECT_EQ(result.FaceSources[i], int(source.Faces.size()) - 1);
          EXPECT_EQ(result.Faces[i].Orientation(), source.Faces.back().Orientation());
          EXPECT_EQ(BRep_ImplicitTool::Value(result.Faces[i], gp_Pnt(1, 2, 3)),
                    BRep_ImplicitTool::Value(source.Faces.back(), gp_Pnt(1, 2, 3)));
          TopLoc_Location location;
          EXPECT_TRUE(BRep_Tool::Triangulation(result.Faces[i], location).IsNull());
        }
      EXPECT_EQ(implicit, 1);
      for (const auto& shell : result.Shells)
      {
        EXPECT_EQ(shell.ShapeType(), TopAbs_SHELL);
        EXPECT_TRUE(shell.Closed());
      }
      std::ostringstream after;
      BRepClass3d_ImplicitBoxIO::Write(source, after);
      EXPECT_EQ(before.str(), after.str());
    }
}

TEST(BRepClass3d_ImplicitBox, NativeShellsSplitDisconnectedComponentsAtCheckedEmptyPlanes)
{
  for (auto operation : {BRepClass3d_ImplicitBox::Operation::Common,
                         BRepClass3d_ImplicitBox::Operation::BoxMinusField})
  {
    const auto source = BRepClass3d_ImplicitBox::Build(
      field(),
      BRepPrimAPI_MakeBox(gp_Pnt(-.31, -.37, 1.25), gp_Pnt(12.93, 2.91, 4)).Shape(),
      operation);
    const auto result = BRepClass3d_ImplicitBox::BuildNativeShells(source, 64);
    ASSERT_TRUE(result.ShellsConstructed) << result.Reason;
    EXPECT_TRUE(result.TopologyValidated);
    EXPECT_EQ(result.Separators, 1);
    EXPECT_GT(result.SeparatorCells, 1u);
    ASSERT_EQ(result.ComponentDomains.size(), 2u);
    auto domains = result.ComponentDomains;
    std::sort(domains.begin(), domains.end(), [](const auto& a, const auto& b) {
      return a[0].X() < b[0].X();
    });
    EXPECT_EQ(domains[0][1].X(), domains[1][0].X());
    EXPECT_GT(domains[0][1].X(), 3);
    EXPECT_LT(domains[0][1].X(), 10);
    int implicit = 0;
    for (const auto& face : result.Faces)
      implicit += BRep_Tool::IsImplicit(face);
    EXPECT_EQ(implicit, 2);
    EXPECT_EQ(result.Shells.size(),
              operation == BRepClass3d_ImplicitBox::Operation::Common ? 1u : 2u);
    if (operation == BRepClass3d_ImplicitBox::Operation::Common)
    {
      const auto limited =
        BRepClass3d_ImplicitBox::BuildNativeShells(source, 64, 16, 4000000, 2000000, 1);
      EXPECT_FALSE(limited.ShellsConstructed);
      EXPECT_TRUE(limited.Shape.IsNull());
      EXPECT_TRUE(limited.Shells.empty());
      EXPECT_FALSE(limited.Reason.empty());
    }
  }
}

TEST(BRepClass3d_ImplicitBox, NativeShellsEnclosedCavityAndEmptyFullRegions)
{
  const auto operand = BRepPrimAPI_MakeBox(gp_Pnt(.3, .3, .3), gp_Pnt(2.3, 2.3, 2.3)).Shape();
  for (auto operation : {BRepClass3d_ImplicitBox::Operation::Common,
                         BRepClass3d_ImplicitBox::Operation::BoxMinusField})
  {
    const auto source = BRepClass3d_ImplicitBox::Build(field(), operand, operation);
    const auto result = BRepClass3d_ImplicitBox::BuildNativeShells(source);
    ASSERT_TRUE(result.ShellsConstructed) << result.Reason;
    EXPECT_EQ(result.ComponentDomains.size(), 1u);
    EXPECT_EQ(result.Shells.size(),
              operation == BRepClass3d_ImplicitBox::Operation::Common ? 2u : 1u);
    for (int region : result.ShellRegions)
      EXPECT_EQ(region, 0);
    EXPECT_EQ(result.Source.BoxEdgeRootCount, 0);
  }
  for (auto operation : {BRepClass3d_ImplicitBox::Operation::Common,
                         BRepClass3d_ImplicitBox::Operation::BoxMinusField})
  {
    const auto source = BRepClass3d_ImplicitBox::Build(field(2), operand, operation);
    const auto result = BRepClass3d_ImplicitBox::BuildNativeShells(source, 4);
    ASSERT_TRUE(result.ShellsConstructed) << result.Reason;
    EXPECT_TRUE(result.ComponentDomains.empty());
    EXPECT_EQ(result.Shells.size(),
              operation == BRepClass3d_ImplicitBox::Operation::Common ? 1u : 0u);
    for (const auto& face : result.Faces)
      EXPECT_FALSE(BRep_Tool::IsImplicit(face));
    for (const auto& shell : result.Shells)
      EXPECT_TRUE(BRepCheck_Analyzer(shell).IsValid());
  }
}

TEST(BRepClass3d_ImplicitBox, NativeSolidsGatedByWholeCapGeometryBothOperations)
{
  for (bool open : {false, true})
    for (auto operation : {BRepClass3d_ImplicitBox::Operation::Common,
                           BRepClass3d_ImplicitBox::Operation::BoxMinusField})
    {
      const auto operand =
        open ? BRepPrimAPI_MakeBox(gp_Pnt(-3.83, -4.07, -2.67), gp_Pnt(4.17, 3.93, 3.37)).Shape()
             : box();
      const auto source =
        BRepClass3d_ImplicitBox::Build(field(open ? 0 : 1.45), operand, operation);
      std::ostringstream before;
      BRepClass3d_ImplicitBoxIO::Write(source, before);
      const auto result = BRepClass3d_ImplicitBox::BuildNativeSolids(source);
      ASSERT_TRUE(result.Constructed)
        << open << " " << int(operation) << " " << result.Reason << " cells "
        << result.Geometry.Cells[0] << "," << result.Geometry.Cells[1] << " segments "
        << result.Geometry.SegmentTests;
      EXPECT_TRUE(result.Validated);
      EXPECT_TRUE(result.Geometry.Resolved);
      EXPECT_LE(result.Geometry.BoundaryBound, .1);
      EXPECT_GT(result.Geometry.Cells[0], 0u);
      EXPECT_GT(result.Geometry.Cells[1], 0u);
      EXPECT_LE(result.MaximumStoredTolerance, .1);
      EXPECT_EQ(result.Solids.size(), 1u);
      for (const auto& solid : result.Solids)
      {
        EXPECT_EQ(solid.ShapeType(), TopAbs_SOLID);
        EXPECT_TRUE(solid.Closed());
      }
      std::ostringstream after;
      BRepClass3d_ImplicitBoxIO::Write(source, after);
      EXPECT_EQ(before.str(), after.str());
    }
}

TEST(BRepClass3d_ImplicitBox, CapGeometryRejectsMissingAndExtraneousNativeCaps)
{
  using Ops                             = BRepClass3d_ImplicitBox;
  const auto                     source = Ops::Build(field(), box(), Ops::Operation::BoxMinusField);
  const std::vector<TopoDS_Face> caps(source.Faces.begin(), source.Faces.end() - 1);
  const auto                     good = Ops::CheckCapGeometry(source, caps);
  ASSERT_TRUE(good.Resolved) << good.Reason;
  const auto missing = Ops::CheckCapGeometry(source, {});
  EXPECT_FALSE(missing.Resolved);
  auto         ghost = source;
  const auto   old   = source.Faces.back();
  const auto   data  = occ::down_cast<BRep_TImplicitFace>(old.TShape());
  TopoDS_Face  changed;
  BRep_Builder builder;
  builder.MakeImplicitFace(changed,
                           new Geom_ImplicitGyroid(10, 2),
                           data->Minimum(),
                           data->Maximum(),
                           TopLoc_Location(),
                           .1);
  for (TopoDS_Iterator it(old, false, false); it.More(); it.Next())
    builder.Add(changed, it.Value());
  changed.Orientation(old.Orientation());
  ghost.Faces.back() = changed;
  builder.MakeCompound(ghost.Shape);
  for (const auto& face : ghost.Faces)
    builder.Add(ghost.Shape, face);
  ghost.GeometricBound    = 0;
  ghost.SeamBoundVerified = true;
  const auto extra        = Ops::CheckCapGeometry(ghost, caps);
  EXPECT_FALSE(extra.Resolved);
  EXPECT_EQ(extra.UnresolvedDirection, 1);
  const auto limited = Ops::CheckCapGeometry(source, caps, .1, 16, 1);
  EXPECT_FALSE(limited.Resolved);
  const auto segmentLimit = Ops::CheckCapGeometry(source, caps, .1, 16, 2000000, 1);
  EXPECT_FALSE(segmentLimit.Resolved);
  EXPECT_EQ(segmentLimit.SegmentTests, 1u);
  EXPECT_EQ(segmentLimit.Reason, "Cap geometry segment budget exhausted");
  EXPECT_THROW(Ops::CheckCapGeometry(source, caps, 0), Standard_DomainError);
}

TEST(BRepClass3d_ImplicitBox, NativeSolidsFullAndEmptyRetainMathematicalRegion)
{
  using Ops = BRepClass3d_ImplicitBox;
  for (auto operation : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
  {
    const auto source = Ops::Build(field(2), box(), operation);
    const auto result = Ops::BuildNativeSolids(source, .1, 4);
    ASSERT_TRUE(result.Validated) << result.Reason;
    EXPECT_TRUE(result.Geometry.Resolved);
    EXPECT_LE(result.Geometry.BoundaryBound, .1);
    EXPECT_EQ(result.Solids.size(), operation == Ops::Operation::Common ? 1u : 0u);
    EXPECT_FALSE(result.Shape.IsNull());
    for (const auto& solid : result.Solids)
      EXPECT_TRUE(BRepCheck_Analyzer(solid).IsValid()); // Only conventional box faces here.
  }
}

TEST(BRepClass3d_ImplicitBox, NativeSolidsSeparateIslandsAndGroupCavityShells)
{
  using Ops = BRepClass3d_ImplicitBox;
  const auto islands =
    Ops::Build(field(),
               BRepPrimAPI_MakeBox(gp_Pnt(-.31, -.37, 1.25), gp_Pnt(12.93, 2.91, 4)).Shape(),
               Ops::Operation::BoxMinusField);
  const auto result = Ops::BuildNativeSolids(islands, .1, 64);
  ASSERT_TRUE(result.Validated) << result.Reason;
  EXPECT_EQ(result.Solids.size(), 2u);
  const auto cavity =
    Ops::Build(field(),
               BRepPrimAPI_MakeBox(gp_Pnt(.3, .3, .3), gp_Pnt(2.3, 2.3, 2.3)).Shape(),
               Ops::Operation::Common);
  const auto hollow = Ops::BuildNativeSolids(cavity);
  ASSERT_TRUE(hollow.Validated) << hollow.Reason;
  ASSERT_EQ(hollow.Solids.size(), 1u);
  int shells = 0;
  for (TopoDS_Iterator it(hollow.Solids[0]); it.More(); it.Next())
    ++shells;
  EXPECT_EQ(shells, 2);
  const auto tight = Ops::BuildNativeSolids(islands, .01, 64);
  EXPECT_FALSE(tight.Constructed);
  EXPECT_FALSE(tight.Validated);
  EXPECT_TRUE(tight.Shape.IsNull());
  const auto budget = Ops::BuildNativeSolids(cavity, .1, 32, 16, 4000000, 2000000, 2000000, 1);
  EXPECT_FALSE(budget.Constructed);
  EXPECT_TRUE(budget.Shape.IsNull());
}

TEST(BRepClass3d_ImplicitBox, NativeSolidsScaledPeriodicTranslationAtFixedTolerance)
{
  using Ops = BRepClass3d_ImplicitBox;
  for (double scale : {.5, 2.})
  {
    // Integer-period coordinate translations preserve the gyroid equation.
    auto point = [scale](double x, double y, double z) {
      return gp_Pnt((x + 20) * scale, (y - 10) * scale, (z + 30) * scale);
    };
    TopoDS_Face face;
    BRep_Builder().MakeImplicitFace(face,
                                    new Geom_ImplicitGyroid(10 * scale, 1.45),
                                    point(-5, -5, -5),
                                    point(16, 16, 6),
                                    TopLoc_Location(),
                                    .1);
    const auto operand = BRepPrimAPI_MakeBox(point(-.31, -.37, 1.25), point(2.93, 2.91, 4)).Shape();
    for (auto operation : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
    {
      const auto source = Ops::Build(face, operand, operation, .1);
      const auto result = Ops::BuildNativeSolids(source, .1);
      ASSERT_TRUE(result.Validated) << scale << " " << int(operation) << " " << result.Reason;
      EXPECT_EQ(result.Solids.size(), 1u);
      EXPECT_LE(result.Geometry.BoundaryBound, .1);
      EXPECT_LE(result.MaximumStoredTolerance, .1);
    }
  }
}

TEST(BRepClass3d_ImplicitBox, NativeSolidTopologyIORestoresActualSharedGraphNotValidation)
{
  using Ops = BRepClass3d_ImplicitBox;
  for (int fixture = 0; fixture < 6; ++fixture)
  {
    const bool open = fixture == 2 || fixture == 3;
    const auto operand =
      fixture == 4   ? BRepPrimAPI_MakeBox(gp_Pnt(-.31, -.37, 1.25), gp_Pnt(12.93, 2.91, 4)).Shape()
      : fixture == 5 ? BRepPrimAPI_MakeBox(gp_Pnt(.3, .3, .3), gp_Pnt(2.3, 2.3, 2.3)).Shape()
      : open ? BRepPrimAPI_MakeBox(gp_Pnt(-3.83, -4.07, -2.67), gp_Pnt(4.17, 3.93, 3.37)).Shape()
             : box();
    const auto operation = fixture == 1 || fixture == 3 || fixture == 4
                             ? Ops::Operation::BoxMinusField
                             : Ops::Operation::Common;
    const auto source    = Ops::Build(field(open ? 0 : 1.45), operand, operation);
    const auto value     = Ops::BuildNativeSolids(source, .1, fixture == 4 ? 64 : 32);
    ASSERT_TRUE(value.Validated) << fixture << " " << value.Reason;
    if (fixture == 0)
      BRepMesh_IncrementalMesh(value.Boundary.Faces[0], .05);
    std::vector<occ::handle<Poly_Triangulation>> caches;
    for (const auto& face : value.Boundary.Faces)
    {
      TopLoc_Location location;
      caches.push_back(BRep_Tool::Triangulation(face, location));
    }
    if (fixture == 0)
      ASSERT_FALSE(caches[0].IsNull());
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> originalParts;
    TopExp::MapShapes(value.Shape, originalParts);
    std::vector<std::array<bool, 3>> flags;
    for (const auto& part : originalParts)
      flags.push_back({part.Closed(), part.Checked(), part.Free()});
    std::ostringstream before;
    BRepClass3d_ImplicitBoxIO::Write(value.Boundary.Source, before);
    std::ostringstream record;
    BRepClass3d_ImplicitSolidIO::Write(value, record);
    Ops::NativeSolids  loaded;
    std::istringstream input(record.str());
    BRepClass3d_ImplicitSolidIO::Read(input, loaded);
    ASSERT_TRUE(loaded.Constructed);
    EXPECT_FALSE(loaded.Validated);
    EXPECT_FALSE(loaded.Geometry.Resolved);
    EXPECT_FALSE(loaded.Boundary.TopologyValidated);
    EXPECT_FALSE(loaded.Boundary.FacesAttached);
    EXPECT_FALSE(loaded.Boundary.Source.SeamBoundVerified);
    EXPECT_EQ(loaded.Reason, "restored_native_solids_unverified");
    EXPECT_EQ(loaded.Solids.size(), value.Solids.size());
    EXPECT_EQ(loaded.Boundary.Shells.size(), value.Boundary.Shells.size());
    EXPECT_EQ(loaded.Boundary.FaceSources, value.Boundary.FaceSources);
    EXPECT_EQ(loaded.Boundary.FaceComponents, value.Boundary.FaceComponents);
    EXPECT_EQ(loaded.Boundary.FaceShells, value.Boundary.FaceShells);
    EXPECT_EQ(loaded.Boundary.ShellRegions, value.Boundary.ShellRegions);
    EXPECT_EQ(loaded.Boundary.ShellEuler, value.Boundary.ShellEuler);
    ASSERT_EQ(loaded.Boundary.ComponentDomains.size(), value.Boundary.ComponentDomains.size());
    for (size_t c = 0; c < value.Boundary.ComponentDomains.size(); ++c)
      for (size_t end = 0; end < 2; ++end)
        for (int axis = 1; axis <= 3; ++axis)
          EXPECT_DOUBLE_EQ(loaded.Boundary.ComponentDomains[c][end].Coord(axis),
                           value.Boundary.ComponentDomains[c][end].Coord(axis));
    for (const auto& solid : loaded.Solids)
    {
      EXPECT_EQ(solid.ShapeType(), TopAbs_SOLID);
      EXPECT_FALSE(solid.Closed());
      EXPECT_FALSE(solid.Checked());
    }
    for (const auto& shell : loaded.Boundary.Shells)
      EXPECT_FALSE(shell.Closed());
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> edges;
    TopExp::MapShapes(loaded.Shape, TopAbs_EDGE, edges);
    for (const auto& edge : edges)
    {
      int uses = 0, balance = 0;
      for (const auto& face : loaded.Boundary.Faces)
        for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
          if (it.Current().IsSame(edge))
          {
            ++uses;
            balance += it.Current().Orientation() == TopAbs_FORWARD ? 1 : -1;
          }
      EXPECT_EQ(uses, 2);
      EXPECT_EQ(balance, 0);
    }
    std::ostringstream saved, after;
    BRepClass3d_ImplicitSolidIO::Write(loaded, saved);
    EXPECT_EQ(saved.str(), record.str());
    BRepClass3d_ImplicitBoxIO::Write(value.Boundary.Source, after);
    EXPECT_EQ(before.str(), after.str());
    if (fixture == 4)
      EXPECT_EQ(loaded.Solids.size(), 2u);
    if (fixture == 5)
    {
      EXPECT_EQ(loaded.Solids.size(), 1u);
      EXPECT_EQ(loaded.Boundary.Shells.size(), 2u);
    }
    size_t partIndex = 0;
    for (const auto& part : originalParts)
    {
      EXPECT_EQ(flags[partIndex],
                (std::array<bool, 3>{part.Closed(), part.Checked(), part.Free()}));
      ++partIndex;
    }
    for (size_t f = 0; f < value.Boundary.Faces.size(); ++f)
    {
      TopLoc_Location location;
      EXPECT_EQ(caches[f], BRep_Tool::Triangulation(value.Boundary.Faces[f], location));
      EXPECT_TRUE(BRep_Tool::Triangulation(loaded.Boundary.Faces[f], location).IsNull());
    }
  }
}

TEST(BRepClass3d_ImplicitBox, NativeSolidTopologyIORejectsCorruptionAtomically)
{
  using Ops = BRepClass3d_ImplicitBox;
  const auto value =
    Ops::BuildNativeSolids(Ops::Build(field(), box(), Ops::Operation::BoxMinusField));
  ASSERT_TRUE(value.Validated) << value.Reason;
  std::ostringstream output;
  BRepClass3d_ImplicitSolidIO::Write(value, output);
  const auto               valid = output.str();
  std::vector<std::string> invalid{valid.substr(0, valid.size() - 1),
                                   valid + "garbage",
                                   "OCCT_IMPLICIT_SOLIDS 99 1 0 0 0 0 .1\n",
                                   "OCCT_IMPLICIT_SOLIDS 1 16777217 0 0 0 0 .1\n"};
  auto                     wire     = valid;
  const std::string        wireTail = " 1 0\nSHELLS\n";
  const auto               at       = wire.find(wireTail);
  ASSERT_NE(at, std::string::npos);
  wire.replace(at, wireTail.size(), " 2 0 0\nSHELLS\n");
  invalid.push_back(wire); // One source wire cannot be assigned twice.
  auto       shell  = valid;
  const auto faceAt = shell.find("FACES\n");
  ASSERT_NE(faceAt, std::string::npos);
  shell.replace(faceAt, std::string("FACES\n0 -1 0").size(), "FACES\n0 -1 1");
  invalid.push_back(shell);
  auto       domain   = valid;
  const auto domainAt = domain.find("COMPONENTS\n") + std::string("COMPONENTS\n").size();
  domain.replace(domainAt, domain.find(' ', domainAt) - domainAt, "-1000");
  invalid.push_back(domain);
  auto destination = value;
  for (const auto& text : invalid)
  {
    std::istringstream input(text);
    EXPECT_THROW(BRepClass3d_ImplicitSolidIO::Read(input, destination), Standard_DomainError);
    EXPECT_TRUE(destination.Validated);
    EXPECT_TRUE(destination.Shape.IsSame(value.Shape));
  }
  auto unsupported                    = value;
  unsupported.Boundary.FaceSources[0] = 99;
  std::ostringstream unchanged;
  unchanged << "untouched";
  EXPECT_THROW(BRepClass3d_ImplicitSolidIO::Write(unsupported, unchanged), Standard_DomainError);
  EXPECT_EQ(unchanged.str(), "untouched");
  Ops::NativeSolids locked;
  BRep_Builder().MakeCompound(locked.Shape);
  locked.Shape.Locked(true);
  std::istringstream input(valid);
  EXPECT_THROW(BRepClass3d_ImplicitSolidIO::Read(input, locked), Standard_DomainError);
  EXPECT_TRUE(locked.Shape.Locked());
  auto loose = value.Boundary.Source.Faces[0];
  loose.Free(true); // Noncanonical ownership must reject, not mutate it.
  std::ostringstream frozen;
  EXPECT_THROW(BRepClass3d_ImplicitSolidIO::Write(value, frozen), Standard_DomainError);
  EXPECT_TRUE(loose.Free());
  EXPECT_TRUE(frozen.str().empty());
  loose.Free(false);
}

TEST(BRepClass3d_ImplicitBox, NativeSolidTopologyIOFullEmptyAndLocale)
{
  struct DecimalComma : std::numpunct<char>
  {
    char do_decimal_point() const override { return ','; }
  };

  using Ops = BRepClass3d_ImplicitBox;
  for (auto operation : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
  {
    const auto value = Ops::BuildNativeSolids(Ops::Build(field(2), box(), operation), .1, 4);
    ASSERT_TRUE(value.Validated) << value.Reason;
    std::ostringstream output;
    output.imbue(std::locale(std::locale::classic(), new DecimalComma));
    const auto precision = output.precision(3);
    BRepClass3d_ImplicitSolidIO::Write(value, output);
    EXPECT_EQ(output.precision(), 3);
    output.precision(precision);
    Ops::NativeSolids  restored;
    std::istringstream input(output.str());
    input.imbue(std::locale(std::locale::classic(), new DecimalComma));
    BRepClass3d_ImplicitSolidIO::Read(input, restored);
    EXPECT_TRUE(restored.Constructed);
    EXPECT_FALSE(restored.Validated);
    EXPECT_EQ(restored.Solids.size(), operation == Ops::Operation::Common ? 1u : 0u);
    for (const auto& face : restored.Boundary.Faces)
      EXPECT_FALSE(BRep_Tool::IsImplicit(face));
  }
}

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

TEST(BRepClass3d_ImplicitBox, StoredTrimGraphCorrespondenceBothOperationsAndRoundtrip)
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
      const auto report = Ops::CheckCorrespondence(source);
      std::cout << "correspondence " << crossing << " " << int(op)
                << " resolved=" << report.Resolved << " reason=" << report.Reason
                << " matched=" << report.MatchedEdges << " bound=" << report.CurveDeviationBound
                << " alignment=" << report.MinimumTangentAlignment << std::endl;
      EXPECT_TRUE(report.Resolved);
      EXPECT_TRUE(report.RegularityResolved);
      EXPECT_GT(report.MatchedVertices, 0);
      EXPECT_EQ(report.MatchedCaps, int(source.Faces.size()) - 1);
      EXPECT_GT(report.MatchedWires, 0);
      EXPECT_GT(report.MinimumTangentAlignment, 0);
      EXPECT_LE(report.ReferenceSeamBound + report.CurveDeviationBound, .1);
      Ops::Boundary      loaded;
      std::istringstream input(record.str());
      BRepClass3d_ImplicitBoxIO::Read(input, loaded);
      EXPECT_FALSE(loaded.SeamBoundVerified);
      EXPECT_TRUE(Ops::CheckCorrespondence(loaded).Resolved);
      EXPECT_FALSE(loaded.SeamBoundVerified);
      if (loaded.Faces.size() > 2)
      {
        std::swap(loaded.Faces[0], loaded.Faces[1]);
        std::swap(loaded.Sources[0], loaded.Sources[1]);
        BRep_Builder builder;
        builder.MakeCompound(loaded.Shape);
        for (const auto& face : loaded.Faces)
          builder.Add(loaded.Shape, face);
        EXPECT_TRUE(Ops::CheckCorrespondence(loaded).Resolved);
      }
      std::ostringstream after;
      BRepClass3d_ImplicitBoxIO::Write(source, after);
      EXPECT_EQ(record.str(), after.str());
    }
}

TEST(BRepClass3d_ImplicitBox, CorrespondenceRefusesChangedMathForgedBoundsAndLimits)
{
  using Ops   = BRepClass3d_ImplicitBox;
  auto source = Ops::Build(field(), box(), Ops::Operation::Common);
  EXPECT_FALSE(Ops::CheckCorrespondence(source, 4).Resolved);
  EXPECT_FALSE(Ops::CheckCorrespondence(source, 512, 1).Resolved);
  EXPECT_THROW(Ops::CheckCorrespondence(source, 3), Standard_DomainError);
  EXPECT_THROW(Ops::CheckCorrespondence(source, 512, 0), Standard_DomainError);
  const auto   old  = source.Faces.back();
  const auto   data = occ::down_cast<BRep_TImplicitFace>(old.TShape());
  TopoDS_Face  changed;
  BRep_Builder builder;
  builder.MakeImplicitFace(changed,
                           new Geom_ImplicitGyroid(10, 1.46),
                           data->Minimum(),
                           data->Maximum(),
                           TopLoc_Location(),
                           .1);
  for (TopoDS_Iterator wire(old); wire.More(); wire.Next())
    builder.Add(changed, wire.Value());
  source.Faces.back() = changed;
  builder.MakeCompound(source.Shape);
  for (const auto& face : source.Faces)
    builder.Add(source.Shape, face);
  source.GeometricBound    = 0;
  source.SeamBoundVerified = true;
  const auto report        = Ops::CheckCorrespondence(source);
  EXPECT_TRUE(report.RegularityResolved);
  EXPECT_FALSE(report.Resolved);
  EXPECT_FALSE(report.Reason.empty());
  const auto caps = Ops::CheckCapAttachment(source);
  EXPECT_FALSE(caps.Resolved);
  EXPECT_FALSE(caps.StoredGraphResolved);
  EXPECT_FALSE(caps.RegionReferenceResolved);
  const auto trim = Ops::CheckTrimAttachment(source);
  EXPECT_FALSE(trim.Resolved);
  EXPECT_FALSE(trim.CapRegionsAttached);
}

TEST(BRepClass3d_ImplicitBox, CorrespondenceAcceptsResamplingAndCyclicWireStarts)
{
  using Ops = BRepClass3d_ImplicitBox;
  auto source =
    Ops::Build(field(0),
               BRepPrimAPI_MakeBox(gp_Pnt(-3.83, -4.07, -2.67), gp_Pnt(4.17, 3.93, 3.37)).Shape(),
               Ops::Operation::Common);
  BRep_Builder builder;
  auto         edge = TopoDS::Edge(
    TopExp_Explorer(source.Faces.back(), TopAbs_EDGE).Current().Oriented(TopAbs_FORWARD));
  TopLoc_Location location;
  double          first, last;
  const auto old = occ::down_cast<Geom_BSplineCurve>(BRep_Tool::Curve(edge, location, first, last));
  std::vector<gp_Pnt> points{old->Pole(1), gp_Pnt((old->Pole(1).XYZ() + old->Pole(2).XYZ()) * .5)};
  for (int i = 2; i <= old->NbPoles(); ++i)
    points.push_back(old->Pole(i));
  const int                  count = int(points.size());
  NCollection_Array1<gp_Pnt> poles(1, count);
  NCollection_Array1<double> knots(1, count);
  NCollection_Array1<int>    mults(1, count);
  for (int i = 1; i <= count; ++i)
  {
    poles(i) = points[size_t(i - 1)];
    knots(i) = i - 1;
    mults(i) = i == 1 || i == count ? 2 : 1;
  }
  builder.UpdateEdge(edge, new Geom_BSplineCurve(poles, knots, mults, 1), .1);
  builder.Range(edge, 0, count - 1);
  TopoDS_Vertex a, b;
  TopExp::Vertices(edge, a, b);
  builder.UpdateVertex(TopoDS::Vertex(a.Oriented(TopAbs_FORWARD)),
                       0,
                       edge,
                       BRep_Tool::Tolerance(a));
  builder.UpdateVertex(TopoDS::Vertex(b.Oriented(TopAbs_REVERSED)),
                       count - 1,
                       edge,
                       BRep_Tool::Tolerance(b));
  for (size_t i = 0; i + 1 < source.Faces.size(); ++i)
  {
    bool used = false;
    for (TopExp_Explorer e(source.Faces[i], TopAbs_EDGE); e.More(); e.Next())
      used |= e.Current().IsSame(edge);
    if (!used)
      continue;
    const auto                   plane = BRepAdaptor_Surface(source.Faces[i]).Plane();
    NCollection_Array1<gp_Pnt2d> uv(1, count);
    for (int p = 1; p <= count; ++p)
    {
      double u, v;
      ElSLib::Parameters(plane, points[size_t(p - 1)], u, v);
      uv(p).SetCoord(u, v);
    }
    builder.UpdateEdge(edge, new Geom2d_BSplineCurve(uv, knots, mults, 1), source.Faces[i], .1);
    builder.Range(edge, source.Faces[i], 0, count - 1);
  }
  builder.SameRange(edge, true);
  builder.SameParameter(edge, true);
  const auto report = Ops::CheckCorrespondence(source);
  EXPECT_TRUE(report.Resolved) << report.Reason;

  const auto  original = source.Faces.back();
  const auto  data     = occ::down_cast<BRep_TImplicitFace>(original.TShape());
  TopoDS_Face changed;
  builder.MakeImplicitFace(changed,
                           data->Field(),
                           data->Minimum(),
                           data->Maximum(),
                           TopLoc_Location(),
                           .1);
  for (TopoDS_Iterator wire(original); wire.More(); wire.Next())
  {
    std::vector<TopoDS_Shape> edges;
    for (TopoDS_Iterator e(wire.Value()); e.More(); e.Next())
      edges.push_back(e.Value());
    if (edges.size() > 1)
      std::rotate(edges.begin(), edges.begin() + 1, edges.end());
    TopoDS_Wire rotated;
    builder.MakeWire(rotated);
    for (const auto& e : edges)
      builder.Add(rotated, e);
    builder.Add(changed, rotated);
  }
  source.Faces.back() = changed;
  builder.MakeCompound(source.Shape);
  for (const auto& face : source.Faces)
    builder.Add(source.Shape, face);
  const auto rotated = Ops::CheckCorrespondence(source);
  EXPECT_TRUE(rotated.Resolved) << rotated.Reason;
  const auto attached = Ops::CheckTrimAttachment(source);
  EXPECT_TRUE(attached.Resolved) << attached.Reason;
}

TEST(BRepClass3d_ImplicitBox, CorrespondenceChecksMultipleRegionsAndCapParameterGeometry)
{
  using Ops       = BRepClass3d_ImplicitBox;
  const auto wide = BRepPrimAPI_MakeBox(gp_Pnt(-.31, -.37, 1.25), gp_Pnt(12.93, 2.91, 4)).Shape();
  for (auto operation : {Ops::Operation::Common, Ops::Operation::BoxMinusField})
  {
    const auto source = Ops::Build(field(), wide, operation);
    const auto report = Ops::CheckCorrespondence(source);
    EXPECT_TRUE(report.Resolved) << report.Reason;
    EXPECT_EQ(report.MatchedWires, 2);
  }
  auto       bad  = Ops::Build(field(), box(), Ops::Operation::BoxMinusField);
  const auto cap  = bad.Faces.front();
  const auto edge = TopoDS::Edge(TopExp_Explorer(cap, TopAbs_EDGE).Current());
  double     first, last;
  const auto original =
    occ::down_cast<Geom2d_BSplineCurve>(BRep_Tool::CurveOnSurface(edge, cap, first, last));
  auto changed = occ::down_cast<Geom2d_BSplineCurve>(original->Copy());
  auto point   = changed->Pole(2);
  point.SetX(point.X() + 1.e-6);
  changed->SetPole(2, point);
  BRep_Builder().UpdateEdge(edge, changed, cap, .1);
  // This small interior p-curve deviation is within conventional edge tolerance,
  // but must NOT inherit the tighter fresh-trace correspondence certificate.
  ASSERT_TRUE(BRepCheck_Analyzer(cap).IsValid());
  const auto report = Ops::CheckCorrespondence(bad);
  EXPECT_FALSE(report.Resolved);
  EXPECT_EQ(report.Reason, "cap_curve_parameter_mismatch");
}
