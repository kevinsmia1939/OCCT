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

namespace
{
TopoDS_Face field(double iso = 1.45)
{
  TopoDS_Face face;
  BRep_Builder().MakeImplicitFace(face,
                                  new Geom_ImplicitGyroid(10, iso),
                                  gp_Pnt(-5, -5, -5),
                                  gp_Pnt(6, 6, 6),
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
  EXPECT_THROW(Ops::Build(field(0), box(), Ops::Operation::Common), Standard_NotImplemented);
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
