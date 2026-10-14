// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_SignedDistance.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <Standard_Failure.hxx>
#include <gtest/gtest.h>
#include <limits>

TEST(BRepClass3d_SignedDistanceTest, InteriorIsNotZero)
{
  const auto                 aBox = BRepPrimAPI_MakeBox(2.0, 2.0, 2.0).Shape();
  BRepClass3d_SignedDistance aField(aBox);
  EXPECT_NEAR(aField.Value(gp_Pnt(1, 1, 1)), -1.0, 1.0e-7);
  EXPECT_NEAR(aField.Value(gp_Pnt(3, 1, 1)), 1.0, 1.0e-7);
  EXPECT_DOUBLE_EQ(aField.Value(gp_Pnt(2, 1, 1)), 0.0);
}

TEST(BRepClass3d_SignedDistanceTest, RejectsInvalidInput)
{
  EXPECT_THROW((void)BRepClass3d_SignedDistance(TopoDS_Shape()), Standard_Failure);
  const auto aBox = BRepPrimAPI_MakeBox(2.0, 2.0, 2.0).Shape();
  EXPECT_THROW(BRepClass3d_SignedDistance(aBox, 0.0), Standard_Failure);
  EXPECT_THROW(BRepClass3d_SignedDistance(aBox.Reversed()), Standard_Failure);
}

TEST(BRepClass3d_SignedDistanceTest, RepeatedQueriesResetMovingOperand)
{
  const auto                 aBox = BRepPrimAPI_MakeBox(2.0, 2.0, 2.0).Shape();
  BRepClass3d_SignedDistance aField(aBox);
  for (int anIndex = 0; anIndex < 20; ++anIndex)
  {
    EXPECT_NEAR(aField.Value(gp_Pnt(1, 1, 1)), -1.0, 1.0e-7);
    EXPECT_NEAR(aField.Value(gp_Pnt(5, 1, 1)), 3.0, 1.0e-7);
    EXPECT_DOUBLE_EQ(aField.Value(gp_Pnt(0, 1, 1)), 0.0);
  }
  EXPECT_THROW(aField.Value(gp_Pnt(std::numeric_limits<double>::infinity(), 0, 0)),
               Standard_Failure);
  EXPECT_NEAR(aField.Value(gp_Pnt(1, 1, 1)), -1.0, 1.0e-7);
}
