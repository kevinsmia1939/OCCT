// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepClass3d_SignedDistance.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRep_Builder.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_Failure.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Compound.hxx>
#include <cmath>

BRepClass3d_SignedDistance::BRepClass3d_SignedDistance(const TopoDS_Shape& theSolid,
                                                       const double        theTolerance)
    : myTolerance(theTolerance)
{
  if (!std::isfinite(theTolerance) || theTolerance <= 0.0)
  {
    throw Standard_DomainError("Signed distance requires a finite positive tolerance");
  }
  if (theSolid.IsNull() || theSolid.ShapeType() != TopAbs_SOLID
      || !BRepCheck_Analyzer(theSolid).IsValid())
  {
    throw Standard_DomainError("Signed distance requires one valid B-Rep solid");
  }
  myClassifier.Load(theSolid);
  myClassifier.PerformInfinitePoint(myTolerance);
  if (myClassifier.State() != TopAbs_OUT)
  {
    throw Standard_DomainError("Signed distance requires a finite outward-oriented solid");
  }

  // DistShapeShape(solid, interior vertex) reports an inner solution at zero.
  // Measure against the boundary faces instead, preserving cavity boundaries.
  BRep_Builder    aBuilder;
  TopoDS_Compound aBoundary;
  aBuilder.MakeCompound(aBoundary);
  bool hasFaces = false;
  for (TopExp_Explorer anExplorer(theSolid, TopAbs_FACE); anExplorer.More(); anExplorer.Next())
  {
    aBuilder.Add(aBoundary, anExplorer.Current());
    hasFaces = true;
  }
  if (!hasFaces)
  {
    throw Standard_DomainError("Signed distance solid has no boundary faces");
  }
  myBoundary = aBoundary;
}

double BRepClass3d_SignedDistance::Value(const gp_Pnt& thePoint)
{
  if (!std::isfinite(thePoint.X()) || !std::isfinite(thePoint.Y()) || !std::isfinite(thePoint.Z()))
  {
    throw Standard_DomainError("Signed distance point must be finite");
  }
  BRepExtrema_DistShapeShape aDistance(myBoundary,
                                       BRepBuilderAPI_MakeVertex(thePoint).Shape(),
                                       myTolerance);
  if (!aDistance.IsDone() || aDistance.NbSolution() == 0 || !std::isfinite(aDistance.Value()))
  {
    throw Standard_Failure("Cannot resolve B-Rep boundary distance");
  }
  const double aValue = aDistance.Value();
  myClassifier.Perform(thePoint, myTolerance);
  const auto aState = myClassifier.State();
  if (aState == TopAbs_ON || aValue <= myTolerance)
  {
    return 0.0;
  }
  if (aState == TopAbs_IN)
  {
    return -aValue;
  }
  if (aState == TopAbs_OUT)
  {
    return aValue;
  }
  throw Standard_Failure("Cannot classify B-Rep signed distance point");
}
