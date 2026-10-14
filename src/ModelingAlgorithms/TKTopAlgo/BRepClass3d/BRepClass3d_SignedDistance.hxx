// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#ifndef _BRepClass3d_SignedDistance_HeaderFile
#define _BRepClass3d_SignedDistance_HeaderFile

#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Pnt.hxx>

//! Numerical signed distance to the boundary of one finite, valid solid.
//! Negative values denote material, positive values the exterior or a cavity.
//! This evaluator bridges a conventional B-Rep to field algorithms without
//! changing its topology or claiming that arbitrary implicit faces are B-Reps.
//! The input shape must not be mutated during the lifetime of the evaluator.
//! Instances retain classification/distance state and must not be shared across threads.
class BRepClass3d_SignedDistance
{
public:
  //! Rejects null, invalid, non-solid and unbounded/reversed solids and invalid tolerance.
  Standard_EXPORT BRepClass3d_SignedDistance(const TopoDS_Shape& theSolid,
                                             const double        theTolerance = 1.0e-7);
  BRepClass3d_SignedDistance(const BRepClass3d_SignedDistance&)            = delete;
  BRepClass3d_SignedDistance& operator=(const BRepClass3d_SignedDistance&) = delete;

  //! Computes nearest boundary distance and uses the solid classifier for sign.
  //! Points classified ON or within tolerance of the boundary return zero.
  //! Raises Standard_Failure if distance or classification cannot be resolved.
  Standard_EXPORT double Value(const gp_Pnt& thePoint);

private:
  TopoDS_Shape                myBoundary;
  BRepClass3d_SolidClassifier myClassifier;
  BRepExtrema_DistShapeShape  myDistance;
  double                      myTolerance;
};

#endif
