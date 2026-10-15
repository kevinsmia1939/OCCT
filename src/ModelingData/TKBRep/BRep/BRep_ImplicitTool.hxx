// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#ifndef _BRep_ImplicitTool_HeaderFile
#define _BRep_ImplicitTool_HeaderFile
#include <TopoDS_Face.hxx>
#include <Bnd_Box.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <vector>

//! Field-aware queries on experimental native implicit faces.
//! Points/directions/distances are in the placed face frame, not field-local UV.
//! Does NOT implement trimming, BOP, generic meshing, validation or exchange.
class BRep_ImplicitTool
{
public:
  struct RayHit
  {
    double Distance;
    gp_Pnt Point;
    double Residual;
    double BracketWidth;
  };

  Standard_EXPORT static double  Value(const TopoDS_Face&, const gp_Pnt&);
  Standard_EXPORT static gp_Vec  Gradient(const TopoDS_Face&, const gp_Pnt&);
  Standard_EXPORT static gp_Vec  Normal(const TopoDS_Face&, const gp_Pnt&);
  Standard_EXPORT static Bnd_Box Bounds(const TopoDS_Face&);
  Standard_EXPORT static bool    IsInDomain(const TopoDS_Face&, const gp_Pnt&);

  //! Sampled sign brackets, NOT root completeness/tangency certification.
  //! Direction is normalized. Zero distance tolerance uses the face tolerance.
  //! Refinement/merge distances have model-length units; residual has field units.
  Standard_EXPORT static std::vector<RayHit> IntersectRay(const TopoDS_Face&,
                                                          const gp_Pnt& theOrigin,
                                                          const gp_Vec& theDirection,
                                                          double        theMinimumDistance,
                                                          double        theMaximumDistance,
                                                          int           theSamples           = 512,
                                                          double        theDistanceTolerance = 0.0);
};
#endif
