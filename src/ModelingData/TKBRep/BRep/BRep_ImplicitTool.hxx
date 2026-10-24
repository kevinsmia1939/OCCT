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
#include <string>

//! Field-aware queries on experimental native implicit faces.
//! Points/directions/distances are in the placed face frame, not field-local UV.
//! Honors explicit native plane-patch predicates and closed domain-exclusion
//! boxes in domain/ray queries. Value/Gradient still query the underlying field;
//! IsInDomain alone is not an on-surface or distance-tolerant membership test.
//! Does NOT implement general trimming, BOP, generic meshing or validation.
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

  struct Projection
  {
    bool Converged = false;
    gp_Pnt Point, First, Last;
    double Residual = 0, Stationarity = 0, Distance = 0, BracketWidth = 0;
    int Iterations = 0;
    std::string Reason = "iteration_limit";
  };

  //! Bounded local regular normal foot, NOT a globally nearest point search.
  //! Numerical KKT solve with a positive tangent second variation; acceptance
  //! needs a nearby sign bracket wholly in the explicit query domain. Honors
  //! domain boxes/exclusions/patch predicates, NOT arbitrary added trim wires.
  //! Ordinary floating reserves, not directed interval certification. Failure
  //! returns Converged=false; invalid settings/nonfinite evaluation throw.
  Standard_EXPORT static Projection Project(const TopoDS_Face&, const gp_Pnt& theSeed,
                                              double theDistanceTolerance = 0.1,
                                              double theMaximumDistance = 1.0,
                                              int theMaximumIterations = 32);

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
