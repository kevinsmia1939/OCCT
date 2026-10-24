// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRep_ImplicitCylinderSection_HeaderFile
#define _BRep_ImplicitCylinderSection_HeaderFile
#include <BRep_ImplicitPlaneSection.hxx>

//! Derivative-bounded gyroid section on a finite open cylindrical chart.
//! This query does not create trim topology, join periodic seams or build solids.
class BRep_ImplicitCylinderSection
{
public:
  struct Settings
  {
    gp_Ax2 Frame; //!< Axis is Z; angle zero is X.
    double Radius = 1;
    double AngleMinimum = 0, AngleMaximum = 1; //!< Radians, span strictly below 2*pi.
    double HeightMinimum = -1, HeightMaximum = 1;
    double DistanceTolerance = .1;
    int MinimumCells = 4, MaximumCells = 512; //!< Powers of two, 4..512.
    void (*WorkCheckpoint)() = nullptr;
  };
  using Result = BRep_ImplicitPlaneSection::Result;
  //! Exact gyroid on a childless unrestricted native face. The whole conservative
  //! cylinder enclosure must be inside its computational domain. Zero samples,
  //! tangencies and unproved cells return unresolved, never guessed topology.
  //! GeometricBound is a model-length chart-cell diameter bound, not residual.
  Standard_EXPORT static Result Trace(const TopoDS_Face&, const Settings&);
};
#endif
