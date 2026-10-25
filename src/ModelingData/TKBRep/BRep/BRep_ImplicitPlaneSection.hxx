// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRep_ImplicitPlaneSection_HeaderFile
#define _BRep_ImplicitPlaneSection_HeaderFile
#include <TopoDS_Face.hxx>
#include <gp_Ax2.hxx>
#include <gp_Pnt.hxx>
#include <vector>

//! Native gyroid/plane-rectangle section with derivative-bound cell acceptance.
//! Polylines are tolerance-bearing approximations, not exact Geom_Curves or
//! trimmed faces. Unresolved cells prohibit downstream trim construction.
class BRep_ImplicitPlaneSection
{
public:
  struct Settings
  {
    gp_Ax2 Frame; //!< World-coordinate orthonormal plane frame.
    double UMinimum = -1, UMaximum = 1, VMinimum = -1, VMaximum = 1;
    double DistanceTolerance = 0.1;              //!< Model length, NOT gyroid field units.
    int    MinimumCells = 4, MaximumCells = 512; //!< Powers of two, 4..512.
    //! Explicit topology policy: recognize exact lattice nodes by identity,
    //! with coordinate-monotone cell cones. Identity field placement and
    //! coordinate-plane frames only; rounded affine identities do not count.
    //! Consumers MUST handle PointContacts separately from one-dimensional curves.
    bool AllowExactLatticeNodes = false;
    //! Optional synchronous caller-thread checkpoint; may throw to cancel.
    //! Never retained in results. Avoids a dependency on higher-level work controls.
    void (*WorkCheckpoint)() = nullptr;
  };

  struct Curve
  {
    std::vector<gp_Pnt> Points; //!< A closed curve repeats its first point.
    bool                Closed = false;
  };

  struct Result
  {
    bool               Resolved   = false; //!< ALL cells passed empty/regular-arc/isolated-corner bounds.
    int                Resolution = 0; //!< May be a non-dyadic fallback grid, within MaximumCells.
    int                UnresolvedCells = 0, PotentialCriticalCells = 0, AmbiguousCells = 0;
    int                GraphDefects    = 0; //!< Branching/interior-ended fragment vertices.
    double             GeometricBound  = 0; //!< Max accepted cell diagonal; partial if unresolved.
    double             MaximumResidual = 0; //!< Field units, diagnostic only.
    double             MaximumRootBracket = 0; //!< Model length.
    int                ExactLatticeNodeCount = 0;
    std::vector<gp_Pnt> PointContacts; //!< Certified isolated rectangle-corner zeros, not curves.
    std::vector<Curve> Curves;
  };

  //! Supports exact Geom_ImplicitGyroid fields on childless BRep_TImplicitFace.
  //! Rectangle must lie entirely inside the computational domain. Placements
  //! are supported. Errors reject invalid/unsupported inputs; valid but
  //! unresolved geometry returns Resolved=false, never an empty-success claim.
  Standard_EXPORT static Result Trace(const TopoDS_Face&, const Settings&);
  struct CellArc
  {
    int U = 0, V = 0;
    size_t First = 0, Last = 0; //!< Shared Vertices indices, zero-based.
  };
  struct GraphResult
  {
    Result Section;
    std::vector<gp_Pnt> Vertices;
    std::vector<CellArc> Arcs; //!< Accepted regular-cell arcs, not display-mesh segments.
    double VertexRadius = 0; //!< Conservative global root-bracket halfwidth + coordinate reserve.
  };
  //! Retain accepted cell/root identities for constrained-domain tracing.
  //! Unresolved queries publish no graph; no native trim/solid qualification.
  Standard_EXPORT static GraphResult TraceGraph(const TopoDS_Face&, const Settings&);
};
#endif
