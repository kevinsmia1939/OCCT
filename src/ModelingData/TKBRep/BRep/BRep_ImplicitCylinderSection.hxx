// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRep_ImplicitCylinderSection_HeaderFile
#define _BRep_ImplicitCylinderSection_HeaderFile
#include <BRep_ImplicitPlaneSection.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Pnt2d.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <array>

//! Derivative-bounded gyroid section on a finite open cylindrical chart.
//! Curves remain transient approximations, not native trim topology or solids.
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
  //! Full-turn query with one canonical shared chart seam (no distance welding).
  //! AngleMaximum must equal AngleMinimum + 2*pi. Seam roots are indexed once;
  //! only height boundaries may contain open ends. Existing Settings/Result
  //! layouts and the open Trace contract are unchanged.
  Standard_EXPORT static Result TracePeriodic(const TopoDS_Face&, const Settings&);
  struct BoundaryRoot
  {
    TopoDS_Vertex Vertex; //!< Owned canonical vertex for reuse by side/end-cap trims.
    gp_Pnt First, Last; //!< Root bracket on the exact circular height boundary.
    double AngleFirst = 0, AngleLast = 0;
    int HeightBoundary = 0; //!< 0 bottom, 1 top.
  };
  struct BoundaryResult
  {
    Result Section;
    std::vector<BoundaryRoot> Roots; //!< Ordered by height boundary then chart angle.
    //! One pair per section curve; 1-based Roots indices. Closed curves use 0,0.
    std::vector<std::array<size_t, 2>> CurveRoots;
    bool BoundaryIncidenceVerified = false;
  };
  //! Retain the exact root/endpoint identity from one periodic tracing pass.
  //! Does not build end-cap regions, native trim edges, shells or solids.
  //! Unresolved results publish no roots or endpoint identities.
  Standard_EXPORT static BoundaryResult TracePeriodicBoundary(const TopoDS_Face&, const Settings&);
  struct ChartVertex
  {
    gp_Pnt Point;
    gp_Pnt2d UV; //!< Canonical angle (radians) and axial height.
    bool Seam = false; //!< Logical shared first/last chart-column root.
    TopoDS_Vertex Vertex; //!< Non-null only for cap or chart-seam roots.
  };
  struct ChartArc
  {
    size_t First = 0, Last = 0;
    int U = 0, V = 0; //!< Accepted source-cell indices, not inferred from points.
  };
  struct PeriodicGraph
  {
    BoundaryResult Boundary;
    std::vector<ChartVertex> Vertices;
    std::vector<ChartArc> Arcs;
    //! Node IDs in exactly the order of each Section curve's Points.
    std::vector<std::vector<size_t>> CurveVertices;
    //! Continuous angle/height coordinates; seam branches follow source cells.
    std::vector<std::vector<gp_Pnt2d>> CurveUV;
    bool Verified = false; //!< Provenance only; not native side-face validity.
  };
  //! Retains periodic chart provenance without proximity welding. No graph is
  //! published after an unresolved pass. Does not construct side faces/solids.
  Standard_EXPORT static PeriodicGraph TracePeriodicGraph(const TopoDS_Face&, const Settings&);
  struct RingArc
  {
    TopoDS_Edge Edge; //!< Analytic circle, outward cap orientation (bottom reversed).
    size_t LowerRoot = 0, UpperRoot = 0; //!< Root IDs; 0 denotes synthetic root-free/chart-seam ends.
    bool Negative = false; //!< F<0 between root brackets; endpoint uncertainty is retained.
  };
  struct RingResult
  {
    BoundaryResult Boundary;
    std::array<std::vector<RingArc>, 2> Rings;
    bool Verified = false; //!< Ring incidence/sign intervals only; NOT cap/solid validation.
  };
  //! Construct both height-circle partitions using the retained root vertices.
  //! Returns no edges unless every interval and root incidence passes.
  //! Does not construct cap regions, cylindrical faces, shells or solids.
  Standard_EXPORT static RingResult BuildBoundaryRings(const TopoDS_Face&, const Settings&);
  struct ChartSkeleton
  {
    PeriodicGraph Graph;
    RingResult Rings;
    std::array<TopoDS_Vertex, 2> SeamEnds;
    occ::handle<Geom_CylindricalSurface> Surface;
    std::vector<TopoDS_Edge> SeamEdges; //!< Increasing height, paired angular p-curves.
    bool Verified = false;
  };
  //! One shared side/ring snapshot. Wrap arcs split at the chart seam, using
  //! a single vertex per cap for both parameter branches. No side faces yet.
  Standard_EXPORT static ChartSkeleton BuildChartSkeleton(const TopoDS_Face&, const Settings&);
  struct DiskResult
  {
    BoundaryResult Boundary;
    std::array<Result, 2> Disks;
    std::array<std::vector<std::array<size_t, 2>>, 2> CurveRoots;
    bool Verified = false; //!< Constrained contour/root incidence, NOT native cap/solid validity.
  };
  //! Trace both planar disks using accepted Cartesian cells and shared circle
  //! root brackets. Includes interior closed contours; never proximity-welds.
  //! Uncertain cut-cell attachment rejects; no partial disk curves on failure.
  Standard_EXPORT static DiskResult TraceEndCapDisks(const TopoDS_Face&, const Settings&);
  struct CapSkeleton
  {
    DiskResult Disks;
    RingResult Rings;
    bool Verified = false;
  };
  //! One tracing snapshot: rings and disk endpoints share actual vertex handles.
  Standard_EXPORT static CapSkeleton BuildCapSkeleton(const TopoDS_Face&, const Settings&);
  struct ChartCapSkeleton
  {
    ChartSkeleton Chart;
    CapSkeleton Caps;
    bool Verified = false;
  };
  //! Side chart, seamed rings and both cap disks from one boundary snapshot.
  Standard_EXPORT static ChartCapSkeleton BuildChartCapSkeleton(const TopoDS_Face&, const Settings&);
private:
  static DiskResult traceEndCapDisks(const TopoDS_Face&, const Settings&, const BoundaryResult&);
};
#endif
