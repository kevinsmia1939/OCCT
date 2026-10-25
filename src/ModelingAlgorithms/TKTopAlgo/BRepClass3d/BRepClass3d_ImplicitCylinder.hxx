// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRepClass3d_ImplicitCylinder_HeaderFile
#define _BRepClass3d_ImplicitCylinder_HeaderFile
#include <BRep_ImplicitCylinderSection.hxx>
#include <string>
#include <TopoDS_Wire.hxx>

//! Restricted native cap-region assembly; NOT a cylinder Boolean/solid builder.
class BRepClass3d_ImplicitCylinder
{
public:
  struct Caps
  {
    TopoDS_Face Field; //!< Authoritative equation/domain snapshot, not replaced by fitted edges.
    BRep_ImplicitCylinderSection::Settings Parameters; //!< Work checkpoint never retained.
    bool Negative = true;
    size_t MaximumNestingTests = 2000000;
    BRep_ImplicitCylinderSection::CapSkeleton Source;
    std::array<std::vector<TopoDS_Face>, 2> Faces;
    std::array<std::vector<TopoDS_Edge>, 2> FieldEdges;
    bool TopologyValidated = false;
    double ContourBound = 0;
    std::string Reason;
  };
  //! Negative selects F<=0; positive selects F>=0. Planar faces use analytic
  //! circles and bounded degree-1 B-spline field seams, not mesh triangles.
  //! Failure returns no faces/field edges. Source field geometry stays immutable.
  //! Cap topology/contour bounds do not certify cylindrical side faces or solids.
  Standard_EXPORT static Caps BuildCaps(const TopoDS_Face&,
    const BRep_ImplicitCylinderSection::Settings&, bool negative = true,
    size_t maximumNestingTests = 2000000);
  struct ChartCaps
  {
    BRep_ImplicitCylinderSection::ChartSkeleton Chart;
    Caps Regions;
  };
  //! Actual planar cap regions reuse the seamed side-chart ring edges and roots.
  //! Does not construct cylinder-side faces, shells or solids.
  Standard_EXPORT static ChartCaps BuildChartCaps(const TopoDS_Face&,
    const BRep_ImplicitCylinderSection::Settings&, bool negative = true,
    size_t maximumNestingTests = 2000000);
  struct SideBoundary
  {
    ChartCaps Source;
    std::vector<TopoDS_Edge> FieldEdges;
    bool Verified = false;
    double SurfaceDeviationBound = 0;
    std::string Reason;
  };
  //! Split field contours only at logical seam/cap vertices. Degree-1 3D and
  //! UV curves retain a curvature-based cylinder deviation allowance. Does not
  //! yet certify material cycles, side faces, implicit attachments or solids.
  Standard_EXPORT static SideBoundary BuildSideBoundary(const TopoDS_Face&,
    const BRep_ImplicitCylinderSection::Settings&, bool negative = true,
    size_t maximumNestingTests = 2000000);
  struct Sides
  {
    SideBoundary Boundary;
    std::vector<TopoDS_Face> Faces;
    bool TopologyValidated = false;
    std::string Reason;
  };
  //! Retained cylindrical material regions. No implicit-face attachment or
  //! closed-solid qualification follows from ordinary side-face validity.
  Standard_EXPORT static Sides BuildSides(const TopoDS_Face&,
    const BRep_ImplicitCylinderSection::Settings&, bool negative = true,
    size_t maximumNestingTests = 2000000);
  struct TrimBoundary
  {
    Sides Source;
    std::vector<TopoDS_Wire> Wires; //!< Opposite conventional-face incidence.
    std::vector<int> EdgeWires; //!< Side edges, bottom edges, then top edges.
    bool IncidenceVerified = false;
    std::string Reason;
  };
  //! Qualify shared conventional-edge incidence and construct oriented field
  //! boundary wires. Does not assign wires to implicit surface components or
  //! attach implicit faces; no shell/solid proof follows from this graph.
  Standard_EXPORT static TrimBoundary BuildTrimBoundary(const TopoDS_Face&,
    const BRep_ImplicitCylinderSection::Settings&, bool negative = true,
    size_t maximumNestingTests = 2000000);
};
#endif
