// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRepClass3d_ImplicitCylinder_HeaderFile
#define _BRepClass3d_ImplicitCylinder_HeaderFile
#include <BRep_ImplicitCylinderSection.hxx>
#include <string>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopAbs_State.hxx>

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
  struct AxialSheet
  {
    TrimBoundary Boundary;
    bool MonotonicityVerified = false;
    int Direction = 0; //!< Sign of the field derivative along the cylinder axis.
    double DerivativeLowerBound = 0; //!< Positive bound on its signed magnitude.
    size_t Cells = 0;
    std::string Reason;
  };
  //! Adaptive Hessian cones over the entire finite cylinder prove at most one
  //! zero per axial fibre. This does not establish existence, connectedness,
  //! projected trim correspondence, implicit-face attachment or solid validity.
  Standard_EXPORT static AxialSheet CheckAxialSheet(const TopoDS_Face&,
    const BRep_ImplicitCylinderSection::Settings&, bool negative = true,
    size_t maximumNestingTests = 2000000, int maximumDepth = 24,
    size_t maximumCells = 200000);
  struct ProjectedRegions
  {
    AxialSheet Source;
    std::vector<TopoDS_Face> Faces; //!< Conventional planar reference regions.
    std::vector<int> WireRegions; //!< Native boundary wire -> reference region.
    bool ProxyTopologyValidated = false;
    bool MathematicalCorrespondenceVerified = false;
    size_t CorrespondenceCells = 0;
    std::string Reason;
  };
  //! Project cap seams affinely and side seams as analytic circular arcs after
  //! checking monotone angular correspondence. Qualifies the projected native
  //! proxy only; mathematical correspondence requires a separate field gate.
  Standard_EXPORT static ProjectedRegions BuildProjectedRegions(const TopoDS_Face&,
    const BRep_ImplicitCylinderSection::Settings&, bool negative = true,
    size_t maximumNestingTests = 2000000, int maximumDepth = 24,
    size_t maximumCells = 200000);
  //! Sufficient full-disk graph correspondence: proves opposite strict signs
  //! over both end disks and one complete circular side boundary. Partial and
  //! multi-sheet regions remain unqualified; failure retains proxy diagnostics.
  Standard_EXPORT static ProjectedRegions CheckProjectedCorrespondence(const TopoDS_Face&,
    const BRep_ImplicitCylinderSection::Settings&, bool negative = true,
    size_t maximumNestingTests = 2000000, int maximumDepth = 24,
    size_t maximumCells = 200000);
  struct NativeFaces
  {
    ProjectedRegions Reference;
    std::vector<TopoDS_Face> Faces;
    bool FacesAttached = false;
    std::string Reason;
  };
  //! Actual cylinder-restricted implicit face for qualified full-disk sheets.
  //! Shared trim incidence is checked; shells/solids and general partial-sheet
  //! or multi-sheet component attachment remain separate requirements.
  Standard_EXPORT static NativeFaces BuildImplicitFaces(const TopoDS_Face&,
    const BRep_ImplicitCylinderSection::Settings&, bool negative = true,
    size_t maximumNestingTests = 2000000, int maximumDepth = 24,
    size_t maximumCells = 200000);
  struct NativeShell
  {
    NativeFaces Source;
    TopoDS_Shell Shell;
    std::vector<TopoDS_Face> Faces;
    bool TopologyValidated = false;
    int Euler = 0;
    size_t Vertices = 0, Edges = 0;
    std::string Reason;
  };
  //! Closed oriented shell for qualified full-disk sheets. Includes periodic
  //! seam occurrences and vertex-link checks; does not qualify solid geometry,
  //! volume, persistence or ordinary OCCT consumers of implicit faces.
  Standard_EXPORT static NativeShell BuildNativeShell(const TopoDS_Face&,
    const BRep_ImplicitCylinderSection::Settings&, bool negative = true,
    size_t maximumNestingTests = 2000000, int maximumDepth = 24,
    size_t maximumCells = 200000);
  struct ShellTolerance
  {
    NativeShell Source;
    bool ToleranceVerified = false;
    double ContourBound = 0, SurfaceDeviationBound = 0, FieldWireBound = 0;
    double MaximumFaceTolerance = 0, MaximumEdgeTolerance = 0, MaximumVertexTolerance = 0;
    std::string Reason;
  };
  //! Audit declared model-length tolerances and combined field-wire bounds.
  //! This is not a global shell-geometry, material-volume or solid certificate.
  Standard_EXPORT static ShellTolerance CheckShellTolerance(const TopoDS_Face&,
    const BRep_ImplicitCylinderSection::Settings&, bool negative = true,
    size_t maximumNestingTests = 2000000, int maximumDepth = 24,
    size_t maximumCells = 200000);
  struct ShellGeometry
  {
    ShellTolerance Source;
    bool BoundaryVerified = false;
    double MathematicalToNativeBound = 0, NativeToMathematicalBound = 0, BoundaryBound = 0;
    size_t Cells = 0, SegmentTests = 0;
    std::string Reason;
  };
  //! Bidirectional boundary-set bounds for the qualified full-disk material.
  //! Covers the entire cylindrical chart using field cones and native trim
  //! witnesses. Exact field/cap correspondence and stored wire bounds are
  //! inherited. Does not construct a solid or qualify volume/general consumers.
  Standard_EXPORT static ShellGeometry CheckShellGeometry(const TopoDS_Face&,
    const BRep_ImplicitCylinderSection::Settings&, bool negative = true,
    size_t maximumNestingTests = 2000000, int maximumDepth = 24,
    size_t maximumCells = 200000, int maximumGeometryDepth = 24,
    size_t maximumGeometryCells = 200000, size_t maximumSegmentTests = 20000000);
  struct MaterialPoint
  {
    TopAbs_State State = TopAbs_UNKNOWN; //!< Mathematical material, with query tolerance.
    double DistanceLowerBound = 0, DistanceUpperBound = 0;
    double FieldValueAtDomainProjection = 0;
    bool BoundaryAllowanceSeparated = false;
  };
  //! Classify the unchanged owned full-disk recipe in model-length tolerance.
  //! Uses spatial distance bounds, not a dimensionless |F| tolerance. UNKNOWN
  //! remains explicit when bounds do not separate the requested distance band.
  //! No ordinary OCCT solid-classifier or generic consumer support is implied.
  Standard_EXPORT static MaterialPoint ClassifyMaterial(const ShellGeometry&,
    const gp_Pnt&, double pointTolerance = 1.e-7);
  struct VolumeSettings
  {
    double AbsoluteTolerance = .01; //!< Model length cubed, independent of trim tolerance.
    int MaximumDepth = 24;
    size_t MaximumCells = 200000, MaximumRootTests = 10000000;
    void (*WorkCheckpoint)() = nullptr; //!< Never retained in results.
  };
  struct VolumeBounds
  {
    bool Verified = false;
    double Lower = 0, Upper = 0, CylinderVolume = 0, AbsoluteTolerance = 0;
    size_t Cells = 0, RootTests = 0;
    std::string Reason;
  };
  //! Adaptive disk/fibre integration of the mathematical material. Conservative
  //! area and Taylor-cone root intervals; no mesh or generic BRep mass property.
  Standard_EXPORT static VolumeBounds BoundMaterialVolume(const ShellGeometry&, const VolumeSettings&);
  struct NativeSolid
  {
    ShellGeometry Source;
    VolumeBounds Volume;
    TopoDS_Solid Solid;
    bool SolidValidated = false;
    std::string Reason;
  };
  //! Actual native solid only after topology, boundary geometry and bounded
  //! positive material volume qualify. Validity is the explicit implicit-aware
  //! contract, not generic OCCT consumer/serialization support.
  Standard_EXPORT static NativeSolid BuildNativeSolid(const ShellGeometry&, const VolumeSettings&);
  //! Structural shared-topology validation only. Does not restore certificates.
  Standard_EXPORT static void ValidateStoredShellTopology(const TopoDS_Shell&,
    const std::vector<TopoDS_Face>&);
};
#endif
