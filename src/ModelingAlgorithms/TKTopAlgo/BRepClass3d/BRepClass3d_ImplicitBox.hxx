// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRepClass3d_ImplicitBox_HeaderFile
#define _BRepClass3d_ImplicitBox_HeaderFile
#include <BRep_ImplicitPlaneSection.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <vector>
#include <array>
#include <string>

//! Experimental hybrid boundary constructor, NOT generic BOP or a solid validator.
class BRepClass3d_ImplicitBox
{
public:
  enum class Operation
  {
    Common,
    BoxMinusField
  };

  enum class CornerPolicy
  {
    Reject,
    Regularized //!< Closure of 3-D interior; isolated closed-inequality contacts recorded separately.
  };

  struct Boundary
  {
    TopoDS_Compound          Shape;
    std::vector<TopoDS_Face> Faces;
    std::vector<int>         Sources; //!< -1 = implicit, otherwise input box face traversal index.
    std::vector<BRep_ImplicitPlaneSection::Result> Sections;
    int                                            SeamCount = 0;
    int    BoxEdgeRootCount  = 0;     //!< Shared roots, each used by two adjacent face contours.
    double GeometricBound    = 0;     //!< Seam bound ONLY, not a surface/solid certificate.
    bool   SeamBoundVerified = false; //!< True only after tracing, false after topology I/O.
    CornerPolicy Policy = CornerPolicy::Reject;
    int ExactCornerCount = 0, CornerContourEnds = 0, PlanarCornerContacts = 0;
    std::vector<gp_Pnt> ExcludedCornerContacts; //!< Not vertices of regularized boundary topology.
  };

  enum class PointState
  {
    Inside,
    Outside,
    NearBoundary,
    Unresolved
  };

  struct PointClassification
  {
    PointState State              = PointState::Unresolved;
    double     DistanceLowerBound = 0; //!< Model length, NOT field residual/distance approximation.
    bool       HasBoundaryWitness = false;
    gp_Pnt WitnessFirst, WitnessLast; //!< A cap point or root bracket containing a boundary point.
    double WitnessDistanceUpperBound = 0;
    double WitnessBracketWidth       = 0;
  };

  struct Regularity
  {
    bool Resolved = false; //!< Numerical regularity/transversality only, NOT solid validity.
    std::array<bool, 4>   Complete{}; //!< Dimensions 0=corners, 1=edges, 2=faces, 3=volume.
    std::array<size_t, 4> Visited{}, Excluded{}, RegularCells{};
    std::array<double, 4>
           MinimumDerivative{};     //!< F/model-length; 0 if no accepted potential-zero cell.
    double HessianBound        = 0; //!< Global operator bound, F/model-length^2.
    int    UnresolvedDimension = -1;
    gp_Pnt UnresolvedMinimum, UnresolvedMaximum;
  };

  enum class CornerKind
  {
    Unresolved,
    Positive,
    Negative,
    PositiveCone,
    NegativeCone,
    SplitCone
  };

  struct CornerGerm
  {
    CornerKind Kind = CornerKind::Unresolved;
    bool ExactZero = false;
    gp_Pnt Anchor, Minimum, Maximum;
    std::array<double, 3> InwardDerivativeLower{}, InwardDerivativeUpper{};
    double Radius = 0; //!< Local neighborhood edge length, NOT an approximation/error bound.
    std::string Reason;
  };

  struct CornerModel
  {
    bool Resolved = false; //!< Local corner neighborhoods ONLY, not strict stratified regularity.
    std::array<CornerGerm, 8> Corners;
    int ExactZeroCount = 0, SplitCount = 0;
    double CoordinateReserve = 0, HessianBound = 0;
  };

  struct CornerMembership
  {
    bool Resolved = false, OnField = false;
    bool ClosedRetained = false, RegularizedRetained = false;
    std::string Reason;
  };

  //! Isolate mesh-free box-corner neighborhoods with uniform-sign F or exact
  //! quarter-period zeros and nonzero signed coordinate derivatives throughout.
  //! Does NOT relax CheckRegularity or enable Build/native solid construction.
  //! Engineering reserves, not directed intervals. Nonrecognized near-zeros,
  //! tangent coordinate directions and exhausted spatial precision stay unresolved.
  Standard_EXPORT static CornerModel CheckCornerGerms(const TopoDS_Face&,
                                                       double theMaximumRadius = 0.1);

  //! Local membership, not a distance-tolerant whole-solid classifier. Reports
  //! closed inequalities separately from closure-of-interior regularization:
  //! an isolated zero corner belongs only to the former. Outside the checked
  //! neighborhood or uncertain field sign returns unresolved, never Outside.
  Standard_EXPORT static CornerMembership ClassifyCornerGerm(const TopoDS_Face&,
    int theCorner, const gp_Pnt&, Operation, double theMaximumRadius = 0.1);

  struct Correspondence
  {
    bool        Resolved = false, RegularityResolved = false;
    size_t      Comparisons = 0, SampleTests = 0;
    int         MatchedEdges = 0, MatchedVertices = 0, MatchedCaps = 0, MatchedWires = 0;
    double      MatchingTolerance = 0, CurveDeviationBound = 0, ReferenceSeamBound = 0;
    double      MinimumTangentAlignment = 0, CapParameterDeviationBound = 0;
    std::string Reason;
  };

  struct TopologyWitness
  {
    bool HomotopyResolved = false, ReferenceManifold = false;
    std::string Reason;
    std::array<size_t, 4> Visited{}, Excluded{}, Transverse{};
    double MinimumHomotopyDerivative = 0, MaximumSimplexDiameter = 0;
    int Resolution = 0, UnresolvedDimension = -1;
    gp_Pnt UnresolvedMinimum, UnresolvedMaximum;
    // Transient nodal-linear reference ONLY, not native topology or a .1mm mesh.
    std::vector<gp_Pnt> Nodes;
    std::vector<std::array<int, 2>> GridEdges; //!< Grid edge IDs; (i,i) is one exact shared corner node.
    std::vector<std::array<int, 3>> Triangles;
    std::vector<std::vector<int>> BoundaryLoops;
    std::vector<int> ComponentEuler;
    CornerPolicy Policy = CornerPolicy::Reject;
    int FixedCornerCount = 0;
    std::vector<gp_Pnt> PointContacts; //!< Isolated zero-set contacts, NOT surface mesh vertices.
  };

  struct RegionTopologyWitness
  {
    TopologyWitness Surface;
    bool ClosedReferenceManifold = false;
    std::string Reason;
    std::vector<gp_Pnt> Nodes; //!< Surface root nodes first, then retained cap grid vertices.
    std::vector<std::array<int, 3>> Triangles;
    std::vector<int> TriangleSides; //!< -1 implicit, otherwise 2*axis+(high side ? 1 : 0).
    std::vector<int> TriangleComponents; //!< Closed-shell reference component, snapshot-local.
    std::vector<int> ComponentEuler;
    std::vector<double> ComponentSignedVolumes; //!< Numerical linear-reference diagnostics ONLY.
    int RegionCount = 0; //!< Connected retained nodal-linear volume regions, not shell count.
    std::vector<int> ComponentRegions; //!< Adjacent retained volume region per closed shell.
    std::vector<double> RegionVolumes; //!< Sum of oriented reference shells per region.
  };

  struct CapAttachment
  {
    bool Resolved = false, StoredGraphResolved = false, RegionReferenceResolved = false;
    std::string Reason;
    size_t WitnessTests = 0, SegmentTests = 0;
    int ReferenceCapCount = 0, ReferenceShellCount = 0, ReferenceRegionCount = 0;
    double RequiredClearance = 0;
    double MatchingTolerance = 0;
    // Source cap order; reference labels are transient, NOT persistent names.
    std::vector<int> ReferenceCaps, ReferenceShells, ReferenceRegions;
    std::vector<gp_Pnt> Witnesses;
    std::vector<double> Clearances;
    RegionTopologyWitness Reference;
    std::vector<int> ReferenceTriangleCaps, ReferenceCapSides;
  };

  struct TrimAttachment
  {
    bool Resolved=false, CapRegionsAttached=false, ComplementCapsAttached=false;
    std::string Reason;
    size_t Comparisons=0, CurveSamples=0;
    int ReferenceArcCount=0, ReferenceWireCount=0, ReferenceSurfaceComponentCount=0;
    std::vector<int> SeamReferenceArcs; //!< Unique implicit-edge traversal order.
    std::vector<int> WireReferenceLoops, WireSurfaceComponents; //!< Native direct wire order.
    CapAttachment RetainedCaps;
    std::vector<int> SurfaceNodeComponents, SurfaceComponentShells, SurfaceComponentRegions;
    std::vector<int> SurfaceComponentEuler;
  };

  struct NativeShells
  {
    bool FacesAttached=false, ShellsConstructed=false, TopologyValidated=false;
    std::string Reason;
    size_t SeparatorCells=0;
    int Separators=0;
    TopoDS_Compound Shape;
    std::vector<TopoDS_Shell> Shells;
    std::vector<TopoDS_Face> Faces;
    std::vector<int> FaceSources, FaceComponents, FaceShells, ShellRegions, ShellEuler;
    std::vector<std::array<gp_Pnt,2>> ComponentDomains;
    Boundary Source; //!< Owned authoritative recipe/topology snapshot; not a mesh definition.
  };

  struct CapGeometry
  {
    bool Resolved=false;
    std::string Reason;
    double MathematicalToCapsBound=0,CapsToMathematicalBound=0,EdgeVertexDeviation=0,BoundaryBound=0;
    std::array<size_t,2> Cells{};
    size_t SegmentTests=0;
    int UnresolvedSide=-1,UnresolvedDirection=-1;
    gp_Pnt UnresolvedMinimum,UnresolvedMaximum;
  };

  //! Whole planar-cap coverage in both directions, using actual lifted p-curves,
  //! Bounded F-value exclusions and explicit retained mathematical cap witnesses.
  //! Includes 3D edge/vertex residuals. No mesh or stored verification flags.
  //! Numerical engineering bounds, NOT directed intervals or topology validation.
  Standard_EXPORT static CapGeometry CheckCapGeometry(const Boundary&,
    const std::vector<TopoDS_Face>& theCaps,double theTolerance=0.1,
    int theMaximumDepth=16,size_t theMaximumCells=2000000,size_t theMaximumSegmentTests=8000000);

  struct NativeSolids
  {
    bool Constructed=false,Validated=false;
    std::string Reason;
    double Tolerance=0.1,MaximumStoredTolerance=0;
    NativeShells Boundary;
    CapGeometry Geometry;
    TopoDS_Compound Shape;
    std::vector<TopoDS_Solid> Solids; //!< Retained material-region order, not shell order.
  };

  //! Fresh native face/shell/topology construction, whole cap/edge/vertex
  //! geometric coverage and bounded stored tolerances gate actual native solids.
  //! Validation is restricted implicit-aware numerical tolerance validity,
  //! NOT generic UV BRepCheck/BOP qualification or a formal interval certificate.
  Standard_EXPORT static NativeSolids BuildNativeSolids(const Boundary&,
    double theTolerance=0.1,int theResolution=32,int theMaximumDepth=16,
    size_t theMaximumCells=4000000,size_t theMaximumTests=2000000,
    size_t theMaximumSeparatorCells=2000000,size_t theMaximumGeometryCells=2000000,
    size_t theMaximumSegmentTests=8000000);

  //! Attach whole field faces, partitioning disconnected components only at
  //! checked F/L-zero-free axis planes, then assemble genuine native shells.
  //! Validate face attachment, oriented incidence, cyclic vertex links, Euler
  //! and component grouping. Never promotes reference triangles to CAD faces.
  //! NOT geometric/solid validation or generic UV-based consumer support.
  Standard_EXPORT static NativeShells BuildNativeShells(const Boundary&,
    int theResolution=32,int theMaximumDepth=16,size_t theMaximumCells=4000000,
    size_t theMaximumTests=2000000,size_t theMaximumSeparatorCells=2000000);

  //! Dual retained/complement cap-region identities and ordered box-edge-root
  //! anchors identify contour arcs. Require a bijection and equal full directed
  //! cyclic native/reference wires. This attaches trim boundaries ONLY, not
  //! the entire implicit face or a native shell/solid; no validity promotion.
  Standard_EXPORT static TrimAttachment CheckTrimAttachment(const Boundary&,
    int theResolution=32, int theMaximumDepth=16, size_t theMaximumCells=4000000,
    size_t theMaximumTests=2000000);

  //! Attach stored conventional cap REGIONS to the closed homotopy reference
  //! using strict common F/L retained points, source face classification and
  //! clearance beyond the independently retraced seam uncertainty. Requires
  //! stored graph correspondence first. NOT a trim-wire/whole-face attachment
  //! or native shell/solid construction/validation certificate.
  Standard_EXPORT static CapAttachment CheckCapAttachment(const Boundary&,
    int theResolution = 32, int theMaximumDepth = 16, size_t theMaximumCells = 4000000,
    size_t theMaximumWitnessTests = 2000000);

  //! Complete the homotopy reference with exact nodal-linear box caps sharing
  //! root/grid-vertex identities. Check closed manifold incidences/vertex links
  //! and report component grouping/oriented volume diagnostics. NOT stored
  //! native trim/cap correspondence or native shell/solid construction.
  Standard_EXPORT static RegionTopologyWitness CheckRegionTopology(const TopoDS_Face&,
    int theResolution = 32, int theMaximumDepth = 16, size_t theMaximumCells = 4000000);

  //! Check the whole BOX zero set against a common Freudenthal nodal-linear
  //! reference. On tetrahedra and true box-boundary strata, exclude homotopy
  //! zeros or prove grad(L).projected-grad(F)>0 for (1-t)F+tL. The continuous
  //! field-gradient direction increases every incident linear piece, including
  //! across internal simplex boundaries. Engineering reserves,
  //! NOT directed intervals. The reference is not the stored trim graph or
  //! preview mesh; no native-shell/solid validation or accuracy claim follows.
  Standard_EXPORT static TopologyWitness CheckTopology(const TopoDS_Face&,
                                                       int theResolution = 32,
                                                       int theMaximumDepth = 16,
                                                       size_t theMaximumCells = 4000000);
  //! Explicit regularized corner reference. Uses an odd grid <= requested
  //! resolution, exact shared corner IDs and separate isolated zero contacts.
  //! Whole-surface homotopy/reference only; no native attachment or solids.
  Standard_EXPORT static TopologyWitness CheckTopology(const TopoDS_Face&,
    int theResolution, int theMaximumDepth, size_t theMaximumCells, CornerPolicy);

  //! Re-trace a PRIVATE canonical box graph; match edge geometry, shared
  //! vertices, cap grouping and directed cyclic wires bijectively. Source
  //! graph/provenance never replaced. This is graph correspondence, NOT
  //! a surface isotopy/reach or native-solid validation certificate.
  Standard_EXPORT static Correspondence CheckCorrespondence(const Boundary&,
                                                            int    theMaximumCells       = 512,
                                                            size_t theMaximumComparisons = 2000000);

  //! Mesh-independent, whole-domain regularity of an identity-placed native
  //! gyroid's box region. Projected gradients are bounded away from zero in
  //! every possible-zero volume/face/edge cell; corners must exclude zero.
  //! Ignores stored trim wires: NOT trim correspondence or topology equivalence.
  //! Engineering floating reserves, NOT directed intervals. Limits unresolved.
  Standard_EXPORT static Regularity CheckRegularity(const TopoDS_Face& theField,
                                                    int                theMaximumDepth = 14,
                                                    size_t             theMaximumCells = 2000000);

  //! Classify the authoritative mathematical region of an unchanged Build or
  //! structurally restored result, not a mesh, approximate cap polygon, or
  //! UV-based solid classifier. Restored seam bounds are not re-certified.
  //! Inside/Outside require a conservative boundary-distance lower bound >
  //! the requested model-length tolerance. NearBoundary requires a cap point
  //! or sign bracket wholly inside that tolerance ball AND the box domain.
  //! Otherwise returns Unresolved, including critical/ambiguous queries.
  //! Ordinary floating-point bounds/reserves, not directed interval arithmetic.
  //! This does not validate boundary topology, prove nonemptiness/components,
  //! or turn the boundary compound into a TopoDS_Solid.
  Standard_EXPORT static PointClassification Classify(const Boundary& theBoundary,
                                                      const gp_Pnt&   thePoint,
                                                      double          theDistanceTolerance = 0.1);

  //! F<=0 inside; Common is box intersect field, BoxMinusField retains F>=0.
  //! Exact gyroid, childless identity-placed forward native face; a valid,
  //! forward, axis-aligned conventional box solid wholly inside its domain.
  //! Transverse box-edge roots are isolated once and shared across faces.
  //! Resolved open/closed contours form multiple/nested planar cap regions.
  //! Tangencies, grid/corner degeneracies and unresolved topology reject.
  //! Seams are shared degree-1 B-splines, with matching planar p-curves.
  //! Returns a boundary compound, NEVER a TopoDS_Solid. Standard UV-based
  //! BRepCheck/BOP/meshing/serialization must not consume the hybrid compound.
  Standard_EXPORT static Boundary Build(const TopoDS_Face&  theField,
                                        const TopoDS_Shape& theBox,
                                        Operation           theOperation,
                                        double              theTolerance    = 0.1,
                                        int                 theMaximumCells = 512);

  //! Explicit corner-aware boundary construction, not a native solid validator.
  //! Regularized material differs from closed inequalities at isolated contacts.
  //! Legacy reference/solid consumers reject this policy
  //! until their corner-aware counterparts are implemented; no silent promotion.
  Standard_EXPORT static Boundary Build(const TopoDS_Face&, const TopoDS_Shape&,
    Operation, double theTolerance, int theMaximumCells, CornerPolicy);
};
#endif
