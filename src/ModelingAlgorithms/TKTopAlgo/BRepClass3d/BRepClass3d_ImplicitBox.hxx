// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRepClass3d_ImplicitBox_HeaderFile
#define _BRepClass3d_ImplicitBox_HeaderFile
#include <BRep_ImplicitPlaneSection.hxx>
#include <TopoDS_Compound.hxx>
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
    std::vector<std::array<int, 2>> GridEdges; //!< Origin grid-edge IDs for transient root identity.
    std::vector<std::array<int, 3>> Triangles;
    std::vector<std::vector<int>> BoundaryLoops;
    std::vector<int> ComponentEuler;
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
};
#endif
