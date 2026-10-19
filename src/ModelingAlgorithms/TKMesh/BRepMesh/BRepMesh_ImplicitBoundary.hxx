// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRepMesh_ImplicitBoundary_HeaderFile
#define _BRepMesh_ImplicitBoundary_HeaderFile
#include <BRepClass3d_ImplicitBox.hxx>
#include <array>
#include <vector>

//! Shared cap constraints and an experimental trim-conforming field mesher.
//! Never sends an implicit face through the UV-based mesher. A candidate
//! mesh is NOT evidence of surface coverage, embeddedness or a native solid.
class BRepMesh_ImplicitBoundary
{
public:
  struct Triangle
  {
    std::array<int, 3> Nodes; //!< Zero-based canonical indices, outward cap orientation.
    int                Cap;   //!< Zero-based index in the native Boundary::Faces.
  };

  struct Result;

  struct Coverage
  {
    bool   Resolved = false; //!< Both directed proximity bounds, NOT topology/solid validity.
    double BoundaryToMeshBound = 0, MeshToBoundaryBound = 0;
    size_t SpatialCells = 0, TrianglePatches = 0;
    gp_Pnt UnresolvedMinimum, UnresolvedMaximum;
    int    UnresolvedTriangle = -1;
    size_t RegularizedFieldCells = 0, RegularizedCapCells = 0;
    size_t CornerWitnesses = 0, ExclusionSlabRejections = 0;
  };

  //! Adaptive whole-box/face coverage and triangle-to-CLIPPED-region witnesses.
  //! Ignores input metadata flags. Numerical reserves, NOT directed intervals,
  //! topology equivalence, contact checks or native-solid validation.
  //! Regularized corners use fresh derivative cones to remove isolated zeros
  //! from the field face/side closure; witness slabs respect face exclusions.
  Standard_EXPORT static Coverage CheckCoverage(const BRepClass3d_ImplicitBox::Boundary&,
                                                const Result&,
                                                double theTolerance      = 0.1,
                                                int    theMaximumDepth   = 14,
                                                size_t theMaximumCells   = 2000000,
                                                size_t theMaximumPatches = 1000000);

  struct MeshAudit
  {
    enum class State
    {
      Clear,
      Contact,
      Limit
    } Status             = State::Limit;
    int    FirstTriangle = -1, SecondTriangle = -1;
    size_t BroadTests = 0, PairTests = 0;
  };

  //! Exact binary64 triangle-contact audit of transient buffers, NOT the
  //! authoritative implicit field or a native solid validator. Triangle IDs
  //! enumerate caps first, then implicit triangles. Shared node IDs define
  //! allowed common simplices; unrelated coincident nodes are contacts.
  //! Clear means no intersections beyond these common simplices. Limit is
  //! unresolved, never clear. Requires coordinate exponent spread <=512.
  Standard_EXPORT static MeshAudit Audit(const Result&,
                                         size_t theMaximumPairTests  = 2000000,
                                         size_t theMaximumBroadTests = 50000000);

  struct Result
  {
    std::vector<gp_Pnt>           Nodes;
    std::vector<Triangle>         CapTriangles;
    std::vector<std::vector<int>> Seams; //!< Native unique edge order; curve-forward samples.
    std::vector<std::vector<int>> Wires; //!< Effective implicit orientation; repeated closing node.
    double SeamGeometricBound  = 0;      //!< Input provenance ONLY, not a surface mesh bound.
    bool   SeamBoundVerified   = false; //!< Preserves input tracing provenance, never re-certifies.
    bool   NativeSurfaceMeshed = false; //!< True for a Tessellate candidate, not Prepare.
    bool   SolidValidated      = false; //!< Never a solid validator.
    std::vector<std::array<int, 3>> ImplicitTriangles; //!< Transient field-surface candidate.
    bool   MeshClosed = false; //!< Paired edges/cyclic vertex links, NOT solid validation.
    bool   SurfaceCoverageVerified = false; //!< Sampling can miss components; always false.
    double FieldDistanceBound      = 0;     //!< Triangle-to-UNTRIMMED-field bound only.
    int    Cells                   = 0;     //!< Accepted longest-axis grid resolution.
    BRepClass3d_ImplicitBox::CornerPolicy Policy = BRepClass3d_ImplicitBox::CornerPolicy::Reject;
    int FixedCornerCount = 0; //!< Exact physical grid-corner identities, not welded nodes.
  };

  //! Mesh a validated private topology snapshot, preserving source caches.
  //! Canonical nodes use topology identity and shared edge-sample ordinal,
  //! NEVER distance-welding unrelated geometry. Every native polyline knot
  //! must occur in cap boundary samples. Cap incidences must be paired except
  //! for the native seams, where effective wire use is exactly opposite.
  //! Limits include a conservative pre-mesh sum of edge poles. Unsupported,
  //! missing/inconsistent samples or exhausted budgets reject atomically.
  Standard_EXPORT static Result Prepare(const BRepClass3d_ImplicitBox::Boundary&,
                                        double theCapDeflection    = 0.1,
                                        int    theMaximumNodes     = 200000,
                                        int    theMaximumTriangles = 400000);

  //! Experimental trim-conforming candidate: uniform tetrahedral extraction
  //! with a constrained collar attached to the exact stored trim samples.
  //! Checks paired/opposite mesh incidence and a triangle-wide one-sided
  //! distance bound to the UNTRIMMED gyroid. Does NOT certify coverage,
  //! embeddedness, clipped-field Hausdorff error or native-solid validity.
  //! Ambiguous loops, collars, critical fields and exhausted limits reject.
  Standard_EXPORT static Result Tessellate(const BRepClass3d_ImplicitBox::Boundary&,
                                           double theTolerance        = 0.1,
                                           int    theMaximumCells     = 256,
                                           int    theMaximumNodes     = 200000,
                                           int    theMaximumTriangles = 400000);
};
#endif
