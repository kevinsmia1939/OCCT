// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRepMesh_ImplicitBoundary_HeaderFile
#define _BRepMesh_ImplicitBoundary_HeaderFile
#include <BRepClass3d_ImplicitBox.hxx>
#include <array>
#include <vector>

//! Shared boundary constraints for a FUTURE native implicit surface mesher.
//! Only conventional planar caps are triangulated. Never sends an implicit
//! face through the UV-based mesher and never claims a complete mesh/solid.
class BRepMesh_ImplicitBoundary
{
public:
  struct Triangle
  {
    std::array<int, 3> Nodes; //!< Zero-based canonical indices, outward cap orientation.
    int                Cap;   //!< Zero-based index in the native Boundary::Faces.
  };

  struct Result
  {
    std::vector<gp_Pnt>           Nodes;
    std::vector<Triangle>         CapTriangles;
    std::vector<std::vector<int>> Seams; //!< Native unique edge order; curve-forward samples.
    std::vector<std::vector<int>> Wires; //!< Effective implicit orientation; repeated closing node.
    double SeamGeometricBound  = 0;      //!< Input provenance ONLY, not a surface mesh bound.
    bool   SeamBoundVerified   = false; //!< Preserves input tracing provenance, never re-certifies.
    bool   NativeSurfaceMeshed = false; //!< Preparation only, always false.
    bool   SolidValidated      = false; //!< Never a solid validator.
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
};
#endif
