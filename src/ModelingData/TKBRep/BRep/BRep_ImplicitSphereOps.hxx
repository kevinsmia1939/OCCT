// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRep_ImplicitSphereOps_HeaderFile
#define _BRep_ImplicitSphereOps_HeaderFile
#include <TopoDS_Face.hxx>
#include <TopoDS_Compound.hxx>
#include <Poly_Triangulation.hxx>
#include <gp_Pln.hxx>
#include <vector>

//! Restricted analytic native-face operations, NOT BRepAlgoAPI/BOP dispatch.
//! Inputs must be forward, childless exact sphere fields with domains enclosing
//! the complete sphere. Arbitrary fields, nested trims, tangency/coincidence and
//! separations within the face tolerance uncertainty band explicitly reject.
class BRep_ImplicitSphereOps
{
public:
  enum class Operation
  {
    Common,
    Fuse,
    Cut
  };

  struct Result
  {
    TopoDS_Compound          Boundary; //!< Faces only, NOT a validated TopoDS_Solid.
    std::vector<TopoDS_Face> Faces;
    std::vector<int>         Sources; //!< 0 = left operand, 1 = right operand.
  };

  struct Mesh
  {
    occ::handle<Poly_Triangulation> Triangulation; //!< Local field coordinates.
    double GeometricBound; //!< Conservative parameter-interpolation distance bound.
  };

  //! Keep the negative signed-distance side of a world-coordinate plane.
  //! Requires a transverse circle, further than tolerance from either pole.
  Standard_EXPORT static TopoDS_Face Trim(const TopoDS_Face&, const gp_Pln&);
  //! Analytic two-sphere boundary split, classification and source provenance.
  //! Transverse results share one exact circle edge with opposite effective use.
  Standard_EXPORT static Result Boolean(const TopoDS_Face&, const TopoDS_Face&, Operation);
  //! Bounded display mesh for full sphere / spherical patch. No generic BRepMesh.
  //! Returned bound applies to this analytic parameterization, not general fields.
  //! Rejects over-budget requests (>120000 vertices or >240000 triangles).
  Standard_EXPORT static Mesh Tessellate(const TopoDS_Face&, double theDeflection);
};
#endif
