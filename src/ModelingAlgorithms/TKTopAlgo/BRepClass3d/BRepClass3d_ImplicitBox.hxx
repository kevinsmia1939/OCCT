// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRepClass3d_ImplicitBox_HeaderFile
#define _BRepClass3d_ImplicitBox_HeaderFile
#include <BRep_ImplicitPlaneSection.hxx>
#include <TopoDS_Compound.hxx>
#include <vector>

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
    double GeometricBound = 0; //!< Seam bound ONLY, not a surface/solid certificate.
  };

  //! F<=0 inside; Common is box intersect field, BoxMinusField retains F>=0.
  //! Exact gyroid, childless identity-placed forward native face; a valid,
  //! forward, axis-aligned conventional box solid wholly inside its domain.
  //! All box edges must be certifiably away from F=0. Each planar face may
  //! contain at most one resolved closed contour (no edge-crossing/nesting).
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
