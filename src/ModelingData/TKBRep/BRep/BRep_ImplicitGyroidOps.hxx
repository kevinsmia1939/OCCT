// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRep_ImplicitGyroidOps_HeaderFile
#define _BRep_ImplicitGyroidOps_HeaderFile
#include <BRep_ImplicitPlaneSection.hxx>

//! Restricted split seam prototype, NOT generic BOP or a solid constructor.
class BRep_ImplicitGyroidOps
{
public:
  struct Split
  {
    TopoDS_Face Negative, Positive; //!< Authoritative field + opposite plane masks.
    BRep_ImplicitPlaneSection::Result Section;
  };

  //! Identity-placed, childless native gyroid, z strictly inside its domain.
  //! Requires a completely resolved section with CLOSED loops only. The two
  //! patches share degree-1 B-spline seam edges, geometrically bounded
  //! approximations, NOT exact implicit intersection curves. No planar cap,
  //! domain-side wires, mesh, generic Boolean history/solid validity is supplied.
  //! Open/unresolved sections throw before any result is returned.
  Standard_EXPORT static Split SplitZ(const TopoDS_Face&,
                                      double theZ,
                                      double theDistanceTolerance = 0.1,
                                      int    theMaximumCells      = 512);
};
#endif
