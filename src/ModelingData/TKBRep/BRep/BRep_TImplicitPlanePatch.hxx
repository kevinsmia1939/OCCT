// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRep_TImplicitPlanePatch_HeaderFile
#define _BRep_TImplicitPlanePatch_HeaderFile
#include <BRep_TImplicitFace.hxx>
#include <gp_Pln.hxx>

//! Field zero set restricted to negative signed plane distance. Explicit
//! implicit consumers only; no UV p-curves or generic BRepCheck validity.
//! Geometry location remains identity; move the face instance instead.
class BRep_TImplicitPlanePatch : public BRep_TImplicitFace
{
public:
  Standard_EXPORT BRep_TImplicitPlanePatch(const occ::handle<Geom_ImplicitField>&,
                                           const gp_Pnt& theMinimum,
                                           const gp_Pnt& theMaximum,
                                           const gp_Pln&,
                                           double theTolerance);

  const gp_Pln& Plane() const { return myPlane; }

  Standard_EXPORT occ::handle<TopoDS_TShape> EmptyCopy() const override;
  Standard_EXPORT void                       DumpJson(Standard_OStream&, int = -1) const override;
  DEFINE_STANDARD_RTTIEXT(BRep_TImplicitPlanePatch, BRep_TImplicitFace)
private:
  gp_Pln myPlane;
};
#endif
