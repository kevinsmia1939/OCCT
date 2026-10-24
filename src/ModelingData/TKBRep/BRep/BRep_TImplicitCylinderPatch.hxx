// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRep_TImplicitCylinderPatch_HeaderFile
#define _BRep_TImplicitCylinderPatch_HeaderFile
#include <BRep_TImplicitFace.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Ax2.hxx>

//! Field zero set restricted to a closed finite cylinder in field-local space.
//! A query-domain restriction, NOT sewn trim topology, a cap or a solid.
class BRep_TImplicitCylinderPatch : public BRep_TImplicitFace
{
public:
  Standard_EXPORT BRep_TImplicitCylinderPatch(const occ::handle<Geom_ImplicitField>&,
    const gp_Pnt&, const gp_Pnt&, const gp_Ax2&, double radius,
    double heightMinimum, double heightMaximum, double tolerance);
  const gp_Ax2& Frame() const { return myFrame; }
  double Radius() const { return myRadius; }
  double HeightMinimum() const { return myHeightMinimum; }
  double HeightMaximum() const { return myHeightMaximum; }
  Standard_EXPORT bool Contains(const gp_Pnt&) const;
  //! Childless primitive only; no silent replacement of pre-existing restrictions.
  //! Frame/radius/heights are field-local; both source placements are retained.
  Standard_EXPORT static TopoDS_Face Restrict(const TopoDS_Face&, const gp_Ax2&,
    double radius, double heightMinimum, double heightMaximum);
  Standard_EXPORT occ::handle<TopoDS_TShape> EmptyCopy() const override;
  Standard_EXPORT void DumpJson(Standard_OStream&, int = -1) const override;
  DEFINE_STANDARD_RTTIEXT(BRep_TImplicitCylinderPatch, BRep_TImplicitFace)
private:
  gp_Ax2 myFrame;
  double myRadius, myHeightMinimum, myHeightMaximum;
};
#endif
