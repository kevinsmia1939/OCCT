// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRep_TImplicitSpherePatch_HeaderFile
#define _BRep_TImplicitSpherePatch_HeaderFile
#include <BRep_TImplicitFace.hxx>
#include <Geom_ImplicitSphere.hxx>
#include <gp_Pln.hxx>

//! Experimental spherical patch: sphere zero set AND signed Plane distance <= 0.
//! Its circle wire is 3D topology, without UV p-curves. Only explicit implicit
//! consumers support this contract; it is NOT a generic BRepCheck-valid face.
//! Geometry location must remain identity; placement uses the face instance.
class BRep_TImplicitSpherePatch : public BRep_TImplicitFace
{
public:
  Standard_EXPORT BRep_TImplicitSpherePatch(const occ::handle<Geom_ImplicitSphere>&,
                                            const gp_Pln&,
                                            double theTolerance);

  const gp_Pln& Plane() const { return myPlane; }

  Standard_EXPORT occ::handle<TopoDS_TShape> EmptyCopy() const override;
  Standard_EXPORT void                       DumpJson(Standard_OStream&, int = -1) const override;
  DEFINE_STANDARD_RTTIEXT(BRep_TImplicitSpherePatch, BRep_TImplicitFace)
private:
  gp_Pln myPlane;
};
#endif
