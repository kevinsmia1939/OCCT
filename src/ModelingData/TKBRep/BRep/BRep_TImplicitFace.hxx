// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#ifndef _BRep_TImplicitFace_HeaderFile
#define _BRep_TImplicitFace_HeaderFile
#include <TopoDS_Shape.hxx>
#include <BRep_TFace.hxx>
#include <Geom_ImplicitField.hxx>

//! Experimental native field-backed face storage, without a fake UV surface.
//! The finite domain clips queries; it is NOT sewn trimming topology or caps.
//! Only the explicit BRep_ImplicitTool operations are supported at this stage.
class BRep_TImplicitFace : public BRep_TFace
{
public:
  Standard_EXPORT BRep_TImplicitFace(const occ::handle<Geom_ImplicitField>& theField,
                                     const gp_Pnt&                          theMinimum,
                                     const gp_Pnt&                          theMaximum,
                                     double                                 theTolerance);

  const occ::handle<Geom_ImplicitField>& Field() const { return myField; }

  const gp_Pnt& Minimum() const { return myMinimum; }

  const gp_Pnt& Maximum() const { return myMaximum; }

  Standard_EXPORT occ::handle<TopoDS_TShape> EmptyCopy() const override;
  Standard_EXPORT void                       DumpJson(Standard_OStream&, int = -1) const override;
  DEFINE_STANDARD_RTTIEXT(BRep_TImplicitFace, BRep_TFace)
private:
  occ::handle<Geom_ImplicitField> myField;
  gp_Pnt                          myMinimum, myMaximum;
};
#endif
