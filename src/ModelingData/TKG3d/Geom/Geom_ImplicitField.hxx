// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#ifndef _Geom_ImplicitField_HeaderFile
#define _Geom_ImplicitField_HeaderFile
#include <Standard_Transient.hxx>
#include <Standard_Type.hxx>
#include <Standard_OStream.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <limits>

//! Immutable continuous scalar geometry F(x,y,z)=0. No Geom_Surface/UV map.
//! Root-bracket guarantees require continuity; gradients may be undefined at singularities.
//! Implementations must be safe for concurrent const evaluation.
//! Values/gradients have field units, not necessarily length units.
class Geom_ImplicitField : public Standard_Transient
{
public:
  virtual double Value(const gp_Pnt& thePoint) const    = 0;
  virtual gp_Vec Gradient(const gp_Pnt& thePoint) const = 0;

  virtual double LipschitzBound() const { return std::numeric_limits<double>::infinity(); }

  Standard_EXPORT virtual void DumpJson(Standard_OStream&, int = -1) const;
  DEFINE_STANDARD_RTTIEXT(Geom_ImplicitField, Standard_Transient)
};
#endif
