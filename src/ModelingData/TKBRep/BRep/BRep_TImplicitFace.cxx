// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRep_TImplicitFace.hxx>
#include <Standard_DomainError.hxx>
#include <cmath>

IMPLEMENT_STANDARD_RTTIEXT(BRep_TImplicitFace, BRep_TFace)

BRep_TImplicitFace::BRep_TImplicitFace(const occ::handle<Geom_ImplicitField>& f,
                                       const gp_Pnt&                          low,
                                       const gp_Pnt&                          high,
                                       double                                 tol)
    : myField(f),
      myMinimum(low),
      myMaximum(high)
{
  if (f.IsNull() || !std::isfinite(tol) || tol <= 0)
    throw Standard_DomainError(
      "Implicit face requires field and positive finite geometric tolerance");
  for (int i = 1; i <= 3; ++i)
    if (!std::isfinite(low.Coord(i)) || !std::isfinite(high.Coord(i))
        || low.Coord(i) >= high.Coord(i))
      throw Standard_DomainError("Implicit face requires finite non-empty domain");
  Tolerance(tol);
  NaturalRestriction(false);
}

occ::handle<TopoDS_TShape> BRep_TImplicitFace::EmptyCopy() const
{
  occ::handle<BRep_TImplicitFace> copy =
    new BRep_TImplicitFace(myField, myMinimum, myMaximum, Tolerance());
  copy->Location(Location());
  return copy;
}

void BRep_TImplicitFace::DumpJson(Standard_OStream& os, int depth) const
{
  OCCT_DUMP_TRANSIENT_CLASS_BEGIN(os)
  OCCT_DUMP_BASE_CLASS(os, depth, BRep_TFace)
  OCCT_DUMP_FIELD_VALUES_DUMPED(os, depth, myField.get())
  OCCT_DUMP_FIELD_VALUES_DUMPED(os, depth, &myMinimum)
  OCCT_DUMP_FIELD_VALUES_DUMPED(os, depth, &myMaximum)
}
