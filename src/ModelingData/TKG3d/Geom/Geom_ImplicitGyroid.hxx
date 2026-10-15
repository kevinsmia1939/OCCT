// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#ifndef _Geom_ImplicitGyroid_HeaderFile
#define _Geom_ImplicitGyroid_HeaderFile
#include <Geom_ImplicitField.hxx>

//! Dimensionless gyroid level set. Value/tolerance are NOT geometric distances.
class Geom_ImplicitGyroid : public Geom_ImplicitField
{
public:
  Standard_EXPORT        Geom_ImplicitGyroid(double thePeriod, double theIsoValue = 0.0);
  Standard_EXPORT double Value(const gp_Pnt& thePoint) const override;
  Standard_EXPORT gp_Vec Gradient(const gp_Pnt& thePoint) const override;
  Standard_EXPORT double LipschitzBound() const override;

  double Period() const { return myPeriod; }

  double IsoValue() const { return myIsoValue; }
  DEFINE_STANDARD_RTTIEXT(Geom_ImplicitGyroid, Geom_ImplicitField)
private:
  double myPeriod, myIsoValue, myWaveNumber;
};
#endif
