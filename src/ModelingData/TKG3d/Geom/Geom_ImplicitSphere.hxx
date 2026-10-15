// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#ifndef _Geom_ImplicitSphere_HeaderFile
#define _Geom_ImplicitSphere_HeaderFile
#include <Geom_ImplicitField.hxx>

//! Exact signed-distance sphere, negative inside; immutable parameters.
class Geom_ImplicitSphere : public Geom_ImplicitField
{
public:
  Standard_EXPORT        Geom_ImplicitSphere(const gp_Pnt& theCenter, double theRadius);
  Standard_EXPORT double Value(const gp_Pnt& thePoint) const override;
  Standard_EXPORT gp_Vec Gradient(const gp_Pnt& thePoint) const override;

  double LipschitzBound() const override { return 1.0; }

  const gp_Pnt& Center() const { return myCenter; }

  double Radius() const { return myRadius; }
  DEFINE_STANDARD_RTTIEXT(Geom_ImplicitSphere, Geom_ImplicitField)
private:
  gp_Pnt myCenter;
  double myRadius;
};
#endif
