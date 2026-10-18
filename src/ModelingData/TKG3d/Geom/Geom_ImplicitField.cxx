// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <Geom_ImplicitField.hxx>
#include <Geom_ImplicitSphere.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Standard_DomainError.hxx>
#include <cmath>

IMPLEMENT_STANDARD_RTTIEXT(Geom_ImplicitField, Standard_Transient)
IMPLEMENT_STANDARD_RTTIEXT(Geom_ImplicitSphere, Geom_ImplicitField)
IMPLEMENT_STANDARD_RTTIEXT(Geom_ImplicitGyroid, Geom_ImplicitField)

void Geom_ImplicitField::DumpJson(Standard_OStream& os, int) const
{
  os << "{\"type\":\"" << DynamicType()->Name() << "\"}";
}

namespace
{
void checkPoint(const gp_Pnt& p)
{
  if (!std::isfinite(p.X()) || !std::isfinite(p.Y()) || !std::isfinite(p.Z()))
    throw Standard_DomainError("Implicit geometry requires finite coordinates");
}

double checked(double v)
{
  if (!std::isfinite(v))
    throw Standard_DomainError("Implicit geometry evaluation overflow");
  return v;
}
} // namespace

Geom_ImplicitSphere::Geom_ImplicitSphere(const gp_Pnt& c, double r)
    : myCenter(c),
      myRadius(r)
{
  checkPoint(c);
  if (!std::isfinite(r) || r <= 0)
    throw Standard_DomainError("Implicit sphere requires positive finite radius");
}

double Geom_ImplicitSphere::Value(const gp_Pnt& p) const
{
  checkPoint(p);
  return checked(std::hypot(p.X() - myCenter.X(), p.Y() - myCenter.Y(), p.Z() - myCenter.Z())
                 - myRadius);
}

gp_Vec Geom_ImplicitSphere::Gradient(const gp_Pnt& p) const
{
  checkPoint(p);
  gp_Vec g(myCenter, p);
  double n = checked(std::hypot(g.X(), g.Y(), g.Z()));
  return n > 0 ? g / n : gp_Vec(0, 0, 0);
}

Geom_ImplicitGyroid::Geom_ImplicitGyroid(double p, double c)
    : myPeriod(p),
      myIsoValue(c),
      myWaveNumber(0)
{
  if (!std::isfinite(p) || p <= 0 || !std::isfinite(c))
    throw Standard_DomainError("Implicit gyroid requires valid finite period and level");
  myWaveNumber = 2 * std::acos(-1.0) / p;
  if (!std::isfinite(myWaveNumber) || myWaveNumber <= 0)
    throw Standard_DomainError("Implicit gyroid requires valid finite period and level");
}

double Geom_ImplicitGyroid::Value(const gp_Pnt& p) const
{
  checkPoint(p);
  double x = checked(myWaveNumber * p.X()), y = checked(myWaveNumber * p.Y()),
         z = checked(myWaveNumber * p.Z());
  return checked(std::sin(x) * std::cos(y) + std::sin(y) * std::cos(z) + std::sin(z) * std::cos(x)
                 - myIsoValue);
}

gp_Vec Geom_ImplicitGyroid::Gradient(const gp_Pnt& p) const
{
  checkPoint(p);
  double x = checked(myWaveNumber * p.X()), y = checked(myWaveNumber * p.Y()),
         z = checked(myWaveNumber * p.Z());
  return gp_Vec(checked(myWaveNumber * (std::cos(x) * std::cos(y) - std::sin(z) * std::sin(x))),
                checked(myWaveNumber * (std::cos(y) * std::cos(z) - std::sin(x) * std::sin(y))),
                checked(myWaveNumber * (std::cos(z) * std::cos(x) - std::sin(y) * std::sin(z))));
}

double Geom_ImplicitGyroid::LipschitzBound() const
{
  return 2 * std::sqrt(3.0) * myWaveNumber;
}

bool Geom_ImplicitGyroid::LatticeJet(const gp_Pnt& p, double& value, gp_Vec& gradient) const
{
  checkPoint(p);
  int sine[3], cosine[3];
  constexpr int s[4] = {0, 1, 0, -1}, c[4] = {1, 0, -1, 0};
  int exponent = 0;
  const double period = std::frexp(myPeriod, &exponent);
  for (int axis = 0; axis < 3; ++axis)
  {
    const double coordinate = p.Coord(axis + 1);
    const double multiple = std::round((coordinate / myPeriod) * 4);
    if (!std::isfinite(multiple) || std::abs(multiple) > 1000000
        || (multiple == 0 && coordinate != 0))
      return false;
    const double four = 4 * std::scalbn(coordinate, -exponent);
    if (!std::isfinite(four) || std::fma(-multiple, period, four) != 0)
      return false;
    // Normalize by an exact power of two first: nonzero candidate coordinates
    // are comparable to period, and a nonzero residual is normal, not a tiny
    // original-coordinate product that could underflow to a false identity.
    const int quadrant = (int(multiple) % 4 + 4) % 4;
    sine[axis] = s[quadrant];
    cosine[axis] = c[quadrant];
  }
  const double candidate = sine[0] * cosine[1] + sine[1] * cosine[2]
                         + sine[2] * cosine[0] - myIsoValue;
  const gp_Vec derivative(myWaveNumber * (cosine[0] * cosine[1] - sine[2] * sine[0]),
                          myWaveNumber * (cosine[1] * cosine[2] - sine[0] * sine[1]),
                          myWaveNumber * (cosine[2] * cosine[0] - sine[1] * sine[2]));
  value = candidate;
  gradient = derivative;
  return true;
}
