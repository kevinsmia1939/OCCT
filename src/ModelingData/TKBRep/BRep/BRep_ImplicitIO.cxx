// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRep_ImplicitIO.hxx>
#include <BRep_Builder.hxx>
#include <BRep_TImplicitFace.hxx>
#include <Geom_ImplicitSphere.hxx>
#include <Geom_ImplicitGyroid.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_LockedShape.hxx>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <typeinfo>

namespace
{
constexpr size_t MaxBytes = 65536;

void require(bool condition)
{
  if (!condition)
    throw Standard_DomainError("BRep_ImplicitIO: invalid schema-1 record or stream failure");
}

void tag(std::istream& in, const char* expected)
{
  std::string token;
  require(bool(in >> token) && token == expected);
}

double number(std::istream& in)
{
  double value = 0;
  require(bool(in >> value) && std::isfinite(value));
  return value;
}

gp_Pnt point(std::istream& in)
{
  // Do not depend on function-argument evaluation order.
  const double x = number(in), y = number(in), z = number(in);
  return gp_Pnt(x, y, z);
}

void writePoint(std::ostream& out, const gp_Pnt& p)
{
  require(std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()));
  out << p.X() << ' ' << p.Y() << ' ' << p.Z() << '\n';
}

gp_Trsf readTransform(std::istream& in)
{
  double a[12];
  for (double& v : a)
    v = number(in);
  double magnitude = 0;
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      magnitude = std::max(magnitude, std::abs(a[r * 4 + c]));
  require(magnitude > 0);
  double b[12];
  for (int i = 0; i < 12; ++i)
    b[i] = i % 4 == 3 ? a[i] : a[i] / magnitude;
  // Validate BEFORE gp_Trsf arithmetic (Release OCCT disables some checks).
  // Normalization also avoids overflow/underflow in its determinant calculation.
  double norm = 0;
  for (int r = 0; r < 3; ++r)
    norm += b[r * 4] * b[r * 4];
  for (int c = 0; c < 3; ++c)
    for (int d = c; d < 3; ++d)
    {
      double dot = 0;
      for (int r = 0; r < 3; ++r)
        dot += b[r * 4 + c] * b[r * 4 + d];
      require(std::abs(dot - (c == d ? norm : 0)) <= 1.e-12 * norm);
    }
  gp_Trsf t;
  t.SetValues(b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11]);
  t.SetScaleFactor(t.ScaleFactor() * magnitude);
  // SetValues orthogonalizes its input: never silently accept a shear or
  // nonuniform scale. Compare each linear coefficient against the input scale.
  const double scale = std::abs(t.ScaleFactor());
  require(std::isfinite(scale) && scale > 0 && std::isfinite(1.0 / scale));
  for (int r = 1; r <= 3; ++r)
    for (int c = 1; c <= 4; ++c)
    {
      const double v = t.Value(r, c);
      require(std::isfinite(v));
      const double bound = c == 4 ? std::max(1.0, std::abs(v)) : scale;
      require(std::abs(v - a[(r - 1) * 4 + c - 1]) <= 1.e-12 * bound);
    }
  return t;
}

void writeTransform(std::ostream& out, const TopLoc_Location& location)
{
  const gp_Trsf& t = location.Transformation();
  for (int r = 1; r <= 3; ++r)
    for (int c = 1; c <= 4; ++c)
    {
      require(std::isfinite(t.Value(r, c)));
      out << t.Value(r, c) << (c == 4 ? '\n' : ' ');
    }
}
} // namespace

void BRep_ImplicitIO::Write(const TopoDS_Face& face, Standard_OStream& stream)
{
  const auto data = occ::down_cast<BRep_TImplicitFace>(face.TShape());
  if (data.IsNull() || typeid(*data) != typeid(BRep_TImplicitFace) || TopoDS_Iterator(face).More())
    throw Standard_NotImplemented("BRep_ImplicitIO: requires a childless native implicit face");

  std::ostringstream out;
  out.imbue(std::locale::classic());
  out << std::setprecision(std::numeric_limits<double>::max_digits10);
  out << "OCCT_IMPLICIT_FACE 1\n";
  const auto& field = data->Field();
  // Exact C++ type, not IsKind(): derived classes may override the equation.
  if (typeid(*field) == typeid(Geom_ImplicitSphere))
  {
    const auto* sphere = static_cast<const Geom_ImplicitSphere*>(field.get());
    out << "SPHERE ";
    writePoint(out, sphere->Center());
    out << sphere->Radius() << '\n';
  }
  else if (typeid(*field) == typeid(Geom_ImplicitGyroid))
  {
    const auto* gyroid = static_cast<const Geom_ImplicitGyroid*>(field.get());
    out << "GYROID " << gyroid->Period() << ' ' << gyroid->IsoValue() << '\n';
  }
  else
    throw Standard_NotImplemented("BRep_ImplicitIO: unsupported field type");

  out << "DOMAIN ";
  writePoint(out, data->Minimum());
  writePoint(out, data->Maximum());
  out << "TOLERANCE " << data->Tolerance() << '\n';
  out << "GEOMETRY_LOCATION\n";
  writeTransform(out, data->Location());
  out << "INSTANCE_LOCATION\n";
  writeTransform(out, face.Location());
  out << "ORIENTATION " << int(face.Orientation()) << "\nEND\n";

  // Use the same validation as the reader before touching the caller's stream.
  TopoDS_Face        check;
  std::istringstream preflight(out.str());
  Read(preflight, check);
  const std::string record = out.str();
  stream.write(record.data(), std::streamsize(record.size()));
  require(bool(stream));
}

void BRep_ImplicitIO::Read(Standard_IStream& stream, TopoDS_Face& face)
{
  if (!face.IsNull() && face.Locked())
    throw TopoDS_LockedShape("BRep_ImplicitIO: destination is locked");
  std::string record;
  char        ch;
  while (stream.get(ch))
  {
    require(record.size() < MaxBytes);
    record.push_back(ch);
  }
  require(stream.eof() && !stream.bad());
  std::istringstream in(record);
  in.imbue(std::locale::classic());
  tag(in, "OCCT_IMPLICIT_FACE");
  tag(in, "1");
  std::string kind;
  require(bool(in >> kind));
  occ::handle<Geom_ImplicitField> field;
  if (kind == "SPHERE")
  {
    const gp_Pnt center = point(in);
    field               = new Geom_ImplicitSphere(center, number(in));
  }
  else if (kind == "GYROID")
  {
    const double period = number(in), iso = number(in);
    field = new Geom_ImplicitGyroid(period, iso);
  }
  else
    throw Standard_NotImplemented("BRep_ImplicitIO: unsupported field type");
  tag(in, "DOMAIN");
  const gp_Pnt low = point(in), high = point(in);
  tag(in, "TOLERANCE");
  const double tolerance = number(in);
  tag(in, "GEOMETRY_LOCATION");
  const TopLoc_Location geometry(readTransform(in));
  tag(in, "INSTANCE_LOCATION");
  const TopLoc_Location instance(readTransform(in));
  tag(in, "ORIENTATION");
  std::string orientation;
  require(bool(in >> orientation));
  require(orientation == "0" || orientation == "1" || orientation == "2" || orientation == "3");
  tag(in, "END");
  in >> std::ws;
  require(in.eof());

  TopoDS_Face result;
  BRep_Builder().MakeImplicitFace(result, field, low, high, geometry, tolerance);
  result.Location(instance, false); // Uniform scale/reflection is supported.
  result.Orientation(TopAbs_Orientation(orientation[0] - '0'));
  face = result;
}
