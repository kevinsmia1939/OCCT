// Copyright (c) 2026 Kavin Teenakul
// This file is part of Open CASCADE Technology software library.
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#ifndef _BRep_ImplicitIO_HeaderFile
#define _BRep_ImplicitIO_HeaderFile
#include <Standard_IStream.hxx>
#include <Standard_OStream.hxx>
#include <Standard_Macro.hxx>
#include <TopoDS_Face.hxx>

//! Experimental, versioned single-field-face interchange, NOT BRepTools/OCAF.
//! Schema 1 stores sphere/gyroid parameters, local domain, raw tolerance,
//! geometry/instance transformations and orientation. Meshes and shape flags
//! are transient; location datum sharing is not preserved. Storage subtypes,
//! child topology and
//! custom fields are rejected rather than silently lost. Units are caller-owned.
class BRep_ImplicitIO
{
public:
  //! Preflight before writing; unsupported input leaves the stream untouched.
  //! Throws Standard_NotImplemented for unsupported fields/topology and
  //! Standard_DomainError for invalid data or stream failure.
  //! Standard stream exceptions propagate if enabled by the caller.
  Standard_EXPORT static void Write(const TopoDS_Face&, Standard_OStream&);

  //! Read one complete record through EOF (maximum 64 KiB), rejecting trailing
  //! data and unknown schemas. On failure the destination face is unchanged.
  //! A locked destination is rejected. Triangulations are never restored.
  Standard_EXPORT static void Read(Standard_IStream&, TopoDS_Face&);
};
#endif
