// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRepClass3d_ImplicitBoxIO_HeaderFile
#define _BRepClass3d_ImplicitBoxIO_HeaderFile
#include <BRepClass3d_ImplicitBox.hxx>
#include <Standard_IStream.hxx>
#include <Standard_OStream.hxx>

//! Explicit schema-1/2 hybrid topology I/O. Not generic BRepTools/OCAF/STEP.
//! Schema 2 retains regularized corner policy and checked contact identities.
//! Persists the authoritative primitive and shared conventional cap/edge/wire
//! graph, never a substitute UV surface, display mesh, pointer or recipe rebuild.
//! Structural validation is NOT whole-solid/geometric-error certification.
class BRepClass3d_ImplicitBoxIO
{
public:
  //! Preflight into a bounded temporary record before touching output. Supports
  //! only unchanged native gyroid/axis-box boundaries, <=128 faces, <=16 MiB.
  Standard_EXPORT static void Write(const BRepClass3d_ImplicitBox::Boundary&, Standard_OStream&);
  //! Complete-record read; rejects unknown schemas, malformed/trailing/oversized
  //! data. Destination unchanged on failure. Shared identities/p-curves survive.
  //! Stored approximation metadata is NOT re-certified: SeamBoundVerified=false.
  //! No meshes, shape flags or original tracing diagnostics are restored.
  Standard_EXPORT static void Read(Standard_IStream&, BRepClass3d_ImplicitBox::Boundary&);
};
#endif
