// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRepClass3d_ImplicitSolidIO_HeaderFile
#define _BRepClass3d_ImplicitSolidIO_HeaderFile
#include <BRepClass3d_ImplicitBox.hxx>
#include <Standard_IStream.hxx>
#include <Standard_OStream.hxx>

//! Versioned persistence of actual restricted native solid topology.
//! Reuses the authoritative shared source cap/edge/wire graph, with explicit
//! field-component domains, wire assignments and face/shell/material relations.
//! Never rebuilds numerically, persists meshes or restores validation flags.
//! Schema 1 keeps strict sources. Schema 2 retains regularized corner sources
//! and derives component restrictions canonically from source/domain metadata.
class BRepClass3d_ImplicitSolidIO
{
public:
  //! Preflight a complete bounded record before touching output. Supported
  //! factory topology only; no arbitrary solid/field serialization or STEP.
  Standard_EXPORT static void Write(const BRepClass3d_ImplicitBox::NativeSolids&,
                                    Standard_OStream&);
  //! Atomic complete-record read, structurally checked but NOT mathematically
  //! revalidated. Constructed=true, Validated=false, shell/solid Closed=false.
  //! Shared native topology is restored, not regenerated from a Boolean recipe.
  Standard_EXPORT static void Read(Standard_IStream&, BRepClass3d_ImplicitBox::NativeSolids&);
};
#endif
