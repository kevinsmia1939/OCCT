// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRepClass3d_ImplicitWork_HeaderFile
#define _BRepClass3d_ImplicitWork_HeaderFile

#include <Standard_Macro.hxx>

//! Optional per-thread execution control for the implicit algorithms.
//! A scope propagates through the existing nested query/build call graph without
//! changing its geometric inputs, tolerance policies or persisted results.
//! Controls must outlive their scope. Callbacks run on the algorithm's thread;
//! cancellation shared with another thread must use synchronized storage.
class BRepClass3d_ImplicitWork
{
public:
  enum class Phase
  {
    Copying,
    Topology,
    CapAttachment,
    TrimAttachment,
    ShellAssembly,
    Geometry,
    SolidAssembly,
    Meshing,
    Coverage,
    ContactAudit
  };

  class Control
  {
  public:
    virtual ~Control() = default;
    virtual bool IsCancelled() const = 0;
    virtual void EnterPhase(Phase) {}
  };

  //! Nested scopes restore the previous control even when an algorithm throws.
  //! Controls are never installed in other threads or stored in geometry.
  class Scope
  {
  public:
    Standard_EXPORT explicit Scope(Control&);
    Standard_EXPORT ~Scope();
    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;

  private:
    Control* myPrevious;
  };

  //! Throws Standard_AbortiveTransaction when the current control is cancelled.
  //! Without a scope, these calls have no effect on existing synchronous APIs.
  Standard_EXPORT static void Checkpoint();
  Standard_EXPORT static void EnterPhase(Phase);
};

#endif
