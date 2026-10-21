// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepClass3d_ImplicitWork.hxx>
#include <Standard_AbortiveTransaction.hxx>

namespace
{
thread_local BRepClass3d_ImplicitWork::Control* currentControl = nullptr;
}

BRepClass3d_ImplicitWork::Scope::Scope(Control& control)
  : myPrevious(currentControl)
{
  currentControl = &control;
}

BRepClass3d_ImplicitWork::Scope::~Scope()
{
  currentControl = myPrevious;
}

void BRepClass3d_ImplicitWork::Checkpoint()
{
  if (currentControl && currentControl->IsCancelled())
    throw Standard_AbortiveTransaction("Native implicit work cancelled");
}

void BRepClass3d_ImplicitWork::EnterPhase(Phase phase)
{
  Checkpoint();
  if (currentControl)
    currentControl->EnterPhase(phase);
  Checkpoint();
}
