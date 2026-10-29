// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRepClass3d_ImplicitCylinderIO_HeaderFile
#define _BRepClass3d_ImplicitCylinderIO_HeaderFile
#include <BRepClass3d_ImplicitCylinder.hxx>
#include <Standard_IStream.hxx>
#include <Standard_OStream.hxx>

//! Explicit bounded schema for the qualified three-face cylinder graph.
//! Geometry/topology are stored, not retraced; meshes/callbacks/certificates are not.
class BRepClass3d_ImplicitCylinderIO
{
public:
  struct StoredSolid
  {
    TopoDS_Face Primitive;
    BRep_ImplicitCylinderSection::Settings Parameters;
    bool Negative = true;
    int StoredAxialDirection = 0;
    TopoDS_Solid Shape;
    std::vector<TopoDS_Face> Faces; //!< Cap, cylinder, implicit; shared handles.
    double StoredBoundaryBound = 0, StoredVolumeLower = 0, StoredVolumeUpper = 0, StoredVolumeTolerance = 0;
    bool TopologyValidated = false;
    bool GeometryVerified = false, VolumeVerified = false, SolidValidated = false;
  };
  //! Preflight to <=16 MiB before output. Only the unchanged qualified graph;
  //! supports circles/line/degree-1 spline and <=65536 poles per curve.
  Standard_EXPORT static void Write(const BRepClass3d_ImplicitCylinder::NativeSolid&, Standard_OStream&);
  //! Reject unknown versions, malformed/trailing/oversized records atomically.
  //! Structural topology is checked. Geometry/volume certificates stay false.
  Standard_EXPORT static void Read(Standard_IStream&, StoredSolid&);
};
#endif
