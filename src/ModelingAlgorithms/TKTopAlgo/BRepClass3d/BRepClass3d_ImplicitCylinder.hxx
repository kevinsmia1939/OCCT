// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRepClass3d_ImplicitCylinder_HeaderFile
#define _BRepClass3d_ImplicitCylinder_HeaderFile
#include <BRep_ImplicitCylinderSection.hxx>
#include <string>

//! Restricted native cap-region assembly; NOT a cylinder Boolean/solid builder.
class BRepClass3d_ImplicitCylinder
{
public:
  struct Caps
  {
    TopoDS_Face Field; //!< Authoritative equation/domain snapshot, not replaced by fitted edges.
    BRep_ImplicitCylinderSection::Settings Parameters; //!< Work checkpoint never retained.
    bool Negative = true;
    size_t MaximumNestingTests = 2000000;
    BRep_ImplicitCylinderSection::CapSkeleton Source;
    std::array<std::vector<TopoDS_Face>, 2> Faces;
    std::array<std::vector<TopoDS_Edge>, 2> FieldEdges;
    bool TopologyValidated = false;
    double ContourBound = 0;
    std::string Reason;
  };
  //! Negative selects F<=0; positive selects F>=0. Planar faces use analytic
  //! circles and bounded degree-1 B-spline field seams, not mesh triangles.
  //! Failure returns no faces/field edges. Source field geometry stays immutable.
  //! Cap topology/contour bounds do not certify cylindrical side faces or solids.
  Standard_EXPORT static Caps BuildCaps(const TopoDS_Face&,
    const BRep_ImplicitCylinderSection::Settings&, bool negative = true,
    size_t maximumNestingTests = 2000000);
};
#endif
