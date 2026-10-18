# Source files for BRepClass3d package
set(OCCT_BRepClass3d_FILES_LOCATION "${CMAKE_CURRENT_LIST_DIR}")

set(OCCT_BRepClass3d_FILES
  BRepClass3d.cxx
  BRepClass3d.hxx
  BRepClass3d_BndBoxTree.hxx
  BRepClass3d_BndBoxTree.cxx

  BRepClass3d_Intersector3d.cxx
  BRepClass3d_Intersector3d.hxx
  BRepClass3d_Intersector3d.lxx

  BRepClass3d_SClassifier.cxx
  BRepClass3d_SClassifier.hxx
  BRepClass3d_SolidClassifier.cxx
  BRepClass3d_SolidClassifier.hxx
  BRepClass3d_SignedDistance.cxx
  BRepClass3d_SignedDistance.hxx
  BRepClass3d_ImplicitBox.hxx
  BRepClass3d_ImplicitBox.cxx
  BRepClass3d_ImplicitBoxRegularity.cxx
  BRepClass3d_ImplicitBoxTopology.cxx
  BRepClass3d_ImplicitBoxCaps.cxx
  BRepClass3d_ImplicitBoxShells.cxx
  BRepClass3d_ImplicitBoxGeometry.cxx
  BRepClass3d_ImplicitBoxIO.hxx
  BRepClass3d_ImplicitBoxIO.cxx
  BRepClass3d_ImplicitSolidIO.hxx
  BRepClass3d_ImplicitSolidIO.cxx
  BRepClass3d_SolidExplorer.cxx
  BRepClass3d_SolidExplorer.hxx
)
