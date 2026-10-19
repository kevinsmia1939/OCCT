// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepMesh_ImplicitBoundary.hxx>
#include <BRepClass3d_ImplicitBoxIO.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepMesh_Context.hxx>
#include <BRepMesh_EdgeDiscret.hxx>
#include <BRepMesh_Deflection.hxx>
#include <BRepMesh_ShapeTool.hxx>
#include <IMeshTools_CurveTessellator.hxx>
#include <IMeshData_Edge.hxx>
#include <IMeshData_Model.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Poly_PolygonOnTriangulation.hxx>
#include <Poly_Triangulation.hxx>
#include <Standard_DomainError.hxx>
#include <Standard_NotImplemented.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Iterator.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>

namespace
{
using Map = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

void require(bool value, const char* message)
{
  if (!value)
    throw Standard_NotImplemented(message);
}

// The ordinary chordal tessellator may discard short native seam knots using
// the edge's implicit approximation tolerance. These exact planar PWL curves
// need every knot as a constraint, independently of that field-error policy.
class ExactPolyline : public IMeshTools_CurveTessellator
{
public:
  ExactPolyline(const TopoDS_Edge& edge, double maximumStep)
  {
    double first, last;
    myCurve = occ::down_cast<Geom_BSplineCurve>(BRep_Tool::Curve(edge, myLocation, first, last));
    require(!myCurve.IsNull() && myCurve->Degree() == 1, "Unsupported exact cap constraint");
    for (int i = 1; i < myCurve->NbPoles(); ++i)
    {
      const int pieces =
        std::max(1, int(std::ceil(myCurve->Pole(i).Distance(myCurve->Pole(i + 1)) / maximumStep)));
      for (int j = 0; j < pieces; ++j)
        myParameters.push_back(myCurve->Knot(i)
                               + (myCurve->Knot(i + 1) - myCurve->Knot(i)) * (double(j) / pieces));
    }
    myParameters.push_back(last);
  }

  int PointsNb() const override { return int(myParameters.size()); }

  bool Value(int index, gp_Pnt& point, double& parameter) const override
  {
    if (index < 1 || index > PointsNb())
      return false;
    parameter = myParameters[size_t(index - 1)];
    point     = myCurve->Value(parameter).Transformed(myLocation.Transformation());
    return true;
  }

private:
  occ::handle<Geom_BSplineCurve> myCurve;
  TopLoc_Location                myLocation;
  std::vector<double>            myParameters;
};

class ExactBoundaryDiscret : public IMeshTools_ModelAlgo
{
public:
  explicit ExactBoundaryDiscret(double maximumStep)
      : myMaximumStep(maximumStep)
  {
  }

protected:
  bool performInternal(const occ::handle<IMeshData_Model>& model,
                       const IMeshTools_Parameters&        settings,
                       const Message_ProgressRange&) override
  {
    if (model.IsNull())
      return false;
    for (int i = 0; i < model->EdgesNb(); ++i)
    {
      const auto& edge = model->GetEdge(i);
      BRepMesh_Deflection::ComputeDeflection(edge, model->GetMaxSize(), settings);
      for (int p = 0; p < edge->PCurvesNb(); ++p)
        BRepMesh_ShapeTool::CheckAndUpdateFlags(edge, edge->GetPCurve(p));
      BRepMesh_EdgeDiscret::Tessellate3d(edge,
                                         new ExactPolyline(edge->GetEdge(), myMaximumStep),
                                         true);
      BRepMesh_EdgeDiscret::Tessellate2d(edge, true);
    }
    return true;
  }

private:
  double myMaximumStep;
};

using Key = std::pair<int, int>;

Key key(int a, int b)
{
  return std::minmax(a, b);
}

struct Use
{
  int Count = 0, Balance = 0;
};

void use(std::map<Key, Use>& uses, int a, int b)
{
  require(a != b, "Degenerate mesh boundary segment");
  auto& value = uses[key(a, b)];
  ++value.Count;
  value.Balance += a < b ? 1 : -1;
}
} // namespace

BRepMesh_ImplicitBoundary::Result BRepMesh_ImplicitBoundary::Prepare(
  const BRepClass3d_ImplicitBox::Boundary& source,
  double                                   deflection,
  int                                      maximumNodes,
  int                                      maximumTriangles)
{
  if (!std::isfinite(deflection) || deflection < 1.e-7 || maximumNodes < 3 || maximumNodes > 1000000
      || maximumTriangles < 1 || maximumTriangles > 2000000)
    throw Standard_DomainError("Invalid native mesh constraint settings");
  // Existing bounded codec validates the supported representation and gives
  // a wholly private shared cap/wire graph without mutating source caches.
  std::ostringstream record;
  BRepClass3d_ImplicitBoxIO::Write(source, record);
  std::istringstream                input(record.str());
  BRepClass3d_ImplicitBox::Boundary copy;
  BRepClass3d_ImplicitBoxIO::Read(input, copy);
  Result result;
  result.Policy = copy.Policy;
  result.FixedCornerCount = copy.ExactCornerCount;
  result.SeamGeometricBound = source.GeometricBound;
  result.SeamBoundVerified  = source.SeamBoundVerified;
  BRep_Builder    builder;
  TopoDS_Compound caps;
  builder.MakeCompound(caps);
  for (size_t i = 0; i + 1 < copy.Faces.size(); ++i)
    builder.Add(caps, copy.Faces[i]);
  Map edges, vertices, nativeEdges;
  TopExp::MapShapes(caps, TopAbs_EDGE, edges);
  TopExp::MapShapes(caps, TopAbs_VERTEX, vertices);
  TopExp::MapShapes(copy.Faces.back(), TopAbs_EDGE, nativeEdges);
  size_t cost      = 0;
  double magnitude = 1;
  for (const auto& shape : edges)
  {
    TopLoc_Location location;
    double          first, last;
    const auto      curve = occ::down_cast<Geom_BSplineCurve>(
      BRep_Tool::Curve(TopoDS::Edge(shape), location, first, last));
    require(!curve.IsNull() && curve->Degree() == 1, "Unsupported mesh constraint curve");
    ++cost;
    for (int i = 1; i < curve->NbPoles(); ++i)
    {
      const double pieces =
        std::ceil(curve->Pole(i).Distance(curve->Pole(i + 1)) / (deflection * .25));
      require(std::isfinite(pieces) && pieces <= maximumNodes,
              "Native boundary sampling budget exhausted");
      cost += size_t(std::max(1., pieces));
      require(cost <= size_t(maximumNodes), "Native boundary sampling budget exhausted");
    }
    require(cost <= size_t(maximumNodes), "Native boundary pole budget exhausted");
    for (int p = 1; p <= curve->NbPoles(); ++p)
      for (int axis = 1; axis <= 3; ++axis)
        magnitude = std::max(magnitude, std::abs(curve->Pole(p).Coord(axis)));
  }
  const double          reserve = 8192 * std::numeric_limits<double>::epsilon() * magnitude;
  IMeshTools_Parameters settings;
  settings.Deflection               = deflection;
  settings.InParallel               = false;
  settings.InternalVerticesMode     = false;
  settings.ControlSurfaceDeflection = false;
  auto context                      = occ::handle<BRepMesh_Context>(new BRepMesh_Context);
  context->SetEdgeDiscret(new ExactBoundaryDiscret(deflection * .25));
  BRepMesh_IncrementalMesh mesher;
  mesher.SetShape(caps); // conventional caps ONLY
  mesher.ChangeParameters() = settings;
  mesher.Perform(context);
  require(mesher.IsDone() && mesher.GetStatusFlags() == 0, "Conventional cap meshing failed");
  std::vector<int>                 vertexNodes(size_t(vertices.Extent()), -1);
  std::vector<std::vector<int>>    edgeNodes(size_t(edges.Extent()));
  std::vector<std::vector<double>> edgeParameters(size_t(edges.Extent()));
  auto                             add = [&](const gp_Pnt& point) {
    require(result.Nodes.size() < size_t(maximumNodes), "Native mesh node budget exhausted");
    for (int axis = 1; axis <= 3; ++axis)
      require(std::isfinite(point.Coord(axis)), "Nonfinite cap mesh point");
    result.Nodes.push_back(point);
    return int(result.Nodes.size() - 1);
  };
  auto vertexNode = [&](const TopoDS_Vertex& vertex) {
    const int index = vertices.FindIndex(vertex) - 1;
    require(index >= 0, "Missing canonical cap vertex");
    auto& node = vertexNodes[size_t(index)];
    if (node < 0)
      node = add(BRep_Tool::Pnt(vertex));
    return node;
  };
  std::map<Key, Use> capUses;
  for (size_t cap = 0; cap + 1 < copy.Faces.size(); ++cap)
  {
    const auto&     face = copy.Faces[cap];
    TopLoc_Location meshLocation;
    const auto      mesh = BRep_Tool::Triangulation(face, meshLocation);
    require(!mesh.IsNull() && mesh->NbTriangles() > 0, "Missing cap triangulation");
    require(mesh->NbNodes() <= maximumNodes, "Cap mesh node budget exhausted");
    std::vector<int> indices(size_t(mesh->NbNodes() + 1), -1);
    Map              faceEdges;
    TopExp::MapShapes(face, TopAbs_EDGE, faceEdges);
    for (const auto& shape : faceEdges)
    {
      const auto edge  = TopoDS::Edge(shape.Oriented(TopAbs_FORWARD));
      const int  index = edges.FindIndex(edge) - 1;
      require(index >= 0, "Missing shared cap edge");
      const auto polygon = BRep_Tool::PolygonOnTriangulation(edge, mesh, meshLocation);
      require(!polygon.IsNull() && polygon->HasParameters() && polygon->NbNodes() >= 2,
              "Missing cap boundary sample parameters");
      TopLoc_Location curveLocation;
      double          first, last;
      const auto      curve =
        occ::down_cast<Geom_BSplineCurve>(BRep_Tool::Curve(edge, curveLocation, first, last));
      const double parameterReserve = 256 * std::numeric_limits<double>::epsilon()
                                      * (1 + std::max(std::abs(first), std::abs(last)));
      require(std::abs(polygon->Parameter(1) - first) <= parameterReserve
                && std::abs(polygon->Parameter(polygon->NbNodes()) - last) <= parameterReserve,
              "Cap edge samples are not curve-forward");
      auto&      canonical  = edgeNodes[size_t(index)];
      auto&      parameters = edgeParameters[size_t(index)];
      const bool initial    = canonical.empty();
      require(initial || canonical.size() == size_t(polygon->NbNodes()),
              "Adjacent caps disagree on shared edge sampling");
      TopoDS_Vertex a, b;
      TopExp::Vertices(edge, a, b);
      for (int p = 1; p <= polygon->NbNodes(); ++p)
      {
        const int local = polygon->Node(p);
        require(local >= 1 && local <= mesh->NbNodes(), "Invalid cap boundary index");
        const double parameter = polygon->Parameter(p);
        require(std::isfinite(parameter) && (p == 1 || parameter > polygon->Parameter(p - 1)),
                "Unordered cap boundary samples");
        if (initial)
        {
          const int node =
            p == 1 ? vertexNode(a)
            : p == polygon->NbNodes()
              ? vertexNode(b)
              : add(curve->Value(parameter).Transformed(curveLocation.Transformation()));
          canonical.push_back(node);
          parameters.push_back(parameter);
        }
        else
          require(std::abs(parameter - parameters[size_t(p - 1)]) <= parameterReserve,
                  "Adjacent cap sample parameters differ");
        const int  node  = canonical[size_t(p - 1)];
        const auto point = mesh->Node(local).Transformed(meshLocation.Transformation());
        require(point.Distance(result.Nodes[size_t(node)]) <= reserve,
                "Cap mesh boundary differs from canonical topology samples");
        auto& mapped = indices[size_t(local)];
        require(mapped < 0 || mapped == node, "Inconsistent topological corner identity");
        mapped = node;
      }
      for (int knot = 1; knot <= curve->NbKnots(); ++knot)
      {
        const double value = curve->Knot(knot);
        auto         found =
          std::lower_bound(parameters.begin(), parameters.end(), value - parameterReserve);
        require(found != parameters.end() && std::abs(*found - value) <= parameterReserve,
                "Cap mesher omitted a native polyline knot");
      }
    }
    for (int node = 1; node <= mesh->NbNodes(); ++node)
      if (indices[size_t(node)] < 0)
        indices[size_t(node)] = add(mesh->Node(node).Transformed(meshLocation.Transformation()));
    for (int triangle = 1; triangle <= mesh->NbTriangles(); ++triangle)
    {
      int a, b, c;
      mesh->Triangle(triangle).Get(a, b, c);
      a = indices.at(size_t(a));
      b = indices.at(size_t(b));
      c = indices.at(size_t(c));
      if (face.Orientation() == TopAbs_REVERSED)
        std::swap(b, c);
      require(a >= 0 && b >= 0 && c >= 0 && a != b && b != c && c != a,
              "Degenerate cap triangle indices");
      require(gp_Vec(result.Nodes[size_t(a)], result.Nodes[size_t(b)])
                  .Crossed(gp_Vec(result.Nodes[size_t(a)], result.Nodes[size_t(c)]))
                  .Magnitude()
                > reserve * reserve,
              "Degenerate cap triangle area");
      require(result.CapTriangles.size() < size_t(maximumTriangles),
              "Native cap triangle budget exhausted");
      result.CapTriangles.push_back({{a, b, c}, int(cap)});
      use(capUses, a, b);
      use(capUses, b, c);
      use(capUses, c, a);
    }
  }
  for (const auto& edge : nativeEdges)
  {
    const int index = edges.FindIndex(edge) - 1;
    require(index >= 0 && !edgeNodes[size_t(index)].empty(), "Native seam has no cap samples");
    result.Seams.push_back(edgeNodes[size_t(index)]);
  }
  std::map<Key, Use> nativeUses;
  for (TopoDS_Iterator wire(copy.Faces.back()); wire.More(); wire.Next())
  {
    std::vector<int> loop;
    for (TopoDS_Iterator edge(wire.Value()); edge.More(); edge.Next())
    {
      const int index = nativeEdges.FindIndex(edge.Value()) - 1;
      require(index >= 0, "Native wire references an unknown seam");
      auto nodes = result.Seams[size_t(index)];
      if (edge.Value().Orientation() == TopAbs_REVERSED)
        std::reverse(nodes.begin(), nodes.end());
      if (!loop.empty())
        require(loop.back() == nodes.front(), "Disconnected native mesh trim wire");
      for (size_t i = 1; i < nodes.size(); ++i)
        use(nativeUses, nodes[i - 1], nodes[i]);
      loop.insert(loop.end(), nodes.begin() + (loop.empty() ? 0 : 1), nodes.end());
    }
    require(loop.size() >= 4 && loop.front() == loop.back(), "Open native mesh trim wire");
    result.Wires.push_back(std::move(loop));
  }
  for (const auto& entry : capUses)
  {
    const auto found = nativeUses.find(entry.first);
    if (found == nativeUses.end())
      require(entry.second.Count == 2 && entry.second.Balance == 0,
              "Unpaired internal cap mesh edge");
    else
      require(entry.second.Count == 1 && found->second.Count == 1
                && entry.second.Balance + found->second.Balance == 0,
              "Cap/native boundary incidences disagree");
  }
  for (const auto& entry : nativeUses)
    require(capUses.count(entry.first) == 1, "Native mesh seam is missing from caps");
  return result;
}
