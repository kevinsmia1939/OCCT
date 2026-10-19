// Copyright (c) 2026 Kavin Teenakul
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepMesh_ImplicitBoundary.hxx>
#include <Standard_DomainError.hxx>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>

namespace
{
static_assert(std::numeric_limits<double>::radix == 2 && std::numeric_limits<double>::digits == 53
                && std::numeric_limits<double>::max_exponent == 1024,
              "Native exact mesh audit requires binary64 double");

// Small signed arbitrary-width integers. All input doubles are converted to
// integers with a common power-of-two scale; vector predicates then use only
// exact integer operations. No system rounding mode or external dependency.
struct Integer
{
  int                   Sign = 0;
  std::vector<uint32_t> Words;
  Integer() = default;

  Integer(double value, int exponent)
  {
    if (value == 0)
      return;
    Sign = value < 0 ? -1 : 1;
    int            e;
    const uint64_t mantissa = uint64_t(std::ldexp(std::frexp(std::abs(value), &e), 53));
    const int      shift    = e - 53 - exponent;
    Words.resize(size_t(shift / 32) + 3, 0);
    const unsigned offset         = unsigned(shift % 32);
    const uint64_t low            = uint64_t(uint32_t(mantissa)) << offset;
    const uint64_t high           = (mantissa >> 32) << offset;
    Words[size_t(shift / 32)]     = uint32_t(low);
    Words[size_t(shift / 32) + 1] = uint32_t(low >> 32) | uint32_t(high);
    Words[size_t(shift / 32) + 2] = uint32_t(high >> 32);
    trim();
  }

  void trim()
  {
    while (!Words.empty() && Words.back() == 0)
      Words.pop_back();
    if (Words.empty())
      Sign = 0;
  }
};

int magnitude(const Integer& a, const Integer& b)
{
  if (a.Words.size() != b.Words.size())
    return a.Words.size() < b.Words.size() ? -1 : 1;
  for (size_t i = a.Words.size(); i-- > 0;)
    if (a.Words[i] != b.Words[i])
      return a.Words[i] < b.Words[i] ? -1 : 1;
  return 0;
}

Integer operator-(Integer a)
{
  a.Sign = -a.Sign;
  return a;
}

Integer operator+(const Integer& a, const Integer& b)
{
  if (!a.Sign)
    return b;
  if (!b.Sign)
    return a;
  Integer r;
  if (a.Sign == b.Sign)
  {
    r.Sign = a.Sign;
    r.Words.resize(std::max(a.Words.size(), b.Words.size()) + 1, 0);
    uint64_t carry = 0;
    for (size_t i = 0; i + 1 < r.Words.size(); ++i)
    {
      const uint64_t sum =
        carry + (i < a.Words.size() ? a.Words[i] : 0ULL) + (i < b.Words.size() ? b.Words[i] : 0ULL);
      r.Words[i] = uint32_t(sum);
      carry      = sum >> 32;
    }
    r.Words.back() = uint32_t(carry);
  }
  else
  {
    const int comparison = magnitude(a, b);
    if (!comparison)
      return r;
    const auto& x = comparison > 0 ? a : b;
    const auto& y = comparison > 0 ? b : a;
    r.Sign        = x.Sign;
    r.Words.resize(x.Words.size());
    int64_t borrow = 0;
    for (size_t i = 0; i < x.Words.size(); ++i)
    {
      const int64_t value = int64_t(x.Words[i]) - (i < y.Words.size() ? y.Words[i] : 0) - borrow;
      r.Words[i]          = uint32_t(value);
      borrow              = value < 0 ? 1 : 0;
    }
  }
  r.trim();
  return r;
}

Integer operator-(const Integer& a, const Integer& b)
{
  return a + (-b);
}

Integer operator*(const Integer& a, const Integer& b)
{
  Integer r;
  if (!a.Sign || !b.Sign)
    return r;
  r.Sign = a.Sign * b.Sign;
  r.Words.resize(a.Words.size() + b.Words.size(), 0);
  for (size_t i = 0; i < a.Words.size(); ++i)
  {
    uint64_t carry = 0;
    for (size_t j = 0; j < b.Words.size(); ++j)
    {
      const uint64_t value = uint64_t(a.Words[i]) * b.Words[j] + r.Words[i + j] + carry;
      r.Words[i + j]       = uint32_t(value);
      carry                = value >> 32;
    }
    r.Words[i + b.Words.size()] = uint32_t(carry);
  }
  r.trim();
  return r;
}

int compare(const Integer& a, const Integer& b)
{
  return (a - b).Sign;
}

using Vector   = std::array<Integer, 3>;
using Triangle = std::array<int, 3>;

Vector operator-(const Vector& a, const Vector& b)
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Vector cross(const Vector& a, const Vector& b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

Integer dot(const Vector& a, const Vector& b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

bool zero(const Vector& a)
{
  return !a[0].Sign && !a[1].Sign && !a[2].Sign;
}

bool separate(const Vector&              axis,
              const Triangle&            a,
              const Triangle&            b,
              const std::vector<Vector>& nodes,
              const std::vector<int>&    shared)
{
  if (zero(axis))
    return false;
  const Vector& origin = nodes[size_t(shared.empty() ? a[0] : shared[0])];
  for (int node : shared)
    if (dot(axis, nodes[size_t(node)] - origin).Sign != 0)
      return false; // The separating plane must contain the ENTIRE common simplex.
  const auto range = [&](const Triangle& t) {
    std::vector<Integer> values;
    for (int node : t)
      if (std::find(shared.begin(), shared.end(), node) == shared.end())
        values.push_back(dot(axis, nodes[size_t(node)] - origin));
    auto low = values.front(), high = low;
    for (const auto& v : values)
    {
      if (compare(v, low) < 0)
        low = v;
      if (compare(v, high) > 0)
        high = v;
    }
    return std::make_pair(low, high);
  };
  const auto ar = range(a), br = range(b);
  if (shared.empty())
    return compare(ar.second, br.first) < 0 || compare(br.second, ar.first) < 0;
  // A separating plane through the common simplex must leave at least one
  // triangle's remaining vertices STRICTLY on one side. This handles adjacent
  // triangles too; simply skipping pairs sharing a node would miss overlaps.
  return (ar.second.Sign < 0 && br.first.Sign >= 0) || (ar.second.Sign <= 0 && br.first.Sign > 0)
         || (br.second.Sign < 0 && ar.first.Sign >= 0)
         || (br.second.Sign <= 0 && ar.first.Sign > 0);
}

bool clear(const Triangle& a, const Triangle& b, const std::vector<Vector>& nodes)
{
  std::vector<int> shared;
  for (int i : a)
    if (std::find(b.begin(), b.end(), i) != b.end())
      shared.push_back(i);
  if (shared.size() == 3)
    return false;
  Vector ae[3], be[3];
  for (int i = 0; i < 3; ++i)
  {
    ae[i] = nodes[size_t(a[(i + 1) % 3])] - nodes[size_t(a[i])];
    be[i] = nodes[size_t(b[(i + 1) % 3])] - nodes[size_t(b[i])];
  }
  const auto an = cross(ae[0], ae[1]), bn = cross(be[0], be[1]);
  if (separate(an, a, b, nodes, shared) || separate(bn, a, b, nodes, shared))
    return true;
  for (int i = 0; i < 3; ++i)
  {
    if (separate(cross(an, ae[i]), a, b, nodes, shared)
        || separate(cross(bn, be[i]), a, b, nodes, shared))
      return true;
    for (int j = 0; j < 3; ++j)
      if (separate(cross(ae[i], be[j]), a, b, nodes, shared))
        return true;
  }
  return false; // Contact, or no separating-axis witness: never a false clear.
}
} // namespace

BRepMesh_ImplicitBoundary::MeshAudit BRepMesh_ImplicitBoundary::Audit(const Result& mesh,
                                                                      size_t maximumPairTests,
                                                                      size_t maximumBroadTests)
{
  if (mesh.Nodes.empty() || mesh.Nodes.size() > 1000000
      || mesh.CapTriangles.size() + mesh.ImplicitTriangles.size() > 2000000 || !maximumPairTests
      || !maximumBroadTests || maximumPairTests > 10000000 || maximumBroadTests > 100000000)
    throw Standard_DomainError("Invalid native mesh audit settings");
  if (mesh.CapTriangles.empty() && mesh.ImplicitTriangles.empty())
    throw Standard_DomainError("Empty native mesh audit");
  const auto validate = [&](const std::array<int, 3>& t) {
    for (int index : t)
      if (index < 0 || size_t(index) >= mesh.Nodes.size())
        throw Standard_DomainError("Invalid native mesh audit index");
  };
  for (const auto& t : mesh.CapTriangles)
    validate(t.Nodes);
  for (const auto& t : mesh.ImplicitTriangles)
    validate(t);
  MeshAudit result;
  int       minimumExponent = std::numeric_limits<int>::max(), maximumExponent = -minimumExponent;
  for (const auto& p : mesh.Nodes)
    for (int axis = 1; axis <= 3; ++axis)
    {
      const double v = p.Coord(axis);
      if (!std::isfinite(v))
        throw Standard_DomainError("Nonfinite native mesh audit point");
      if (v != 0)
      {
        int e;
        std::frexp(std::abs(v), &e);
        e -= 53;
        minimumExponent = std::min(minimumExponent, e);
        maximumExponent = std::max(maximumExponent, e);
      }
    }
  if (minimumExponent == std::numeric_limits<int>::max())
    minimumExponent = 0;
  else if (maximumExponent - minimumExponent > 512)
    return result;
  std::vector<Vector> nodes;
  for (const auto& p : mesh.Nodes)
    nodes.push_back({Integer(p.X(), minimumExponent),
                     Integer(p.Y(), minimumExponent),
                     Integer(p.Z(), minimumExponent)});
  std::vector<std::array<int, 3>> triangles;
  for (const auto& t : mesh.CapTriangles)
    triangles.push_back(t.Nodes);
  triangles.insert(triangles.end(), mesh.ImplicitTriangles.begin(), mesh.ImplicitTriangles.end());
  if (triangles.empty())
    throw Standard_DomainError("Empty native mesh audit");

  struct Box
  {
    double Low[3], High[3];
  };

  std::vector<Box> boxes(triangles.size());
  for (size_t i = 0; i < triangles.size(); ++i)
  {
    const auto& t = triangles[i];
    for (int n : t)
      if (n < 0 || size_t(n) >= nodes.size())
        throw Standard_DomainError("Invalid native mesh audit index");
    if (zero(cross(nodes[size_t(t[1])] - nodes[size_t(t[0])],
                   nodes[size_t(t[2])] - nodes[size_t(t[0])])))
    {
      result.Status        = MeshAudit::State::Contact;
      result.FirstTriangle = int(i);
      return result;
    }
    for (int axis = 0; axis < 3; ++axis)
    {
      boxes[i].Low[axis]  = std::min({mesh.Nodes[size_t(t[0])].Coord(axis + 1),
                                      mesh.Nodes[size_t(t[1])].Coord(axis + 1),
                                      mesh.Nodes[size_t(t[2])].Coord(axis + 1)});
      boxes[i].High[axis] = std::max({mesh.Nodes[size_t(t[0])].Coord(axis + 1),
                                      mesh.Nodes[size_t(t[1])].Coord(axis + 1),
                                      mesh.Nodes[size_t(t[2])].Coord(axis + 1)});
    }
  }
  std::vector<size_t> order(triangles.size());
  std::iota(order.begin(), order.end(), 0);
  struct Node
  {
    Box Bounds;
    size_t Begin, End;
    int Left = -1, Right = -1;
  };
  std::vector<Node> tree;
  auto build = [&](auto&& self, size_t begin, size_t end) -> int {
    Node node{boxes[order[begin]], begin, end};
    for (size_t i = begin + 1; i < end; ++i)
      for (int axis = 0; axis < 3; ++axis)
      {
        node.Bounds.Low[axis] = std::min(node.Bounds.Low[axis], boxes[order[i]].Low[axis]);
        node.Bounds.High[axis] = std::max(node.Bounds.High[axis], boxes[order[i]].High[axis]);
      }
    const int index = int(tree.size());
    tree.push_back(node);
    if (end - begin > 8)
    {
      int axis = 0;
      for (int a = 1; a < 3; ++a)
        if (node.Bounds.High[a] - node.Bounds.Low[a] > node.Bounds.High[axis] - node.Bounds.Low[axis])
          axis = a;
      const size_t middle = begin + (end - begin) / 2;
      std::nth_element(order.begin() + begin, order.begin() + middle, order.begin() + end,
        [&](size_t a, size_t b) {
          return .5 * boxes[a].Low[axis] + .5 * boxes[a].High[axis]
                 < .5 * boxes[b].Low[axis] + .5 * boxes[b].High[axis];
        });
      const int left = self(self, begin, middle), right = self(self, middle, end);
      tree[size_t(index)].Left = left;
      tree[size_t(index)].Right = right;
    }
    return index;
  };
  build(build, 0, order.size());
  const auto overlaps = [&](const Box& a, const Box& b) {
    for (int axis = 0; axis < 3; ++axis)
      if (a.High[axis] < b.Low[axis] || b.High[axis] < a.Low[axis]) return false;
    return true;
  };
  auto visit = [&](auto&& self, int a, int b) -> bool {
    if (result.BroadTests == maximumBroadTests) return false;
    ++result.BroadTests; // Count node-box AND leaf-pair box comparisons.
    const auto& x = tree[size_t(a)];
    const auto& y = tree[size_t(b)];
    if (!overlaps(x.Bounds, y.Bounds)) return true;
    if (x.Left < 0 && y.Left < 0)
    {
      for (size_t i = x.Begin; i < x.End; ++i)
        for (size_t j = a == b ? i + 1 : y.Begin; j < y.End; ++j)
        {
          if (result.BroadTests == maximumBroadTests) return false;
          ++result.BroadTests;
          const size_t first = order[i], second = order[j];
          if (!overlaps(boxes[first], boxes[second])) continue;
          if (result.PairTests == maximumPairTests) return false;
          ++result.PairTests;
          if (!clear(triangles[first], triangles[second], nodes))
          {
            result.Status = MeshAudit::State::Contact;
            result.FirstTriangle = int(std::min(first, second));
            result.SecondTriangle = int(std::max(first, second));
            return false;
          }
        }
      return true;
    }
    if (a == b)
      return self(self, x.Left, x.Left) && self(self, x.Left, x.Right)
             && self(self, x.Right, x.Right);
    if (x.Left >= 0 && (y.Left < 0 || x.End - x.Begin >= y.End - y.Begin))
      return self(self, x.Left, b) && self(self, x.Right, b);
    return self(self, a, y.Left) && self(self, a, y.Right);
  };
  if (!visit(visit, 0, 0)) return result;
  result.Status = MeshAudit::State::Clear;
  return result;
}
