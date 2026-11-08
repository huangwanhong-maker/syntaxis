// Geometric graphs: the objects a geometric graph grammar rewrites (docs/theory.md, §2).
//
// A geometric graph is a simple undirected graph whose nodes carry a position in the plane and a symbol: 0 for a terminal (a vertex of
// the data), s > 0 for a nonterminal (a node standing for a pattern of rule s, placed by its frame). Edges carry no labels. Induced
// subgraphs on node sets are what patterns are matched against (Fiser et al. 2016, §3.1).

#pragma once

#include <syntaxis/shape/shape.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <set>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace syntaxis::graphs {

using shape::Point;
using shape::Similarity;

struct Node {
  Point at;
  int symbol = 0;     // 0 terminal; s > 0 the nonterminal of rule s
  Similarity frame;   // for a nonterminal, the map from its rule's canonical frame to the plane
};

struct Graph {
  std::vector<Node> nodes;
  std::vector<std::set<std::size_t>> adjacent;

  std::size_t add(Point at, int symbol = 0, Similarity frame = {}) {
    nodes.push_back({at, symbol, frame});
    adjacent.emplace_back();
    return nodes.size() - 1;
  }
  void connect(std::size_t u, std::size_t v) {
    if (u == v) return;
    if (u >= nodes.size() || v >= nodes.size()) throw std::invalid_argument("syntaxis::graphs::Graph::connect: no such node");
    adjacent[u].insert(v);
    adjacent[v].insert(u);
  }
  bool joined(std::size_t u, std::size_t v) const { return adjacent[u].count(v) > 0; }
  std::size_t size() const { return nodes.size(); }
  std::size_t edge_count() const {
    std::size_t e = 0;
    for (const auto& a : adjacent) e += a.size();
    return e / 2;
  }
  std::vector<std::pair<std::size_t, std::size_t>> edges() const {
    std::vector<std::pair<std::size_t, std::size_t>> out;
    for (std::size_t u = 0; u < nodes.size(); u++)
      for (std::size_t v : adjacent[u])
        if (u < v) out.push_back({u, v});
    return out;
  }
};

// The edges of the subgraph induced on a set of nodes, as pairs of positions in that set.
inline std::vector<std::pair<std::size_t, std::size_t>> induced_edges(const Graph& g, std::span<const std::size_t> set) {
  std::vector<std::pair<std::size_t, std::size_t>> out;
  for (std::size_t i = 0; i < set.size(); i++)
    for (std::size_t j = i + 1; j < set.size(); j++)
      if (g.joined(set[i], set[j])) out.push_back({i, j});
  return out;
}

// The canonical key of the subgraph induced on a set of nodes (positions, symbols and edges, up to similarity).
inline shape::Canonical key_of(const Graph& g, std::span<const std::size_t> set, double tolerance) {
  std::vector<Point> p;
  std::vector<int> label;
  p.reserve(set.size()), label.reserve(set.size());
  for (std::size_t v : set) { p.push_back(g.nodes[v].at); label.push_back(g.nodes[v].symbol); }
  const auto e = induced_edges(g, set);
  return shape::canonical(p, label, e, tolerance);
}

// Whether the subgraph induced on a set is connected.
inline bool connected(const Graph& g, std::span<const std::size_t> set) {
  if (set.empty()) return false;
  std::set<std::size_t> in(set.begin(), set.end()), seen{set[0]};
  std::vector<std::size_t> stack{set[0]};
  while (!stack.empty()) {
    const std::size_t u = stack.back();
    stack.pop_back();
    for (std::size_t v : g.adjacent[u])
      if (in.count(v) && seen.insert(v).second) stack.push_back(v);
  }
  return seen.size() == set.size();
}

// The number of pairs of edges that cross (a plane graph has none).
inline std::size_t crossings(const Graph& g) {
  const auto e = g.edges();
  std::size_t c = 0;
  for (std::size_t i = 0; i < e.size(); i++)
    for (std::size_t j = i + 1; j < e.size(); j++)
      if (shape::segments_cross(g.nodes[e[i].first].at, g.nodes[e[i].second].at, g.nodes[e[j].first].at, g.nodes[e[j].second].at)) c++;
  return c;
}

// ---- examples ---------------------------------------------------------------------------------------------------------

// The Sierpinski graph of a depth: the triangle with corners 0, 1 and e^{iπ/3}, its three corner halves at each depth, nodes shared
// where corners meet (depth 0 is the triangle itself).
inline Graph sierpinski(int depth) {
  Graph g;
  std::vector<std::pair<std::pair<std::int64_t, std::int64_t>, std::size_t>> index;
  const double q = 1e-9;
  const auto node = [&](Point p) {
    const auto k = std::make_pair(static_cast<std::int64_t>(std::llround(p.real() / q)), static_cast<std::int64_t>(std::llround(p.imag() / q)));
    for (const auto& [key, id] : index) if (key == k) return id;
    const std::size_t id = g.add(p);
    index.push_back({k, id});
    return id;
  };
  const std::function<void(Point, Point, Point, int)> go = [&](Point a, Point b, Point c, int d) {
    if (d == 0) {
      const std::size_t u = node(a), v = node(b), w = node(c);
      g.connect(u, v), g.connect(v, w), g.connect(w, u);
      return;
    }
    const Point ab = (a + b) / 2.0, bc = (b + c) / 2.0, ca = (c + a) / 2.0;
    go(a, ab, ca, d - 1), go(ab, b, bc, d - 1), go(ca, bc, c, d - 1);
  };
  go(Point{0, 0}, Point{1, 0}, std::polar(1.0, 3.141592653589793 / 3), depth);
  return g;
}

// A self-similar graph whose copies are joined by edges rather than shared corners: at depth 0 a triangle; at depth d, three copies
// of depth d − 1 shrunk by `ratio` towards the corners of the triangle with corners 0, 1 and e^{iπ/3}, each pair joined by one edge
// between their nearest nodes. A node-replacement grammar reduces it to one node.
inline Graph bridged_triangles(int depth, double ratio = 0.4) {
  Graph g;
  const std::function<std::vector<std::size_t>(Point, Point, Point, int)> go = [&](Point a, Point b, Point c, int d) -> std::vector<std::size_t> {
    if (d == 0) {
      const std::size_t u = g.add(a), v = g.add(b), w = g.add(c);
      g.connect(u, v), g.connect(v, w), g.connect(w, u);
      return {u, v, w};
    }
    const auto toward = [&](Point corner, Point p) { return corner + (p - corner) * ratio; };
    const std::vector<std::size_t> A = go(a, toward(a, b), toward(a, c), d - 1);
    const std::vector<std::size_t> B = go(toward(b, a), b, toward(b, c), d - 1);
    const std::vector<std::size_t> C = go(toward(c, a), toward(c, b), c, d - 1);
    const auto join = [&](const std::vector<std::size_t>& X, const std::vector<std::size_t>& Y) {
      std::size_t bx = X[0], by = Y[0];
      for (std::size_t x : X) for (std::size_t y : Y) if (std::abs(g.nodes[x].at - g.nodes[y].at) < std::abs(g.nodes[bx].at - g.nodes[by].at)) bx = x, by = y;
      g.connect(bx, by);
    };
    join(A, B), join(B, C), join(C, A);
    std::vector<std::size_t> all = A;
    all.insert(all.end(), B.begin(), B.end());
    all.insert(all.end(), C.begin(), C.end());
    return all;
  };
  go(Point{0, 0}, Point{1, 0}, std::polar(1.0, 3.141592653589793 / 3), depth);
  return g;
}

// A grid of n × m nodes, spacing 1, joined to their neighbours.
inline Graph grid(std::size_t n, std::size_t m) {
  Graph g;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < m; j++) g.add(Point{static_cast<double>(j), static_cast<double>(i)});
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < m; j++) {
      if (j + 1 < m) g.connect(i * m + j, i * m + j + 1);
      if (i + 1 < n) g.connect(i * m + j, (i + 1) * m + j);
    }
  return g;
}

// The graph carried by a similarity.
inline Graph transformed(const Graph& g, const Similarity& s) {
  Graph out = g;
  for (auto& nd : out.nodes) { nd.at = s(nd.at); nd.frame = nd.frame.then(s); }
  return out;
}

}  // namespace syntaxis::graphs
