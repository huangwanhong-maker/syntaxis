// Derivations: a grammar's productions applied, to rebuild what was learned or to make something new (docs/theory.md, §3–4).
//
// - replay: the input rebuilt from a learned derivation. Its edges are exactly the input's (every contraction recorded the input
//   edges it swallowed, and the axiom's edges carry the rest). Its positions are the rules' canonical patterns carried by each
//   instance's frame, down to the terminals: equal to the input's where the occurrences of each rule are exactly similar to its first,
//   and within the tolerance otherwise.
// - expand: the grammar applied to a start graph whose nodes may be nonterminals: each nonterminal is replaced by its rule's pattern,
//   placed by the node's frame; the pattern's own edges join its nodes, and each edge that reached the replaced node is joined to the
//   pattern's node nearest to the edge's other end — the relaxed embedding of Fiser et al. (2016), §3.3, which always rewrites and
//   connects to what lies nearest. Nonterminals inside a pattern are placed by their frames relative to it, and expanded in turn.
//
// Limits: the relaxed embedding does not know which node of a pattern an edge reached in the data, so expand rebuilds a learned
// input's nodes but only approximately its edges; replay is exact.

#pragma once

#include <syntaxis/graphs/graphs.hpp>
#include <syntaxis/induction/induction.hpp>

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace syntaxis::derivation {

using graphs::Graph;
using induction::Grammar;
using induction::Learned;
using shape::Point;
using shape::Similarity;

struct Replay {
  Graph graph;                      // node v is the input's node v
  double worst = 0;                 // the largest distance from an input position, relative to the input's extent
};

inline Replay replay(const Learned& L, const Graph& input) {
  Replay r;
  const std::size_t n = input.size();
  std::vector<Point> at(n, Point{0, 0});
  std::vector<char> placed(n, 0);
  // Down from each root: an instance's frame places its rule's nodes, which are its children in canonical order.
  std::vector<std::pair<std::size_t, Similarity>> stack;
  for (std::size_t root : L.roots) stack.push_back({root, L.instances[root].frame});
  while (!stack.empty()) {
    const auto [id, frame] = stack.back();
    stack.pop_back();
    const auto& in = L.instances[id];
    if (in.terminal != induction::none && in.children.empty()) {
      at[in.terminal] = in.symbol == 0 && in.parent == induction::none ? in.at : at[in.terminal];
      if (in.parent == induction::none) placed[in.terminal] = 1;
      continue;
    }
    const auto& rule = L.grammar.rule(in.symbol);
    for (std::size_t k = 0; k < in.children.size(); k++) {
      const std::size_t c = in.children[k];
      const Point p = frame(rule.nodes[k].at);
      const auto& child = L.instances[c];
      if (child.children.empty()) {
        at[child.terminal] = p;
        placed[child.terminal] = 1;
      } else {
        stack.push_back({c, rule.nodes[k].frame.then(frame)});
      }
    }
  }
  double lo_x = 1e300, hi_x = -1e300, lo_y = 1e300, hi_y = -1e300;
  for (const auto& nd : input.nodes) { lo_x = std::min(lo_x, nd.at.real()); hi_x = std::max(hi_x, nd.at.real()); lo_y = std::min(lo_y, nd.at.imag()); hi_y = std::max(hi_y, nd.at.imag()); }
  const double extent = std::max(1e-300, std::max(hi_x - lo_x, hi_y - lo_y));
  for (std::size_t v = 0; v < n; v++) {
    if (!placed[v]) at[v] = input.nodes[v].at;
    r.graph.add(at[v]);
    r.worst = std::max(r.worst, std::abs(at[v] - input.nodes[v].at) / extent);
  }
  for (const auto& in : L.instances)
    for (std::size_t e : in.swallowed) r.graph.connect(L.input_edges[e].first, L.input_edges[e].second);
  for (const auto& list : L.remaining)
    for (std::size_t e : list) r.graph.connect(L.input_edges[e].first, L.input_edges[e].second);
  return r;
}

// Every nonterminal of `start` replaced until only terminals remain (at most max_nodes nodes are made).
inline Graph expand(const Grammar& G, const Graph& start, std::size_t max_nodes = 200000) {
  Graph g = start;
  std::vector<char> gone(g.size(), 0);
  std::vector<std::size_t> todo;
  for (std::size_t v = 0; v < g.size(); v++) if (g.nodes[v].symbol > 0) todo.push_back(v);
  while (!todo.empty()) {
    const std::size_t x = todo.back();
    todo.pop_back();
    const auto& rule = G.rule(g.nodes[x].symbol);
    const Similarity frame = g.nodes[x].frame;
    if (g.size() + rule.nodes.size() > max_nodes) throw std::length_error("syntaxis::derivation::expand: more nodes than allowed");
    std::vector<std::size_t> made;
    for (const auto& rn : rule.nodes) {
      const std::size_t v = g.add(frame(rn.at), rn.symbol, rn.symbol ? rn.frame.then(frame) : Similarity{});
      gone.push_back(0);
      made.push_back(v);
      if (rn.symbol > 0) todo.push_back(v);
    }
    for (const auto& [a, b] : rule.edges) g.connect(made[a], made[b]);
    // The relaxed embedding: each neighbour of x joins the pattern's node nearest to it.
    const std::vector<std::size_t> around(g.adjacent[x].begin(), g.adjacent[x].end());
    for (std::size_t y : around) {
      std::size_t best = made[0];
      for (std::size_t v : made) if (std::abs(g.nodes[v].at - g.nodes[y].at) < std::abs(g.nodes[best].at - g.nodes[y].at)) best = v;
      g.adjacent[y].erase(x);
      g.connect(best, y);
    }
    g.adjacent[x].clear();
    gone[x] = 1;
  }
  // The terminal graph, without the replaced nodes.
  Graph out;
  std::vector<std::size_t> idx(g.size(), induction::none);
  for (std::size_t v = 0; v < g.size(); v++) if (!gone[v]) idx[v] = out.add(g.nodes[v].at, g.nodes[v].symbol, g.nodes[v].frame);
  for (const auto& [u, v] : g.edges()) if (!gone[u] && !gone[v]) out.connect(idx[u], idx[v]);
  return out;
}

}  // namespace syntaxis::derivation
