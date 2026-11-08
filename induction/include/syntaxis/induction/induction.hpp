// Learning a geometric graph grammar from one graph: grammar induction as compression (docs/theory.md, §3).
//
// After Fiser, Benes, Garcia Galicia, Abdul-Massih, Aliaga and Krs (2016), "Learning geometric graph grammars", SCCG '16, §4–6
// (DOI 10.1145/2948628.2948635), written here from that statement. At each level:
// - **Detection.** Subgraphs are grown from every edge one node at a time, through a neighbour, and the induced subgraphs of each
//   size are grouped by their canonical key (syntaxis::shape): a group is an *isogroup*, a pattern and its occurrences, the same up
//   to similarity and within the tolerance. Only the most frequent groups are grown further (the paper's frequency-based pruning).
// - **Selection.** For each group of at least three nodes, a greedy choice of node-disjoint occurrences; the group chosen is the one
//   whose replacement saves most of the description, nodes and edges counted (m occurrences of k nodes and e edges save m(k − 1 + e)
//   and cost a rule of k + e), a description-length criterion; it prefers a triangle to a path of the same nodes, as the paper's
//   convex-hull heuristic prefers bifurcations to chains.
//   Patterns of two nodes are grown from but never chosen: the similarities act transitively on pairs of distinct points, so every
//   edge is similar to every other and a two-node pattern has no shape; shape begins with three points, whose classes form Kendall's
//   two-dimensional shape space (Kendall 1984, Bull. London Math. Soc. 16).
// - **Encoding.** Each chosen occurrence is replaced by one node of a new nonterminal symbol, at its centroid, carrying the
//   similarity from the rule's canonical frame to the occurrence; edges leaving the occurrence now leave the new node. The rule is
//   the pattern in its canonical frame: its nodes (with their symbols, and for nonterminals their frames relative to it) and edges.
// The levels repeat on the smaller graph until nothing repeats; what remains is the axiom.
//
// The grammar is context-free and deterministic, and since each nonterminal has one rule, its language is the input alone: a
// straight-line grammar, the graph form of Sequitur's and Re-Pair's compressions of strings (Nevill-Manning and Witten 1997; Larsson
// and Moffat 2000). The derivation records, for every node ever made, its children in the rule's canonical order, and the terminal
// edges each contraction swallowed, so that the input is rebuilt exactly (syntaxis::derivation::replay).
//
// Limits: the search is greedy, as the paper's is; the grammar found is a good compression, not the shortest.

#pragma once

#include <syntaxis/graphs/graphs.hpp>
#include <syntaxis/shape/shape.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace syntaxis::induction {

using graphs::Graph;
using shape::Point;
using shape::Similarity;

constexpr std::size_t none = std::numeric_limits<std::size_t>::max();

struct RuleNode {
  Point at;            // in the rule's canonical frame
  int symbol = 0;
  Similarity frame;    // a nonterminal's frame relative to the rule's canonical frame
};

struct Rule {
  int symbol = 0;      // the nonterminal it rewrites (rules[symbol − 1])
  std::vector<RuleNode> nodes;
  std::vector<std::pair<std::size_t, std::size_t>> edges;
  std::string key;     // the pattern's canonical key
  std::size_t uses = 0, level = 0;
};

struct Grammar {
  std::vector<Rule> rules;
  Graph axiom;         // what remains when nothing more repeats
  double tolerance = 0;
  const Rule& rule(int symbol) const { return rules.at(static_cast<std::size_t>(symbol - 1)); }
};

// A node of the derivation: a terminal of the input, or a nonterminal made by a contraction.
struct Instance {
  int symbol = 0;
  Point at;
  Similarity frame;
  std::size_t terminal = none;          // the input's node, for a terminal
  std::size_t parent = none;
  std::vector<std::size_t> children;    // in the rule's canonical order
  std::vector<std::size_t> swallowed;   // the input's edges (indices) that became internal when it was made
  std::size_t level = 0;
};

// One level of the process, for replaying it: the graph before, which of its nodes were replaced together, and by what.
struct Step {
  int symbol = 0;
  std::size_t order = 0;
  std::vector<std::vector<std::size_t>> occurrences;  // nodes of `before`, in canonical order
  Graph before;
  std::vector<std::size_t> instance_of;               // instance of each node of `before`
  std::vector<std::size_t> made;                      // the instances made, one per occurrence
};

struct Options {
  double tolerance = 0.06;          // relative to a pattern's size
  std::size_t max_order = 6;        // nodes in a pattern
  std::size_t max_levels = 40;
  double keep = 0.3;                // the share of groups grown further at each order
  std::size_t max_groups = 48;      // and at most this many
  std::size_t max_occurrences = 40000;  // per order, all groups together
  std::size_t min_count = 2;        // a pattern must repeat
  std::size_t min_order = 3;        // the least pattern chosen: similarity is transitive on pairs of points, so an edge has no shape
};

struct Learned {
  Grammar grammar;
  std::vector<Instance> instances;
  std::vector<std::size_t> roots;                         // the axiom's nodes, as instances
  std::vector<std::pair<std::size_t, std::size_t>> input_edges;  // the input's edges, by which `swallowed` and `remaining` count
  std::vector<std::vector<std::size_t>> remaining;        // for each axiom edge (in axiom.edges() order), the input edges it carries
  std::vector<Step> steps;
};

// ---- shared with the syntax analyser ------------------------------------------------------------------------------------

// The state of a graph under contraction: each node's instance, and each edge's input edges.
struct Working {
  Graph graph;
  std::vector<std::size_t> instance_of;
  std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> carries;  // (u < v) → input edges
};

// Greedy node-disjoint occurrences, in the order given.
inline std::vector<std::size_t> disjoint(const std::vector<std::vector<std::size_t>>& occurrences, std::size_t nodes) {
  std::vector<char> used(nodes, 0);
  std::vector<std::size_t> out;
  for (std::size_t k = 0; k < occurrences.size(); k++) {
    bool free = true;
    for (std::size_t v : occurrences[k]) if (used[v]) { free = false; break; }
    if (!free) continue;
    for (std::size_t v : occurrences[k]) used[v] = 1;
    out.push_back(k);
  }
  return out;
}

// Every connected node set of size k, of nodes whose symbols fit the multiset given, whose induced subgraph has the key given: the
// occurrences of a pattern, exhaustively (grown from each node; capped at `cap` sets examined), in lexicographic order.
inline std::vector<std::vector<std::size_t>> occurrences(const Graph& g, const std::string& key, std::size_t k, const std::multiset<int>& symbols, double tolerance,
                                                         std::size_t cap = 400000) {
  std::set<std::vector<std::size_t>> seen;
  std::vector<std::vector<std::size_t>> frontier;
  for (std::size_t v = 0; v < g.size(); v++)
    if (symbols.count(g.nodes[v].symbol)) frontier.push_back({v});
  for (std::size_t size = 1; size < k && !frontier.empty(); size++) {
    std::vector<std::vector<std::size_t>> next;
    for (const auto& s : frontier) {
      std::set<std::size_t> around;
      for (std::size_t v : s) for (std::size_t q : g.adjacent[v]) if (!std::binary_search(s.begin(), s.end(), q)) around.insert(q);
      for (std::size_t q : around) {
        if (!symbols.count(g.nodes[q].symbol)) continue;
        std::vector<std::size_t> t = s;
        t.insert(std::upper_bound(t.begin(), t.end(), q), q);
        std::multiset<int> have;
        for (std::size_t v : t) have.insert(g.nodes[v].symbol);
        bool fits = true;
        for (int sy : have) if (have.count(sy) > symbols.count(sy)) { fits = false; break; }
        if (!fits || !seen.insert(t).second) continue;
        next.push_back(std::move(t));
        if (seen.size() > cap) break;
      }
      if (seen.size() > cap) break;
    }
    frontier = std::move(next);
  }
  std::vector<std::vector<std::size_t>> out;
  for (auto& s : frontier)
    if (s.size() == k && graphs::key_of(g, s, tolerance).key == key) out.push_back(std::move(s));
  std::sort(out.begin(), out.end());
  return out;
}

// Replaces each occurrence (node sets of w.graph) by one node of `symbol`, recording instances. Returns the occurrences in canonical
// order and the instances made; `rule` is filled from the first occurrence when it is empty.
inline std::pair<std::vector<std::vector<std::size_t>>, std::vector<std::size_t>> contract(Working& w, const std::vector<std::vector<std::size_t>>& occurrences, int symbol,
                                                                                         Rule& rule, std::vector<Instance>& instances, double tolerance, std::size_t level) {
  const Graph& g = w.graph;
  std::vector<std::size_t> replaced_by(g.size(), none);
  std::vector<std::vector<std::size_t>> ordered;
  std::vector<std::size_t> made;
  Graph next;
  std::vector<std::size_t> next_instance;
  // Nodes that stay keep their order, first; the new nodes follow.
  std::vector<std::size_t> new_index(g.size(), none);
  std::vector<char> inside(g.size(), 0);
  for (const auto& o : occurrences) for (std::size_t v : o) inside[v] = 1;
  for (std::size_t v = 0; v < g.size(); v++)
    if (!inside[v]) {
      new_index[v] = next.add(g.nodes[v].at, g.nodes[v].symbol, g.nodes[v].frame);
      next_instance.push_back(w.instance_of[v]);
    }
  for (const auto& o : occurrences) {
    const shape::Canonical c = graphs::key_of(g, o, tolerance);
    std::vector<std::size_t> canon(o.size());
    for (std::size_t k = 0; k < o.size(); k++) canon[k] = o[c.order[k]];
    const Similarity to_plane = c.frame.inverse();
    if (rule.nodes.empty()) {
      rule.symbol = symbol, rule.key = c.key, rule.level = level;
      for (std::size_t k = 0; k < canon.size(); k++) {
        const auto& n = g.nodes[canon[k]];
        rule.nodes.push_back({c.coordinates[k], n.symbol, n.symbol ? n.frame.then(c.frame) : Similarity{}});
      }
      std::vector<std::size_t> pos(g.size(), none);
      for (std::size_t k = 0; k < canon.size(); k++) pos[canon[k]] = k;
      for (std::size_t a : canon)
        for (std::size_t b : g.adjacent[a])
          if (pos[b] != none && pos[a] < pos[b]) rule.edges.push_back({pos[a], pos[b]});
      std::sort(rule.edges.begin(), rule.edges.end());
    }
    rule.uses++;
    Point centre = 0;
    for (std::size_t v : canon) centre += g.nodes[v].at;
    centre /= static_cast<double>(canon.size());
    Instance in;
    in.symbol = symbol, in.at = centre, in.frame = to_plane, in.level = level;
    for (std::size_t v : canon) in.children.push_back(w.instance_of[v]);
    const std::size_t id = instances.size();
    for (std::size_t v : canon) instances[w.instance_of[v]].parent = id;
    // The edges inside the occurrence are swallowed.
    for (std::size_t i = 0; i < canon.size(); i++)
      for (std::size_t j = i + 1; j < canon.size(); j++) {
        const auto e = std::minmax(canon[i], canon[j]);
        const auto it = w.carries.find(e);
        if (it != w.carries.end()) in.swallowed.insert(in.swallowed.end(), it->second.begin(), it->second.end());
      }
    instances.push_back(std::move(in));
    const std::size_t node = next.add(centre, symbol, to_plane);
    next_instance.push_back(id);
    for (std::size_t v : canon) new_index[v] = node;
    ordered.push_back(canon);
    made.push_back(id);
  }
  std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> carries;
  for (const auto& [e, list] : w.carries) {
    const std::size_t a = new_index[e.first], b = new_index[e.second];
    if (a == b) continue;
    next.connect(a, b);
    auto& to = carries[std::minmax(a, b)];
    to.insert(to.end(), list.begin(), list.end());
  }
  w.graph = std::move(next);
  w.instance_of = std::move(next_instance);
  w.carries = std::move(carries);
  return {ordered, made};
}

// ---- learning ---------------------------------------------------------------------------------------------------------

namespace detail {

struct Group {
  std::size_t order = 0;
  std::vector<std::vector<std::size_t>> occurrences;  // sorted node sets
};

// The isogroups of every order up to max_order, grown from the edges.
inline std::map<std::string, Group> detect(const Graph& g, const Options& o) {
  std::map<std::string, Group> all;
  std::map<std::string, Group> current;
  for (const auto& [u, v] : g.edges()) {
    std::vector<std::size_t> s{u, v};
    auto& grp = current[graphs::key_of(g, s, o.tolerance).key];
    grp.order = 2;
    grp.occurrences.push_back(s);
  }
  for (std::size_t order = 2;; order++) {
    // Keep the most frequent groups of this order.
    std::vector<std::pair<std::size_t, std::string>> by_count;
    for (const auto& [k, grp] : current) if (grp.occurrences.size() >= o.min_count) by_count.push_back({grp.occurrences.size(), k});
    std::sort(by_count.begin(), by_count.end(), [](const auto& a, const auto& b) { return a.first != b.first ? a.first > b.first : a.second < b.second; });
    const std::size_t keep = std::min(o.max_groups, std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(o.keep * by_count.size()))));
    if (by_count.size() > keep) by_count.resize(keep);
    std::map<std::string, Group> kept;
    for (const auto& [count, k] : by_count) kept[k] = std::move(current[k]);
    for (auto& [k, grp] : kept) all[k] = grp;
    if (order >= o.max_order || kept.empty()) break;
    // Grow each occurrence by one neighbouring node.
    std::map<std::string, Group> grown;
    std::set<std::vector<std::size_t>> seen;
    std::size_t budget = o.max_occurrences;
    for (const auto& [k, grp] : kept)
      for (const auto& occ : grp.occurrences) {
        std::set<std::size_t> around;
        for (std::size_t v : occ) for (std::size_t q : g.adjacent[v]) if (!std::binary_search(occ.begin(), occ.end(), q)) around.insert(q);
        for (std::size_t q : around) {
          std::vector<std::size_t> s = occ;
          s.insert(std::upper_bound(s.begin(), s.end(), q), q);
          if (!seen.insert(s).second) continue;
          if (budget == 0) break;
          budget--;
          auto& to = grown[graphs::key_of(g, s, o.tolerance).key];
          to.order = order + 1;
          to.occurrences.push_back(std::move(s));
        }
      }
    current = std::move(grown);
  }
  return all;
}

}  // namespace detail

inline Learned learn(const Graph& input, const Options& o = {}) {
  Learned L;
  L.grammar.tolerance = o.tolerance;
  Working w;
  w.graph = input;
  for (std::size_t v = 0; v < input.size(); v++) {
    Instance t;
    t.symbol = input.nodes[v].symbol, t.at = input.nodes[v].at, t.frame = input.nodes[v].frame, t.terminal = v;
    L.instances.push_back(t);
    w.instance_of.push_back(v);
  }
  L.input_edges = input.edges();
  for (std::size_t e = 0; e < L.input_edges.size(); e++) w.carries[L.input_edges[e]].push_back(e);
  for (std::size_t level = 1; level <= o.max_levels && w.graph.size() > 1; level++) {
    const auto groups = detail::detect(w.graph, o);
    // The group whose replacement saves most.
    double best = 0;
    std::string best_key;
    std::vector<std::vector<std::size_t>> best_occ;
    for (const auto& [k, grp] : groups) {
      if (grp.order < o.min_order) continue;
      auto occ = grp.occurrences;
      std::sort(occ.begin(), occ.end());
      const auto pick = disjoint(occ, w.graph.size());
      if (pick.size() < o.min_count) continue;
      // The description saved: m occurrences of k nodes and e edges save m(k − 1 + e) and cost a rule of k + e.
      const double e = static_cast<double>(graphs::induced_edges(w.graph, occ[pick[0]]).size());
      const double gain = static_cast<double>(pick.size()) * (static_cast<double>(grp.order - 1) + e) - (static_cast<double>(grp.order) + e);
      if (gain > best || (gain == best && gain > 0 && grp.order > (best_occ.empty() ? 0 : best_occ[0].size()))) {
        best = gain, best_key = k;
        best_occ.clear();
        for (std::size_t i : pick) best_occ.push_back(occ[i]);
      }
    }
    if (best <= 0 || best_occ.empty()) break;
    // Every occurrence of the pattern chosen, found as the syntax analyser finds them, so that parsing the input repeats learning.
    {
      std::multiset<int> symbols;
      for (std::size_t v : best_occ[0]) symbols.insert(w.graph.nodes[v].symbol);
      const auto all = occurrences(w.graph, best_key, best_occ[0].size(), symbols, o.tolerance);
      const auto pick = disjoint(all, w.graph.size());
      if (pick.size() >= best_occ.size()) {
        best_occ.clear();
        for (std::size_t i : pick) best_occ.push_back(all[i]);
      }
    }
    const int symbol = static_cast<int>(L.grammar.rules.size()) + 1;
    Rule rule;
    Step step;
    step.symbol = symbol, step.order = best_occ[0].size(), step.before = w.graph, step.instance_of = w.instance_of;
    auto [ordered, made] = contract(w, best_occ, symbol, rule, L.instances, o.tolerance, level);
    step.occurrences = std::move(ordered), step.made = std::move(made);
    L.grammar.rules.push_back(std::move(rule));
    L.steps.push_back(std::move(step));
  }
  L.grammar.axiom = w.graph;
  L.roots = w.instance_of;
  for (const auto& e : w.graph.edges()) L.remaining.push_back(w.carries[e]);
  return L;
}

}  // namespace syntaxis::induction
