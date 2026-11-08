// Checks of syntaxis/induction: grammars learned by compression.
//
// - The Sierpinski graph of depth 4 (123 nodes) compresses: its rules are used at least twice each, and the description (rules and
//   axiom) is far smaller than the graph. Its sub-triangles share corners, which node replacement cannot express, so its axiom keeps a
//   few nodes; a self-similar graph whose copies are joined by edges reduces, with patterns of three nodes, to its top triangle, by
//   a rule for each depth. A rotated, scaled and moved copy learns rules of the same keys, in
//   the same order (the grammar of a shape does not depend on where it lies).
// - A grid of 8 × 8 compresses too; a graph of random points scarcely does (few rules, a large axiom).
// - The derivation is a forest whose leaves are the input's nodes, each with one parent or none.

#include <syntaxis/induction/induction.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace syntaxis;

namespace {
std::size_t description(const induction::Learned& L) {
  std::size_t d = L.grammar.axiom.size() + L.grammar.axiom.edge_count();
  for (const auto& r : L.grammar.rules) d += r.nodes.size() + r.edges.size();
  return d;
}
}  // namespace

int main() {
  calx::verification::Suite s;

  const graphs::Graph S = graphs::sierpinski(4);
  const induction::Learned L = induction::learn(S);
  std::size_t least = 1000;
  for (const auto& r : L.grammar.rules) least = std::min(least, r.uses);
  s.require("sierpinski(4): learns rules, each used at least twice", !L.grammar.rules.empty() && least >= 2,
            std::to_string(L.grammar.rules.size()) + " rules, least used " + std::to_string(least));
  // The gasket's sub-triangles share their corners, and node replacement contracts node-disjoint occurrences only: it cannot tile
  // the gasket, whose natural grammar is a hyperedge-replacement grammar (planned). What remains is small, not a single node.
  s.require("sierpinski(4): the axiom keeps under a fifth of the nodes (shared corners stop node replacement short of one node)", L.grammar.axiom.size() * 5 < S.size(),
            std::to_string(L.grammar.axiom.size()) + " nodes");
  // With patterns of three nodes the depths appear as levels: a rule for each, used 81, 27, 9 and 3 times; the top triangle occurs
  // once, and a pattern that occurs once saves nothing (m(k − 1) − k < 0), so it is the axiom. With larger patterns the greedy choice
  // prefers pairs of triangles and their bridge, which save more at the first level, and the levels no longer follow the depths.
  const graphs::Graph B = graphs::bridged_triangles(4);
  induction::Options three;
  three.max_order = 3;
  const induction::Learned LB = induction::learn(B, three);
  bool counts = LB.grammar.rules.size() == 4;
  const std::size_t want[4] = {81, 27, 9, 3};
  for (std::size_t k = 0; counts && k < 4; k++) counts = LB.grammar.rules[k].uses == want[k];
  s.require("bridged triangles (depth 4, 243 nodes), patterns of three: a rule per depth used 81, 27, 9, 3 times; the top triangle the axiom", counts && LB.grammar.axiom.size() == 3,
            std::to_string(LB.grammar.rules.size()) + " rules, axiom of " + std::to_string(LB.grammar.axiom.size()));
  const std::size_t whole = S.size() + S.edge_count();
  s.require("sierpinski(4): the description is under a third of the graph's", description(L) * 3 < whole, std::to_string(description(L)) + " of " + std::to_string(whole));
  const graphs::Graph T = graphs::transformed(S, {std::polar(7.0, 2.2), {100, -40}});
  const induction::Learned LT = induction::learn(T);
  bool same = LT.grammar.rules.size() == L.grammar.rules.size();
  for (std::size_t k = 0; same && k < L.grammar.rules.size(); k++) same = LT.grammar.rules[k].key == L.grammar.rules[k].key;
  s.require("sierpinski(4): a rotated, scaled and moved copy learns the same rules, in order", same);

  const graphs::Graph G = graphs::grid(8, 8);
  const induction::Learned LG = induction::learn(G);
  s.require("grid(8, 8): compresses (description under half the graph's)", description(LG) * 2 < G.size() + G.edge_count(),
            std::to_string(description(LG)) + " of " + std::to_string(G.size() + G.edge_count()));

  calx::random::Stream st(5, 2, 0, 0);
  graphs::Graph R;
  for (int i = 0; i < 80; i++) R.add({st.uniform(0, 10), st.uniform(0, 10)});
  for (std::size_t i = 0; i < 80; i++) {
    // Each node to its two nearest.
    std::vector<std::pair<double, std::size_t>> d;
    for (std::size_t j = 0; j < 80; j++) if (j != i) d.push_back({std::abs(R.nodes[i].at - R.nodes[j].at), j});
    std::sort(d.begin(), d.end());
    R.connect(i, d[0].second), R.connect(i, d[1].second);
  }
  const induction::Learned LR = induction::learn(R);
  s.require("random points: scarcely compress (the axiom keeps most of the graph)", LR.grammar.axiom.size() * 2 > R.size(),
            std::to_string(LR.grammar.rules.size()) + " rules, axiom of " + std::to_string(LR.grammar.axiom.size()));

  bool forest = true;
  std::size_t leaves = 0;
  for (std::size_t i = 0; i < L.instances.size(); i++) {
    const auto& in = L.instances[i];
    if (in.children.empty()) leaves++;
    for (std::size_t c : in.children) forest = forest && L.instances[c].parent == i;
  }
  s.require("derivation: a forest whose leaves are the input's 123 nodes", forest && leaves == S.size(), std::to_string(leaves));
  return s.done();
}
