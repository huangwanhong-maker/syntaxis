// The syntax analysis of a geometric graph by a given grammar: which parts of it the grammar's rules derive (docs/theory.md, §5).
//
// Parsing runs the grammar backwards. The rules are taken in the order they were made (each rule's pattern holds only terminals and
// the symbols of earlier rules), and for each, the connected node sets of its pattern's size whose induced subgraph has its pattern's
// canonical key — the same shape, up to similarity and within the grammar's tolerance — are found, chosen node-disjoint, and
// contracted to a node of the rule's symbol, exactly as learning contracted them (syntaxis::induction::contract). Passes over the
// rules repeat until no rule applies. Every reduction is a step, kept for replay; what is reduced forms a parse forest whose leaves
// are the input's nodes; the input is *derived* when what remains is the grammar's axiom (the same key), and *covered* in the share of
// its nodes that some rule absorbed.
//
// Limits: a geometric graph may parse in more than one way (an ambiguous grammar), and the greedy choice of node-disjoint
// occurrences finds one parse, not all of them, nor the best. Finding occurrences enumerates connected node sets up to the pattern's
// size, capped by `max_sets` per rule and pass.

#pragma once

#include <syntaxis/graphs/graphs.hpp>
#include <syntaxis/induction/induction.hpp>

#include <algorithm>
#include <cstddef>
#include <set>
#include <string>
#include <vector>

namespace syntaxis::parsing {

using graphs::Graph;
using induction::Grammar;
using induction::Instance;
using induction::Step;

struct Options {
  std::size_t max_passes = 8;
  std::size_t max_sets = 400000;  // connected node sets examined per rule and pass
};

struct Parse {
  std::vector<Instance> instances;     // the input's nodes first (as terminals), then every node the parse made
  std::vector<Step> steps;
  std::vector<std::size_t> roots;      // what remains, as instances
  Graph remaining;
  std::vector<std::size_t> uses;       // per rule, its applications
  double covered = 0;                  // the share of input nodes inside some reduction
  bool derived = false;                // what remains is the axiom
};



inline Parse parse(const Grammar& G, const Graph& input, const Options& o = {}) {
  Parse P;
  P.uses.assign(G.rules.size(), 0);
  induction::Working w;
  w.graph = input;
  for (std::size_t v = 0; v < input.size(); v++) {
    Instance t;
    t.symbol = input.nodes[v].symbol, t.at = input.nodes[v].at, t.frame = input.nodes[v].frame, t.terminal = v;
    P.instances.push_back(t);
    w.instance_of.push_back(v);
  }
  const auto edges = input.edges();
  for (std::size_t e = 0; e < edges.size(); e++) w.carries[edges[e]].push_back(e);
  for (std::size_t pass = 0; pass < o.max_passes; pass++) {
    bool any = false;
    for (const auto& rule : G.rules) {
      std::multiset<int> symbols;
      for (const auto& rn : rule.nodes) symbols.insert(rn.symbol);
      const auto match = induction::occurrences(w.graph, rule.key, rule.nodes.size(), symbols, G.tolerance, o.max_sets);
      if (match.empty()) continue;
      const auto pick = induction::disjoint(match, w.graph.size());
      std::vector<std::vector<std::size_t>> chosen;
      for (std::size_t i : pick) chosen.push_back(match[i]);
      induction::Rule copy = rule;  // contract fills nothing more: the rule is given
      Step step;
      step.symbol = rule.symbol, step.order = rule.nodes.size(), step.before = w.graph, step.instance_of = w.instance_of;
      auto [ordered, made] = induction::contract(w, chosen, rule.symbol, copy, P.instances, G.tolerance, pass + 1);
      step.occurrences = std::move(ordered), step.made = std::move(made);
      P.uses[static_cast<std::size_t>(rule.symbol - 1)] += chosen.size();
      P.steps.push_back(std::move(step));
      any = true;
    }
    if (!any) break;
  }
  P.remaining = w.graph;
  P.roots = w.instance_of;
  std::size_t absorbed = 0;
  for (std::size_t v = 0; v < input.size(); v++) absorbed += P.instances[v].parent != induction::none;
  P.covered = input.size() ? static_cast<double>(absorbed) / static_cast<double>(input.size()) : 0;
  if (w.graph.size() == G.axiom.size() && w.graph.size() >= 2) {
    std::vector<std::size_t> all(w.graph.size());
    for (std::size_t v = 0; v < all.size(); v++) all[v] = v;
    std::vector<std::size_t> ax(G.axiom.size());
    for (std::size_t v = 0; v < ax.size(); v++) ax[v] = v;
    P.derived = graphs::key_of(w.graph, all, G.tolerance).key == graphs::key_of(G.axiom, ax, G.tolerance).key;
  } else if (w.graph.size() == 1 && G.axiom.size() == 1) {
    P.derived = w.graph.nodes[0].symbol == G.axiom.nodes[0].symbol;
  }
  return P;
}

}  // namespace syntaxis::parsing
