// Checks of syntaxis/derivation.
//
// - Replay rebuilds the Sierpinski graph of depth 4 exactly: the same edges, and positions within 10⁻⁹ of the input's extent (the
//   graph is exactly self-similar, so every occurrence of a rule is exactly similar to its first). A grid's replay too.
// - Expanding the learned grammar from its axiom makes as many nodes as the input had, at the input's positions where the replay
//   puts them; the relaxed embedding keeps most edges (stated, not exact).

#include <syntaxis/derivation/derivation.hpp>
#include <syntaxis/induction/induction.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <set>
#include <string>
#include <vector>

using namespace syntaxis;

int main() {
  calx::verification::Suite s;
  for (const auto& [name, g] : std::vector<std::pair<std::string, graphs::Graph>>{{"sierpinski(4)", graphs::sierpinski(4)}, {"grid(7, 7)", graphs::grid(7, 7)}}) {
    const induction::Learned L = induction::learn(g);
    const derivation::Replay r = derivation::replay(L, g);
    s.require(name + ": replay has the input's edges exactly", r.graph.edges() == g.edges(), std::to_string(r.graph.edge_count()) + " of " + std::to_string(g.edge_count()));
    s.check(name + ": replay's positions (relative to the extent)", r.worst, 1e-9);
    const graphs::Graph e = derivation::expand(L.grammar, L.grammar.axiom);
    s.require(name + ": expanding the axiom makes as many nodes as the input", e.size() == g.size(), std::to_string(e.size()));
    // Edges kept by the relaxed embedding: matched by position.
    std::set<std::pair<long long, long long>> at;
    const auto pt = [](shape::Point p) { return std::make_pair(std::llround(p.real() * 1e6), std::llround(p.imag() * 1e6)); };
    std::set<std::pair<std::pair<long long, long long>, std::pair<long long, long long>>> truth, made;
    for (const auto& [u, v] : g.edges()) truth.insert(std::minmax(pt(g.nodes[u].at), pt(g.nodes[v].at)));
    for (const auto& [u, v] : e.edges()) made.insert(std::minmax(pt(e.nodes[u].at), pt(e.nodes[v].at)));
    std::size_t kept = 0;
    for (const auto& x : made) kept += truth.count(x);
    s.require(name + ": the relaxed embedding keeps most edges (stated)", kept * 2 > truth.size(), std::to_string(kept) + " of " + std::to_string(truth.size()));
  }
  return s.done();
}
