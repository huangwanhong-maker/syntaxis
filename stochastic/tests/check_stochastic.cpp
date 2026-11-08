// Checks of syntaxis/stochastic: probabilistic and context-sensitive grammars from learned ones.
//
// - The bridged triangles' rules of every depth merge into one category, recursive (a production holds nonterminals of its own
//   category): a straight-line grammar becomes a grammar of every depth, and its samples come in several sizes.
// - Each category's probabilities sum to one and are its productions' shares of uses.
// - Samples from the Sierpinski PCFG lie on the gasket: every edge's length is the base's halved a whole number of times.
// - The context-sensitive grammar refuses productions that would cross what is drawn, on a learned grid perturbed into a city of
//   uneven blocks, and its samples cross no more than the context-free grammar's, on average over twenty seeds.
// - Conditional probabilities given a context sum to one, and fall back to the PCFG for a context never seen.

#include <syntaxis/induction/induction.hpp>
#include <syntaxis/stochastic/stochastic.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <set>
#include <string>
#include <vector>

using namespace syntaxis;

int main() {
  calx::verification::Suite s;
  const graphs::Graph B = graphs::bridged_triangles(4);
  induction::Options three;
  three.max_order = 3;
  const induction::Learned LB = induction::learn(B, three);
  const stochastic::Model mb = stochastic::model(LB, 0.12);
  bool recursive = false;
  for (std::size_t c = 0; c < mb.categories.size(); c++)
    for (int r : mb.categories[c].rules)
      for (const auto& n : mb.grammar.rule(r).nodes)
        if (n.symbol > 0 && mb.category_of[static_cast<std::size_t>(n.symbol - 1)] == c) recursive = true;
  s.require("bridged triangles: the rules of every depth merge into one recursive category", mb.categories.size() == 1 && recursive,
            std::to_string(mb.grammar.rules.size()) + " rules, " + std::to_string(mb.categories.size()) + " categories");
  // Sampled from the recursive grammar, the depths vary: graphs of several sizes.
  {
    calx::random::Stream sb(13, 7, 0, 0);
    std::set<std::size_t> sizes;
    for (int r = 0; r < 30; r++) sizes.insert(stochastic::sample(mb, LB.grammar.axiom, sb).graph.size());
    s.require("bridged triangles: samples of the recursive grammar have several sizes (depths drawn)", sizes.size() >= 3, std::to_string(sizes.size()) + " sizes");
  }
  const graphs::Graph S = graphs::sierpinski(4);
  const induction::Learned L = induction::learn(S);
  const stochastic::Model m = stochastic::model(L, 0.12);
  double worst = 0;
  for (const auto& c : m.categories) {
    double sum = 0, uses = 0;
    for (double p : c.probability) sum += p;
    for (auto n : c.count) uses += static_cast<double>(n);
    worst = std::max(worst, std::fabs(sum - 1));
    for (std::size_t k = 0; k < c.count.size() && uses > 0; k++) worst = std::max(worst, std::fabs(c.probability[k] - static_cast<double>(c.count[k]) / uses));
  }
  s.check("each category's probabilities sum to one and are their productions' shares of uses", worst, 1e-12);

  calx::random::Stream st(7, 3, 0, 0);
  bool lattice = true;
  std::size_t edges = 0;
  for (int r = 0; r < 10; r++) {
    const auto smp = stochastic::sample(m, L.grammar.axiom, st);
    for (const auto& [u, v] : smp.graph.edges()) {
      const double len = std::abs(smp.graph.nodes[u].at - smp.graph.nodes[v].at);
      const double k = -std::log2(len);
      lattice = lattice && std::fabs(k - std::round(k)) < 1e-6;
      edges++;
    }
  }
  s.require("sierpinski PCFG: every sampled edge is the base halved a whole number of times", lattice && edges > 0, std::to_string(edges) + " edges");

  // A grid perturbed into uneven blocks, learned at a loose tolerance.
  graphs::Graph city = graphs::grid(9, 9);
  for (auto& nd : city.nodes) nd.at += shape::Point{0.18 * std::sin(nd.at.imag() * 1.7), 0.18 * std::cos(nd.at.real() * 1.3)};
  induction::Options o;
  o.tolerance = 0.12;
  const induction::Learned LC = induction::learn(city, o);
  const stochastic::Model mc = stochastic::model(LC, 0.3);
  double cf = 0, cs = 0;
  std::size_t refused = 0;
  stochastic::SampleOptions pcfg, pcsg;
  pcsg.kind = stochastic::Kind::pcsg;
  for (int r = 0; r < 20; r++) {
    calx::random::Stream a(100 + r, 4, 0, 0), b(100 + r, 4, 0, 0);
    cf += static_cast<double>(stochastic::statistics(stochastic::sample(mc, LC.grammar.axiom, a, pcfg).graph).crossings);
    const auto x = stochastic::sample(mc, LC.grammar.axiom, b, pcsg);
    cs += static_cast<double>(stochastic::statistics(x.graph).crossings);
    refused += x.refused;
  }
  s.require("city: the context-sensitive grammar crosses no more than the context-free one, on average over twenty seeds", cs <= cf,
            "crossings " + std::to_string(cs / 20) + " against " + std::to_string(cf / 20) + "; refused " + std::to_string(refused));
  double cond = 0;
  for (const auto& [ctx, counts] : mc.context) {
    const auto p = stochastic::conditional(mc, ctx.first, ctx.second);
    double sum = 0;
    for (double v : p) sum += v;
    cond = std::max(cond, std::fabs(sum - 1));
  }
  s.check("PCSG: probabilities given each context sum to one", cond, 1e-12);
  const auto unseen = stochastic::conditional(mc, 0, 987654);
  s.require("PCSG: an unseen context falls back to the PCFG", unseen == mc.categories[0].probability);
  return s.done();
}
