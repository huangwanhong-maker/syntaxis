// Probabilistic and context-sensitive geometric graph grammars, from a learned straight-line grammar (docs/theory.md, §4).
//
// A straight-line grammar generates only its input. To have a language, its nonterminals are merged into *categories*: rules whose
// patterns have the same shape at a coarser `merge` tolerance, their symbols set aside, are alternatives of one category. A
// category is a nonterminal of the new grammar, and its alternatives its productions (the merging of Stolcke and Omohundro 1994,
// "Inducing probabilistic grammars by Bayesian model merging", here by geometric equivalence). Scale invariance makes the grammar
// recursive where the data are self-similar: the Sierpinski triangle's rules at every depth share one shape, so its category
// rewrites into copies of itself and generates every depth.
//
// - **PCFG.** Each production's probability is its share of its category's uses in the derivation, the maximum-likelihood estimate
//   from a parsed sample (Chi and Geman 1998, "Estimation of probabilistic context-free grammars", Comput. Linguist. 24(2): such
//   estimates make the grammar consistent, its derivations finite with probability one).
// - **PCSG.** A production's probability given its context, the category of the node's parent, estimated with additive smoothing
//   towards the PCFG (the parent annotation of Johnson 1998, "PCFG models of linguistic tree representations", Comput. Linguist.
//   24(4)); and a context condition on the geometry: a production whose placed edges would cross an edge already drawn is refused,
//   and another drawn. A context-free grammar cannot state either condition; both depend on what surrounds the node rewritten.
// - **sample.** A start graph's nonterminals rewritten by productions drawn from the grammar, placed by the nodes' frames, joined
//   to their neighbours by the relaxed embedding (syntaxis::derivation::expand), until only terminals remain or a depth is reached,
//   beyond which the production with the fewest nonterminals is taken.
// - **statistics.** What a generated graph is like: its degree distribution, the spread of its edges' lengths, its crossings.

#pragma once

#include <syntaxis/graphs/graphs.hpp>
#include <syntaxis/induction/induction.hpp>

#include <calx/random/stream.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace syntaxis::stochastic {

using graphs::Graph;
using induction::Grammar;
using induction::Learned;
using shape::Point;
using shape::Similarity;

struct Category {
  std::string shape;                 // the pattern's key without symbols, at the merge tolerance
  std::vector<int> rules;            // its productions: the learned rules merged into it
  std::vector<double> probability;   // their probabilities (PCFG)
  std::vector<std::size_t> count;    // their uses in the derivation
};

struct Model {
  Grammar grammar;
  std::vector<std::size_t> category_of;  // per rule (index symbol − 1)
  std::vector<Category> categories;
  // Context counts: (category, the parent's category or none) → per production of the category, its uses.
  std::map<std::pair<std::size_t, std::size_t>, std::vector<double>> context;
  double smoothing = 1.0;
};

constexpr std::size_t no_parent = static_cast<std::size_t>(-1);

inline std::string shape_key(const induction::Rule& r, double tolerance) {
  std::vector<Point> p;
  std::vector<int> label(r.nodes.size(), 0);
  for (const auto& n : r.nodes) p.push_back(n.at);
  return shape::canonical(p, label, r.edges, tolerance).key;
}

inline Model model(const Learned& L, double merge, double smoothing = 1.0) {
  Model m;
  m.grammar = L.grammar;
  m.smoothing = smoothing;
  std::map<std::string, std::size_t> by_shape;
  for (const auto& r : m.grammar.rules) {
    const std::string k = shape_key(r, merge);
    auto [it, fresh] = by_shape.emplace(k, m.categories.size());
    if (fresh) m.categories.push_back({k, {}, {}, {}});
    m.category_of.push_back(it->second);
    auto& c = m.categories[it->second];
    c.rules.push_back(r.symbol);
    c.count.push_back(0);
  }
  const auto slot = [&](int symbol) {
    const auto& c = m.categories[m.category_of[static_cast<std::size_t>(symbol - 1)]];
    return static_cast<std::size_t>(std::find(c.rules.begin(), c.rules.end(), symbol) - c.rules.begin());
  };
  for (const auto& in : L.instances) {
    if (in.symbol <= 0 || in.children.empty()) continue;
    const std::size_t cat = m.category_of[static_cast<std::size_t>(in.symbol - 1)];
    m.categories[cat].count[slot(in.symbol)]++;
    std::size_t parent = no_parent;
    if (in.parent != induction::none && L.instances[in.parent].symbol > 0) parent = m.category_of[static_cast<std::size_t>(L.instances[in.parent].symbol - 1)];
    auto& ctx = m.context[{cat, parent}];
    if (ctx.empty()) ctx.assign(m.categories[cat].rules.size(), 0.0);
    ctx[slot(in.symbol)] += 1;
  }
  for (auto& c : m.categories) {
    double total = 0;
    for (auto v : c.count) total += static_cast<double>(v);
    for (auto v : c.count) c.probability.push_back(total > 0 ? static_cast<double>(v) / total : 1.0 / static_cast<double>(c.count.size()));
  }
  return m;
}

// P(production | category, parent) with additive smoothing towards the PCFG.
inline std::vector<double> conditional(const Model& m, std::size_t cat, std::size_t parent) {
  const auto& c = m.categories[cat];
  std::vector<double> p = c.probability;
  const auto it = m.context.find({cat, parent});
  if (it == m.context.end()) return p;
  double n = 0;
  for (double v : it->second) n += v;
  for (std::size_t k = 0; k < p.size(); k++) p[k] = (it->second[k] + m.smoothing * c.probability[k]) / (n + m.smoothing);
  return p;
}

enum class Kind { pcfg, pcsg };

struct SampleOptions {
  Kind kind = Kind::pcfg;
  std::size_t max_depth = 12;
  std::size_t max_nodes = 60000;
  std::size_t tries = 6;   // productions drawn before one that crosses nothing is given up (PCSG)
};

struct Sample {
  Graph graph;
  std::size_t rewrites = 0, refused = 0;  // productions applied; productions refused for crossing (PCSG)
  bool truncated = false;                 // the node budget ended the derivation
};

namespace detail {
inline std::size_t draw(const std::vector<double>& p, calx::random::Stream& s) {
  double total = 0;
  for (double v : p) total += v;
  double u = s.uniform() * total;
  for (std::size_t k = 0; k < p.size(); k++) { if (u < p[k]) return k; u -= p[k]; }
  return p.size() - 1;
}
inline std::size_t nonterminals(const induction::Rule& r) {
  std::size_t n = 0;
  for (const auto& x : r.nodes) n += x.symbol > 0;
  return n;
}
}  // namespace detail

inline Sample sample(const Model& m, const Graph& start, calx::random::Stream& s, const SampleOptions& o = {}) {
  Sample out;
  Graph g = start;
  struct Pending { std::size_t node, depth, parent_cat; };
  std::vector<Pending> todo;
  std::vector<char> gone(g.size(), 0);
  for (std::size_t v = 0; v < g.size(); v++) if (g.nodes[v].symbol > 0) todo.push_back({v, 0, no_parent});
  for (std::size_t head = 0; head < todo.size(); head++) {
    const Pending x = todo[head];
    if (g.size() > o.max_nodes) { out.truncated = true; break; }
    const std::size_t cat = m.category_of[static_cast<std::size_t>(g.nodes[x.node].symbol - 1)];
    const auto& c = m.categories[cat];
    std::vector<double> p = o.kind == Kind::pcsg ? conditional(m, cat, x.parent_cat) : c.probability;
    if (x.depth >= o.max_depth) {
      // Beyond the depth, the production with the fewest nonterminals.
      std::size_t best = 0;
      for (std::size_t k = 1; k < c.rules.size(); k++)
        if (detail::nonterminals(m.grammar.rule(c.rules[k])) < detail::nonterminals(m.grammar.rule(c.rules[best]))) best = k;
      std::fill(p.begin(), p.end(), 0.0);
      p[best] = 1;
    }
    const Similarity frame = g.nodes[x.node].frame;
    std::size_t choice = detail::draw(p, s);
    std::vector<Point> placed;
    if (o.kind == Kind::pcsg) {
      // The context condition: a production whose edges cross an edge already drawn is refused.
      const auto edges = g.edges();
      for (std::size_t t = 0; t < o.tries; t++) {
        const auto& r = m.grammar.rule(c.rules[choice]);
        bool crosses = false;
        for (const auto& [a, b] : r.edges) {
          const Point pa = frame(r.nodes[a].at), pb = frame(r.nodes[b].at);
          for (const auto& [u, v] : edges) {
            if (u == x.node || v == x.node || gone[u] || gone[v]) continue;
            if (shape::segments_cross(pa, pb, g.nodes[u].at, g.nodes[v].at)) { crosses = true; break; }
          }
          if (crosses) break;
        }
        if (!crosses) break;
        out.refused++;
        p[choice] = 0;
        double left = 0;
        for (double v : p) left += v;
        if (left <= 0) break;
        choice = detail::draw(p, s);
      }
    }
    const auto& r = m.grammar.rule(c.rules[choice]);
    std::vector<std::size_t> made;
    for (const auto& rn : r.nodes) {
      const std::size_t v = g.add(frame(rn.at), rn.symbol, rn.symbol ? rn.frame.then(frame) : Similarity{});
      gone.push_back(0);
      made.push_back(v);
      if (rn.symbol > 0) todo.push_back({v, x.depth + 1, cat});
    }
    for (const auto& [a, b] : r.edges) g.connect(made[a], made[b]);
    const std::vector<std::size_t> around(g.adjacent[x.node].begin(), g.adjacent[x.node].end());
    for (std::size_t y : around) {
      std::size_t best = made[0];
      for (std::size_t v : made) if (std::abs(g.nodes[v].at - g.nodes[y].at) < std::abs(g.nodes[best].at - g.nodes[y].at)) best = v;
      g.adjacent[y].erase(x.node);
      g.connect(best, y);
    }
    g.adjacent[x.node].clear();
    gone[x.node] = 1;
    out.rewrites++;
  }
  std::vector<std::size_t> idx(g.size(), induction::none);
  for (std::size_t v = 0; v < g.size(); v++) if (!gone[v] && g.nodes[v].symbol == 0) idx[v] = out.graph.add(g.nodes[v].at);
  for (const auto& [u, v] : g.edges()) if (idx[u] != induction::none && idx[v] != induction::none) out.graph.connect(idx[u], idx[v]);
  return out;
}

struct Statistics {
  std::size_t nodes = 0, edges = 0, crossings = 0;
  std::vector<double> degree;     // the share of nodes of degree 0, 1, …, 7, and 8 or more
  double mean_length = 0, cv_length = 0;
};

inline Statistics statistics(const Graph& g, bool count_crossings = true) {
  Statistics st;
  st.nodes = g.size();
  st.edges = g.edge_count();
  st.degree.assign(9, 0.0);
  for (const auto& a : g.adjacent) st.degree[std::min<std::size_t>(a.size(), 8)] += 1;
  for (double& d : st.degree) d /= std::max<double>(1, static_cast<double>(g.size()));
  double s = 0, ss = 0;
  const auto e = g.edges();
  for (const auto& [u, v] : e) { const double l = std::abs(g.nodes[u].at - g.nodes[v].at); s += l; ss += l * l; }
  if (!e.empty()) {
    st.mean_length = s / static_cast<double>(e.size());
    st.cv_length = std::sqrt(std::max(0.0, ss / static_cast<double>(e.size()) - st.mean_length * st.mean_length)) / st.mean_length;
  }
  if (count_crossings) st.crossings = graphs::crossings(g);
  return st;
}

// Half the L1 distance between two degree distributions (total variation).
inline double degree_distance(const Statistics& a, const Statistics& b) {
  double d = 0;
  for (std::size_t k = 0; k < a.degree.size(); k++) d += std::fabs(a.degree[k] - b.degree[k]);
  return d / 2;
}

}  // namespace syntaxis::stochastic
