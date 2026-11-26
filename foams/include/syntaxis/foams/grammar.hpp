// The stochastic grammar of a foam's stars (docs/theory.md, §7.3–7.5): rules weighed by what they consume.
//
// A grammar learned from a set of events counts their rules and, for each, the match and the context it was applied in. The
// probability of a rule r is read at three levels, each a hierarchical Dirichlet backoff to the one before (MacKay and Peto 1995,
// "A hierarchical Dirichlet language model", Nat. Lang. Eng. 1, 289–308):
//
//   level 0 (free)       P₀(r) = n(r) / (N + γ) for a rule seen, γ / (N + γ) for one not seen ("new");
//   level 1 (match)      P₁(r | m) = (n(m, r) + α₁ P₀(r)) / (n(m) + α₁);
//   level 2 (context)    P₂(r | c) = (n(c, r) + α₂ P₁(r | m)) / (n(c) + α₂), m the match of the context c.
//
// Each is a distribution over the rules seen and the one outcome "a new rule". The escape γ and the concentrations α₁, α₂ are chosen
// by leave-one-out on the events learned from: γ in closed form, γ = S(N − 1)/(N − S) for S the rules seen once (it maximises Σ log P₀
// with each event left out), and α₁, α₂ by a search over a grid of 121 values from 10⁻³ to 10³, each event left out of every count it
// enters. A grammar is judged on events it did not learn from by their mean −log₂ P at each level (bits per event).
//
// *Caveat.* A rule determines its match (the participants it consumes and the lenses it closes are part of it), so level 0 predicts
// what an event touches as well as what it does, and its bits exceed level 1's by the entropy of the matches, whatever the record:
// that difference is no finding. Levels 1 and 2 predict the same thing, the event's action given what it consumes; level 1 is the
// context-free grammar (a production depends on its left side alone) and level 2 the context-sensitive one, and their difference on
// held-out events is the test of whether the surroundings and the states tell what an event does.
//
// Descriptive measures on a set of events (plug-in estimates, in bits, biased low where contexts are many and their events few; the
// held-out and prequential bits are the tests): H(R); H(R | M), how far what an event consumes leaves its action open; H(R | C); and the context
// beyond the match, H(R | M) − H(R | C). (I(R; M) = H(R) − H(R | M) equals H(M), the rule determining its match.) Rates (the event operator's couplings): the events of a rule applied to one
// participant of a kind, over the time that participants of that kind were alive (their segments' durations); rules that consume
// nothing, over the record's span. Comparisons of several grammars on several measures are a partial order by dominance: one lies
// above another when it is at least as good on every measure and better on one.

#pragma once

#include <syntaxis/foams/foams.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace syntaxis::foams {

struct Grammar {
  std::size_t n = 0;                                              // events learned from
  std::map<std::string, std::size_t> rule, match, context, star;  // n(r), n(m), n(c), and the stars
  std::map<std::pair<std::string, std::string>, std::size_t> match_rule, context_rule;  // n(m, r), n(c, r)
  std::map<std::string, std::string> match_of;                    // context ↦ its match
  std::map<std::string, std::size_t> example;                     // rule ↦ the first event it was seen at
  double gamma = 1, alpha1 = 1, alpha2 = 1;
};

namespace detail {
inline double log_grid(std::size_t i) { return std::pow(10.0, -3.0 + 0.05 * static_cast<double>(i)); }  // i = 0..120
inline std::size_t count(const std::map<std::string, std::size_t>& m, const std::string& k) {
  const auto it = m.find(k);
  return it == m.end() ? 0 : it->second;
}
inline std::size_t count(const std::map<std::pair<std::string, std::string>, std::size_t>& m, const std::string& a, const std::string& b) {
  const auto it = m.find({a, b});
  return it == m.end() ? 0 : it->second;
}
}  // namespace detail

// P₀, P₁ and P₂ of a rule (a rule not seen reads as "new").
inline double p0(const Grammar& G, const std::string& r) {
  const double N = static_cast<double>(G.n), c = static_cast<double>(detail::count(G.rule, r));
  return c > 0 ? c / (N + G.gamma) : G.gamma / (N + G.gamma);
}
inline double p1(const Grammar& G, const std::string& r, const std::string& m) {
  return (static_cast<double>(detail::count(G.match_rule, m, r)) + G.alpha1 * p0(G, r)) / (static_cast<double>(detail::count(G.match, m)) + G.alpha1);
}
inline double p2(const Grammar& G, const std::string& r, const std::string& m, const std::string& c) {
  return (static_cast<double>(detail::count(G.context_rule, c, r)) + G.alpha2 * p1(G, r, m)) / (static_cast<double>(detail::count(G.context, c)) + G.alpha2);
}
inline double probability(const Grammar& G, const Star& s, int level) {
  return level <= 0 ? p0(G, s.rule) : level == 1 ? p1(G, s.rule, s.match) : p2(G, s.rule, s.match, s.context);
}

// The grammar of the events at the given vertices (none of them a moment at the ends), its escape and concentrations by
// leave-one-out.
inline Grammar learn(const std::vector<Star>& S, const std::vector<std::size_t>& vertices) {
  Grammar G;
  for (auto v : vertices) {
    detail::require(v < S.size(), "a vertex of the foam");
    const Star& s = S[v];
    detail::require(!s.boundary, "a moment at the ends of the record is no event to learn from");
    G.n++;
    G.rule[s.rule]++, G.match[s.match]++, G.context[s.context]++, G.star[s.star]++;
    G.match_rule[{s.match, s.rule}]++, G.context_rule[{s.context, s.rule}]++;
    G.match_of.emplace(s.context, s.match);
    G.example.emplace(s.rule, v);
  }
  const double N = static_cast<double>(G.n);
  if (G.n < 2) return G;
  // γ by leave-one-out, in closed form.
  std::size_t singles = 0;
  for (const auto& [r, c] : G.rule) singles += c == 1;
  if (singles == 0) G.gamma = 1e-3;
  else if (singles == G.n) G.gamma = 1e6;
  else G.gamma = std::max(1e-3, static_cast<double>(singles) * (N - 1) / (N - static_cast<double>(singles)));
  // P₀ with one event of rule r left out.
  const auto p0_out = [&](const std::string& r) {
    const double c = static_cast<double>(G.rule.at(r));
    return c >= 2 ? (c - 1) / (N - 1 + G.gamma) : G.gamma / (N - 1 + G.gamma);
  };
  const auto p1_out = [&](const std::string& r, const std::string& m, double a1) {
    return (static_cast<double>(G.match_rule.at({m, r})) - 1 + a1 * p0_out(r)) / (static_cast<double>(G.match.at(m)) - 1 + a1);
  };
  double best = -std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i <= 120; i++) {
    const double a = detail::log_grid(i);
    double l = 0;
    for (const auto& [mr, c] : G.match_rule) l += static_cast<double>(c) * std::log(p1_out(mr.second, mr.first, a));
    if (l > best) best = l, G.alpha1 = a;
  }
  best = -std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i <= 120; i++) {
    const double a = detail::log_grid(i);
    double l = 0;
    for (const auto& [cr, c] : G.context_rule) {
      const auto& m = G.match_of.at(cr.first);
      const double lower = p1_out(cr.second, m, G.alpha1);
      l += static_cast<double>(c) * std::log((static_cast<double>(c) - 1 + a * lower) / (static_cast<double>(G.context.at(cr.first)) - 1 + a));
    }
    if (l > best) best = l, G.alpha2 = a;
  }
  return G;
}

// One more event counted in a grammar (its escape and concentrations left as they are).
inline void observe(Grammar& G, const Star& s) {
  detail::require(!s.boundary, "a moment at the ends of the record is no event to learn from");
  G.n++;
  G.rule[s.rule]++, G.match[s.match]++, G.context[s.context]++, G.star[s.star]++;
  G.match_rule[{s.match, s.rule}]++, G.context_rule[{s.context, s.rule}]++;
  G.match_of.emplace(s.context, s.match);
  G.example.emplace(s.rule, s.vertex);
}

// The prequential judgement (Dawid 1984, "Statistical theory: the prequential approach", J. R. Stat. Soc. A 147, 278–292): each
// event, in the order given (the order of time), is predicted by the grammar of every event before it, and then counted; the escape and
// the concentrations are chosen again by leave-one-out on the events before each of the times given. The sum of −log₂ P is the length
// of a code for the record's rules learned as it unfolds (a prequential code in the sense of the minimum description length, Grünwald
// 2007, The Minimum Description Length Principle, MIT Press), and a grammar whose context tells what an event does shortens it. A
// rule is charged once as "new" when it first appears and is known after; which new rule it is costs nothing here, the same for
// every level.
struct Prequential {
  std::vector<std::array<double, 3>> bits;  // per event, in the order given
  std::array<double, 3> total{};
  std::size_t refits = 0;
  Grammar grammar;                          // the counts of every event, with the last escape and concentrations
};

inline Prequential prequential(const Foam& F, const std::vector<Star>& S, const std::vector<std::size_t>& events, const std::vector<double>& refit_at) {
  Prequential P;
  Grammar& G = P.grammar;
  std::vector<std::size_t> past;
  std::size_t next = 0;
  for (std::size_t k = 0; k < events.size(); k++) {
    const auto v = events[k];
    detail::require(k == 0 || F.time[v] >= F.time[events[k - 1]], "the events in the order of time");
    while (next < refit_at.size() && refit_at[next] <= F.time[v]) {
      if (past.size() >= 2) {
        const auto fit = learn(S, past);
        G.gamma = fit.gamma, G.alpha1 = fit.alpha1, G.alpha2 = fit.alpha2;
        P.refits++;
      }
      next++;
    }
    std::array<double, 3> b{};
    for (int lvl = 0; lvl < 3; lvl++) b[static_cast<std::size_t>(lvl)] = -std::log2(probability(G, S[v], lvl));
    for (std::size_t lvl = 0; lvl < 3; lvl++) P.total[lvl] += b[lvl];
    P.bits.push_back(b);
    observe(G, S[v]);
    past.push_back(v);
  }
  return P;
}

// A grammar learned from some events and judged on others.
struct HeldOut {
  std::size_t train = 0, test = 0;
  std::size_t rules_seen = 0, matches_seen = 0, contexts_seen = 0;  // test events whose rule, match, context the training saw
  std::array<double, 3> bits{};                                     // mean −log₂ P per test event, at levels 0, 1, 2
  std::array<std::size_t, 3> parameters{};                          // the counts each level reads: rules; (match, rule); (context, rule)
  Grammar grammar;
};

inline HeldOut held_out(const std::vector<Star>& S, const std::vector<std::size_t>& train, const std::vector<std::size_t>& test) {
  HeldOut H;
  H.grammar = learn(S, train);
  const Grammar& G = H.grammar;
  H.train = train.size(), H.test = test.size();
  H.parameters = {G.rule.size(), G.match_rule.size(), G.context_rule.size()};
  for (auto v : test) {
    const Star& s = S[v];
    detail::require(!s.boundary, "a moment at the ends of the record is no event to judge");
    H.rules_seen += G.rule.count(s.rule) > 0, H.matches_seen += G.match.count(s.match) > 0, H.contexts_seen += G.context.count(s.context) > 0;
    for (int k = 0; k < 3; k++) H.bits[static_cast<std::size_t>(k)] -= std::log2(probability(G, s, k));
  }
  if (!test.empty())
    for (auto& b : H.bits) b /= static_cast<double>(test.size());
  return H;
}

// Plug-in entropies on a set of events, in bits.
struct Entropies {
  std::size_t events = 0, rules = 0, matches = 0, contexts = 0, stars = 0;
  double rule = 0, rule_given_match = 0, rule_given_context = 0;
  double information() const { return rule - rule_given_match; }                // I(R; M), which is H(M)
  double context_beyond() const { return rule_given_match - rule_given_context; }  // what the context adds
};

inline Entropies entropies(const std::vector<Star>& S, const std::vector<std::size_t>& vertices) {
  Entropies E;
  std::map<std::string, std::size_t> r, m, c, st;
  std::map<std::pair<std::string, std::string>, std::size_t> mr, cr;
  for (auto v : vertices) {
    const Star& s = S[v];
    if (s.boundary) continue;
    E.events++;
    r[s.rule]++, m[s.match]++, c[s.context]++, st[s.star]++, mr[{s.match, s.rule}]++, cr[{s.context, s.rule}]++;
  }
  E.rules = r.size(), E.matches = m.size(), E.contexts = c.size(), E.stars = st.size();
  const double N = static_cast<double>(E.events);
  if (N == 0) return E;
  for (const auto& [k, n] : r) E.rule -= n / N * std::log2(n / N);
  for (const auto& [k, n] : mr) E.rule_given_match -= n / N * std::log2(static_cast<double>(n) / static_cast<double>(m[k.first]));
  for (const auto& [k, n] : cr) E.rule_given_context -= n / N * std::log2(static_cast<double>(n) / static_cast<double>(c[k.first]));
  return E;
}

// The match key of an event that consumes one participant of a kind (with no state), and of one that consumes nothing.
inline std::string single_match(const std::string& kind) {
  mc::LabelledGraph G;
  G.vertex = {"m", "x:" + kind};
  G.edges.push_back({0, 1, "at"});
  return detail::key_any(G).code;
}
inline std::string empty_match() {
  mc::LabelledGraph G;
  G.vertex = {"m"};
  return detail::key_any(G).code;
}

// The time participants of each kind were alive (the durations of their segments), and the record's span.
struct Exposure {
  std::map<std::string, double> by_kind;
  double span = 0;
};

inline Exposure exposure(const Foam& F) {
  Exposure X;
  for (const auto& sg : F.W.segments) X.by_kind[F.participant[sg.participant]] += F.time[sg.to] - F.time[sg.from];
  if (F.size() > 0) X.span = F.time.back() - F.time.front();
  return X;
}

// Rates of the rules applied to one participant of a kind (per unit of its alive time) and of the rules that consume nothing (per
// unit of the record's span), from a grammar's counts.
struct Rate {
  std::string rule, consumes;  // consumes: a participant kind, or empty for nothing
  std::size_t events = 0;
  double exposure = 0, rate = 0;
};

inline std::vector<Rate> rates(const Grammar& G, const Exposure& X) {
  std::vector<Rate> out;
  std::map<std::string, std::string> kind_of;
  for (const auto& [kind, t] : X.by_kind) kind_of[single_match(kind)] = kind;
  const std::string none = empty_match();
  for (const auto& [mr, n] : G.match_rule) {
    Rate R;
    R.rule = mr.second, R.events = n;
    if (mr.first == none) {
      R.exposure = X.span;
    } else if (const auto it = kind_of.find(mr.first); it != kind_of.end()) {
      R.consumes = it->second, R.exposure = X.by_kind.at(it->second);
    } else {
      continue;  // a match of several participants: no exposure is counted for it
    }
    R.rate = R.exposure > 0 ? static_cast<double>(n) / R.exposure : 0;
    out.push_back(R);
  }
  return out;
}

// The partial order of candidates by dominance on measures, each lower being better.
struct Dominance {
  std::vector<std::pair<std::size_t, std::size_t>> above;  // (i, j): i at least as good everywhere and better somewhere
  std::vector<std::size_t> undominated;
};

inline Dominance dominance(const std::vector<std::vector<double>>& measures) {
  Dominance D;
  const std::size_t n = measures.size();
  for (std::size_t i = 0; i < n; i++) detail::require(measures[i].size() == measures[0].size(), "the same measures for every candidate");
  std::vector<char> dominated(n, 0);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) {
      if (i == j) continue;
      bool no_worse = true, better = false;
      for (std::size_t k = 0; k < measures[i].size(); k++) {
        if (measures[i][k] > measures[j][k]) no_worse = false;
        if (measures[i][k] < measures[j][k]) better = true;
      }
      if (no_worse && better) D.above.push_back({i, j}), dominated[j] = 1;
    }
  for (std::size_t i = 0; i < n; i++)
    if (!dominated[i]) D.undominated.push_back(i);
  return D;
}

}  // namespace syntaxis::foams
