// Checks of syntaxis/foams: grammars whose words are foams.
//
// - A foam worked by hand: three buildings (two sharing a wall), a street between two junctions; the third building begun in 2010,
//   the first altered in 2015, two events on the street (one meeting its two cross-sections, the junctions present), censuses in
//   2004 and 2012 (the third building's in 2012 only). Its 13 vertices, 21 segments and 10 lenses; 8 kinds of star, 7 rules, 3
//   matches and 5 contexts among its 11 events; κ = 0 at every event and 1 at each end, Σκ = χ = 2; Foata steps of widths 3, 4, 2,
//   1, 1; the entropies H(R) and H(R | M) by hand; the rates per participant-year (the census 5/61 per building-year); its lives
//   (16 symbols, two rules of Re-Pair, the grammar as large as its input), and with the censuses left out (five events out, one
//   rule).
// - Three worldlines pinched at one event that passes through them all: κ = −2 there, Σκ = χ = −3.
// - On 200 random records of encounters: the cells counted from the stars equal the foam's and Σκ = V − E + F; every key unchanged
//   when the participants are renamed; the Foata steps equal morphe's causal order of the trace; at every match and context the
//   probabilities of the rules seen and of a new rule sum to one; γ by its closed form at least as good, by leave-one-out, as any
//   value of a fine grid.
// - A planted context: events whose kind is set by their participant's state, 5% of them flipped. The context-sensitive grammar
//   (level 2) predicts held-out events in under 0.5 bits (the flips cost h(0.05) = 0.29) where the context-free one (level 1) needs
//   nearly a bit; when the kind is drawn apart from the
//   state, level 2 costs at most 0.05 bits more than level 1 (its backoff learns that the context tells nothing).
// - The prequential code worked by hand on three events; on the planted context, level 2's code a third of level 1's or less, and
//   within 2% when the kind is drawn apart from the state.
// - The dominance order of three candidates on two measures; refusals.

#include <syntaxis/foams/foams.hpp>
#include <syntaxis/foams/grammar.hpp>
#include <syntaxis/foams/lives.hpp>

#include <morphe/formal/traces.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace syntaxis::foams;
namespace mt = morphe::topology;

namespace {
template <class F>
bool refused(F f) {
  try {
    f();
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

// The hand foam: participants a, b, c (buildings), s (a street), j0, j1 (junctions); pairs: the party wall (a, b), frontages (a, s),
// (b, s), (c, s), cross-sections (s, j0), (s, j1).
Foam hand() {
  mt::Encounters E;
  E.participants = 6;
  E.present = {{0, 1, 3, 4, 5}, {0}, {1}, {2, 3}, {0}, {1}, {2}, {0}, {3}, {3}, {3, 4, 5}, {3, 4, 5}, {0, 1, 2, 3, 4, 5}};
  E.pairs = {{0, 1, {0, 12}}, {0, 3, {0, 12}}, {1, 3, {0, 12}}, {2, 3, {3, 12}}, {3, 4, {0, 10, 11, 12}}, {3, 5, {0, 10, 11, 12}}};
  std::vector<std::string> kind = {"start", "observation", "observation", "birth", "observation", "observation", "observation", "alteration",
                                   "regime", "regime ended", "regime", "regime ended", "end"};
  std::vector<double> t = {2002, 2004, 2004, 2010, 2012, 2012, 2012, 2015, 2016.2, 2016.21, 2017.5, 2017.51, 2025};
  std::vector<char> b(13, 0);
  b[0] = b[12] = 1;
  return make(E, kind, t, b, {"building", "building", "building", "street", "junction", "junction"},
              {"party-wall", "frontage", "frontage", "frontage", "cross-section", "cross-section"});
}

mt::Encounters random_encounters(std::mt19937& g, std::size_t P, std::size_t V) {
  mt::Encounters E;
  E.participants = P;
  E.present.assign(V, {});
  for (std::size_t v = 0; v < V; v++) {
    for (std::size_t x = 0; x < P; x++)
      if (g() % 10 < 4) E.present[v].push_back(x);
    if (E.present[v].empty()) E.present[v].push_back(g() % P);
  }
  for (std::size_t t = 0; t < P; t++) {
    std::size_t a = g() % P, c = g() % P;
    if (a == c) continue;
    std::vector<std::size_t> common;
    for (std::size_t v = 0; v < V; v++) {
      const auto& p = E.present[v];
      if (std::find(p.begin(), p.end(), a) != p.end() && std::find(p.begin(), p.end(), c) != p.end()) common.push_back(v);
    }
    std::vector<std::size_t> meet;
    for (auto v : common)
      if (g() % 10 < 6) meet.push_back(v);
    if (meet.size() >= 2) E.pairs.push_back({a, c, meet});
  }
  return E;
}
}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- The hand foam ------------------------------------------------------------------------------------------------------------
  {
    const auto F = hand();
    const auto S = stars(F);
    std::set<std::string> star, rule, match, context;
    for (const auto& st : S)
      if (!st.boundary) star.insert(st.star), rule.insert(st.rule), match.insert(st.match), context.insert(st.context);
    s.require("a foam by hand: 13 vertices, 21 segments, 10 lenses", F.size() == 13 && F.W.segments.size() == 21 && F.W.lenses.size() == 10);
    s.require("... 8 kinds of star, 7 rules, 3 matches, 5 contexts among its 11 events",
              star.size() == 8 && rule.size() == 7 && match.size() == 3 && context.size() == 5,
              std::to_string(star.size()) + " " + std::to_string(rule.size()) + " " + std::to_string(match.size()) + " " + std::to_string(context.size()));
    s.require("... the censuses of the two walled buildings one star; the third building's census another star, the same rule",
              S[1].star == S[2].star && S[1].star == S[4].star && S[1].star == S[5].star && S[6].star != S[1].star && S[6].rule == S[1].rule);
    s.require("... the birth begins the third building beside the street through, opens their frontage, and four lenses pass on the street",
              S[3].opens == 1 && S[3].passes == 4 && S[3].ends == 3 && S[3].in[0] == 0 && S[3].out[0] == 1);
    s.require("... the closure closes and opens both cross-sections, three frontages passing on the street",
              S[10].closes == 2 && S[10].opens == 2 && S[10].passes == 3 && S[10].ends == 6);
    bool zero = true;
    for (std::size_t v = 1; v + 1 < S.size(); v++) zero = zero && S[v].kappa == 0;
    const auto R = replay(F, S);
    s.require("... κ = 0 at every event and 1 at each end; Σκ = χ = 13 − 21 + 10 = 2 (Proposition 7.1)",
              zero && S[0].kappa == 1 && S[12].kappa == 1 && R.consistent && R.euler == 2 && R.local == 2 && R.edges == 21 && R.faces == 10);
    const auto T = steps(F);
    s.require("... the trace's Foata steps: depth 5, widths 3, 4, 2, 1, 1", T.depth == 5 && T.width == std::vector<std::size_t>{3, 4, 2, 1, 1});
    std::vector<std::size_t> events;
    for (std::size_t v = 1; v + 1 < S.size(); v++) events.push_back(v);
    const auto H = entropies(S, events);
    const double hR = -(5.0 / 11) * std::log2(5.0 / 11) - 6 * (1.0 / 11) * std::log2(1.0 / 11);
    const double hRM = (6.0 / 11) * (-(5.0 / 6) * std::log2(5.0 / 6) - (1.0 / 6) * std::log2(1.0 / 6)) + (3.0 / 11) * std::log2(3.0) + (2.0 / 11) * 1.0;
    s.check("... H(R) and H(R | M) by hand", std::fabs(H.rule - hR) + std::fabs(H.rule_given_match - hRM), 1e-12);
    const auto G = learn(S, events);
    const auto X = exposure(F);
    double census_rate = -1;
    for (const auto& r : rates(G, X))
      if (r.rule == S[1].rule && r.consumes == "building") census_rate = r.rate;
    s.check("... buildings alive 61 years in all; the census 5 events per 61 building-years", std::fabs(X.by_kind.at("building") - 61) + std::fabs(census_rate - 5.0 / 61), 1e-12);
    const auto L = lives(F, S);
    std::size_t symbols = 0;
    for (const auto& w : L.words) symbols += w.size();
    s.require("... its lives: 16 symbols in 6 words, two rules of Re-Pair, the grammar as large as its input, every life recovered",
              symbols == 16 && L.words.size() == 6 && L.grammar.rules.size() == 2 && L.grammar.size() == 16 && L.grammar.expand() == L.words);
    const auto C = lives(F, S, {"observation"});
    s.require("... with the censuses left out: five events out, five lives, one rule (the junctions' closing and reopening)",
              C.left_out == 5 && C.words.size() == 5 && C.grammar.rules.size() == 1 && C.uses[0] == 2);
    s.require("... a readable account of the birth", describe(F, S[3]) == "birth: building begins, street through; opens frontage", describe(F, S[3]));
  }

  // ---- Worldlines pinched -------------------------------------------------------------------------------------------------------
  {
    // Three participants present at the two ends and together at one event between, no pair: the event has κ = 1 − 3 = −2, the ends
    // −1/2 each, and χ = 3 − 6 + 0 = −3.
    mt::Encounters E;
    E.participants = 3;
    E.present = {{0, 1, 2}, {0, 1, 2}, {0, 1, 2}};
    const auto F = make(E, {"start", "e", "end"}, {0, 1, 2}, {1, 0, 1}, {"k", "k", "k"}, {});
    const auto S = stars(F);
    const auto R = replay(F, S);
    s.require("three worldlines pinched at one event passing through: κ = 1 − 3 = −2, the ends −1/2 each, Σκ = χ = 3 − 6 = −3",
              S[1].kappa == -2 && S[0].kappa == -0.5 && S[2].kappa == -0.5 && R.euler == -3 && R.local == -3);
  }

  // ---- The prequential code by hand -------------------------------------------------------------------------------------------
  {
    // One participant through three events of kinds a, a, b, with γ = α₁ = 1 (no refit): at level 0 the first event is new with
    // probability 1, the second a with 1/2, the third new with 1/3, 0 + 1 + log₂ 3 bits; at level 1 the second is (1 + 1/2)/2 = 3/4 and
    // the third (0 + 1/3)/3 = 1/9, log₂(4/3) + log₂ 9 bits.
    mt::Encounters E;
    E.participants = 1;
    E.present = {{0}, {0}, {0}, {0}, {0}};
    const auto F = make(E, {"start", "a", "a", "b", "end"}, {0, 1, 2, 3, 4}, {1, 0, 0, 0, 1}, {"k"}, {});
    const auto S = stars(F);
    const auto P = prequential(F, S, {1, 2, 3}, {});
    s.check("the prequential code by hand: 1 + log₂ 3 bits at level 0, log₂(4/3) + log₂ 9 at level 1",
            std::fabs(P.total[0] - (1 + std::log2(3.0))) + std::fabs(P.total[1] - (std::log2(4.0 / 3) + std::log2(9.0))), 1e-12);
  }

  // ---- Random records -----------------------------------------------------------------------------------------------------------
  {
    std::mt19937 g(20261126);
    bool cells = true, renamed = true, foata = true, sums = true, gamma_best = true;
    const std::vector<std::string> kinds = {"p", "q", "r"}, pairs = {"u", "w"}, events = {"e", "f", "h"};
    for (int t = 0; t < 200; t++) {
      const std::size_t P = 2 + g() % 5, V = 3 + g() % 12;
      const auto E = random_encounters(g, P, V);
      std::vector<std::string> vk(V), pk(P), qk(E.pairs.size());
      for (auto& k : vk) k = events[g() % 3];
      for (auto& k : pk) k = kinds[g() % 3];
      for (auto& k : qk) k = pairs[g() % 2];
      std::vector<double> time(V);
      std::iota(time.begin(), time.end(), 0.0);
      const auto F = make(E, vk, time, std::vector<char>(V, 0), pk, qk);
      const auto S = stars(F);
      const auto R = replay(F, S);
      cells = cells && R.consistent && std::fabs(R.local - static_cast<double>(R.euler)) < 1e-12 &&
              R.euler == static_cast<long>(V) - static_cast<long>(F.W.segments.size()) + static_cast<long>(F.W.lenses.size());
      // Rename the participants.
      std::vector<std::size_t> pi(P);
      std::iota(pi.begin(), pi.end(), 0);
      std::shuffle(pi.begin(), pi.end(), g);
      mt::Encounters E2;
      E2.participants = P;
      for (const auto& pr : E.present) {
        std::vector<std::size_t> q;
        for (auto x : pr) q.push_back(pi[x]);
        std::shuffle(q.begin(), q.end(), g);
        E2.present.push_back(q);
      }
      for (const auto& pr : E.pairs) E2.pairs.push_back({pi[pr.a], pi[pr.b], pr.meetings});
      std::vector<std::string> pk2(P);
      for (std::size_t x = 0; x < P; x++) pk2[pi[x]] = pk[x];
      const auto F2 = make(E2, vk, time, std::vector<char>(V, 0), pk2, qk);
      const auto S2 = stars(F2);
      for (std::size_t v = 0; v < V; v++)
        renamed = renamed && S[v].star == S2[v].star && S[v].rule == S2[v].rule && S[v].match == S2[v].match && S[v].context == S2[v].context;
      // Foata steps against morphe's causal order of the trace.
      std::vector<morphe::formal::traces::Footprint> fp(V);
      for (std::size_t v = 0; v < V; v++)
        for (auto x : E.present[v]) fp[v].reads.insert(std::to_string(x)), fp[v].writes.insert(std::to_string(x));
      const auto order = morphe::formal::traces::causal_order(fp);
      const auto T = steps(F);
      for (std::size_t v = 0; v < V; v++) foata = foata && T.step[v] == static_cast<long>(order.step[v]);
      // The probabilities sum to one at every match and context.
      std::vector<std::size_t> all(V);
      std::iota(all.begin(), all.end(), 0);
      const auto G = learn(S, all);
      for (const auto& [m, n] : G.match) {
        double sum = p1(G, "a rule never seen", m);
        for (const auto& [r, c] : G.rule) sum += p1(G, r, m);
        sums = sums && std::fabs(sum - 1) < 1e-12;
      }
      for (const auto& [c, n] : G.context) {
        const auto& m = G.match_of.at(c);
        double sum = p2(G, "a rule never seen", m, c);
        for (const auto& [r, k] : G.rule) sum += p2(G, r, m, c);
        sums = sums && std::fabs(sum - 1) < 1e-12;
      }
      // γ by its closed form against a fine grid, by leave-one-out.
      if (G.n >= 2) {
        const double N = static_cast<double>(G.n);
        const auto loo = [&](double gm) {
          double l = 0;
          for (const auto& [r, c] : G.rule) l += static_cast<double>(c) * std::log(c >= 2 ? (c - 1) / (N - 1 + gm) : gm / (N - 1 + gm));
          return l;
        };
        double grid = -1e300;
        for (int i = 0; i <= 400; i++) grid = std::max(grid, loo(std::pow(10.0, -3 + 0.02 * i)));
        gamma_best = gamma_best && loo(G.gamma) >= grid - 1e-9;
      }
    }
    s.require("200 random records: the cells counted from the stars are the foam's, and Σκ = V − E + F", cells);
    s.require("... every key unchanged when the participants are renamed", renamed);
    s.require("... the Foata steps equal morphe's causal order of the trace", foata);
    s.require("... at every match and context, the rules seen and a new rule sum to one", sums);
    s.require("... γ in closed form at least as good by leave-one-out as any value of a fine grid", gamma_best);
  }

  // ---- A planted context --------------------------------------------------------------------------------------------------------
  {
    calx::random::Stream rs(20261126, 1, 2, 3);
    const auto run = [&](bool planted) {
      // Every participant present at the two moments at the ends, as in a record of a window, so that each event between is
      // passed through; 1,200 events, each touching one participant.
      const std::size_t P = 40, V = 1202;
      mt::Encounters E;
      E.participants = P;
      std::vector<std::string> vk(V), state(P);
      std::vector<char> ends(V, 0);
      for (std::size_t x = 0; x < P; x++) state[x] = x % 2 ? "A" : "B";
      std::vector<std::size_t> everyone(P);
      std::iota(everyone.begin(), everyone.end(), 0);
      E.present.push_back(everyone), vk[0] = "start", ends[0] = 1;
      for (std::size_t v = 1; v + 1 < V; v++) {
        const std::size_t x = rs.below(P);
        E.present.push_back({x});
        const bool a = planted ? state[x] == "A" : rs.uniform() < 0.5;
        vk[v] = (rs.uniform() < 0.05) != a ? "e1" : "e2";
      }
      E.present.push_back(everyone), vk[V - 1] = "end", ends[V - 1] = 1;
      std::vector<double> time(V);
      std::iota(time.begin(), time.end(), 0.0);
      auto F = make(E, vk, time, ends, std::vector<std::string>(P, "k"), {});
      F.states.assign(P, {});
      for (std::size_t x = 0; x < P; x++) F.states[x][0] = state[x];
      const auto S = stars(F);
      std::vector<std::size_t> train, test;
      for (std::size_t v = 1; v + 1 < V; v++) (v <= 840 ? train : test).push_back(v);
      return std::make_pair(held_out(S, train, test), prequential(F, S, [&] { std::vector<std::size_t> e; for (std::size_t v = 1; v + 1 < V; v++) e.push_back(v); return e; }(), {200.0, 400.0, 600.0, 800.0, 1000.0}));
    };
    const auto [P, Pq] = run(true);
    const auto [C, Cq] = run(false);
    s.require("a planted context: the kind set by the participant's state; level 2 under 0.5 bits on held-out events, level 1 near a bit",
              P.bits[2] < 0.5 && P.bits[1] > 0.9, std::to_string(P.bits[1]) + " / " + std::to_string(P.bits[2]));
    s.require("... drawn apart from the state, level 2 at most 0.05 bits above level 1",
              C.bits[2] <= C.bits[1] + 0.05, std::to_string(C.bits[1]) + " / " + std::to_string(C.bits[2]));
    s.require("... prequentially (refitted every 200 events), level 2's code a third of level 1's or less when planted; within 2% when not",
              Pq.total[2] < Pq.total[1] / 3 && Cq.total[2] < Cq.total[1] * 1.02,
              std::to_string(Pq.total[1]) + " / " + std::to_string(Pq.total[2]) + "; " + std::to_string(Cq.total[1]) + " / " + std::to_string(Cq.total[2]));
  }

  // ---- Dominance and refusals ---------------------------------------------------------------------------------------------------
  {
    const auto D = dominance({{1.0, 3}, {2.0, 1}, {2.5, 3}});
    s.require("dominance on two measures: the first above the third; the first and second undominated, incomparable",
              D.above == std::vector<std::pair<std::size_t, std::size_t>>{{0, 2}, {1, 2}} && D.undominated == std::vector<std::size_t>{0, 1});
    const auto F = hand();
    const auto S = stars(F);
    s.require("a moment at the ends is no event to learn from", refused([&] { learn(S, {0, 1}); }));
    mt::Encounters E;
    E.participants = 1;
    E.present = {{0}, {0}};
    s.require("times out of order are refused", refused([&] { make(E, {"e", "e"}, {1, 0}, {0, 0}, {"k"}, {}); }));
  }
  return s.done();
}
