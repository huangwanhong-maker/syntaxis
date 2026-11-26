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
#include <syntaxis/foams/grow.hpp>
#include <syntaxis/foams/nets.hpp>

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

  // ---- The process net and its patterns -----------------------------------------------------------------------------------------
  {
    // Three streets and a building over a window: four episodes on the three streets, each a regime and its end, a census of the
    // building before each. With the censuses left out: 8 events and 21 wires; the episode (a regime, its end, three wires) found
    // first, four times; then two episodes in a row, twice; two boxes and three wires left: 15 against 29.
    mt::Encounters E;
    E.participants = 4;
    std::vector<std::string> vk = {"start"};
    std::vector<double> t = {0};
    E.present.push_back({0, 1, 2, 3});
    for (int k = 0; k < 4; k++) {
      E.present.push_back({3}), vk.push_back("observation");
      E.present.push_back({0, 1, 2}), vk.push_back("regime");
      E.present.push_back({0, 1, 2}), vk.push_back("regime ended");
    }
    E.present.push_back({0, 1, 2, 3}), vk.push_back("end");
    t.resize(vk.size());
    std::iota(t.begin(), t.end(), 0.0);
    std::vector<char> ends(vk.size(), 0);
    ends.front() = ends.back() = 1;
    const auto F = make(E, vk, t, ends, {"street", "street", "street", "building"}, {});
    const auto S = stars(F);
    const auto N0 = net(F, S, {"observation"});
    s.require("a process net: the four censuses left out, 8 events, 21 wires, 6 open ends",
              N0.left_out == 4 && N0.events == 8 && N0.wires_made == 21 && N0.open_ends == 6,
              std::to_string(N0.events) + " / " + std::to_string(N0.wires_made) + " / " + std::to_string(N0.open_ends));
    const auto N = compressed(F, S, {"observation"});
    const bool episode = N.rules.size() == 2 && N.rules[0].uses == 4 && N.rules[0].wires.size() == 3 && N.rules[0].width == 3 && N.rules[0].length == 2 &&
                         N.labels[N.rules[0].first] == S[2].rule && N.labels[N.rules[0].second] == S[3].rule;
    s.require("... compressed: the episode (a regime, its end, three wires) four times, then two episodes in a row twice", episode && N.rules[1].uses == 2 && N.rules[1].length == 4,
              std::to_string(N.rules.size()) + " rules");
    s.require("... two boxes and three wires left; the description 15 against 29; the pattern derives regime, end, regime, end",
              N.boxes.size() == 2 && N.wires.size() == 3 && N.size_before() == 29 && N.size_after() == 15 &&
                  derive(N, static_cast<std::uint32_t>(N.terminals + 1)) == std::vector<std::uint32_t>{N.rules[0].first, N.rules[0].second, N.rules[0].first, N.rules[0].second});
    // On random records: every event in one box of the top level; the wires inside the patterns and those left add up; each
    // pattern derives as many events as its instances hold; the same patterns when the participants are renamed.
    std::mt19937 g(7);
    bool covered = true, wires = true, lengths = true, renamed = true;
    for (int trial = 0; trial < 200; trial++) {
      const std::size_t P = 2 + g() % 5, V = 4 + g() % 30;
      const auto R = random_encounters(g, P, V);
      std::vector<std::string> kinds(V), pk(P), qk(R.pairs.size());
      for (auto& k : kinds) k = std::string(1, static_cast<char>('e' + g() % 2));
      for (auto& k : pk) k = std::string(1, static_cast<char>('p' + g() % 2));
      for (auto& k : qk) k = "u";
      std::vector<double> time(V);
      std::iota(time.begin(), time.end(), 0.0);
      const auto Fr = make(R, kinds, time, std::vector<char>(V, 0), pk, qk);
      const auto Sr = stars(Fr);
      const auto Nr = compressed(Fr, Sr);
      std::vector<int> seen(V, 0);
      for (const auto& b : Nr.boxes)
        for (auto v : b.events) seen[v]++;
      covered = covered && std::all_of(seen.begin(), seen.end(), [](int k) { return k == 1; });
      std::size_t inside = 0;
      for (const auto& r : Nr.rules) inside += r.uses * r.wires.size();
      wires = wires && inside + Nr.wires.size() == Nr.wires_made;
      for (std::size_t i = 0; i < Nr.rules.size(); i++) lengths = lengths && derive(Nr, static_cast<std::uint32_t>(Nr.terminals + i)).size() == Nr.rules[i].length;
      // Renamed participants: the same rules, uses and wires.
      std::vector<std::size_t> pi(P);
      std::iota(pi.begin(), pi.end(), 0);
      std::shuffle(pi.begin(), pi.end(), g);
      mt::Encounters R2;
      R2.participants = P;
      for (const auto& pr : R.present) {
        std::vector<std::size_t> q;
        for (auto x : pr) q.push_back(pi[x]);
        R2.present.push_back(q);
      }
      for (const auto& pr : R.pairs) R2.pairs.push_back({pi[pr.a], pi[pr.b], pr.meetings});
      std::vector<std::string> pk2(P);
      for (std::size_t x = 0; x < P; x++) pk2[pi[x]] = pk[x];
      const auto F2 = make(R2, kinds, time, std::vector<char>(V, 0), pk2, qk);
      const auto N2 = compressed(F2, stars(F2));
      renamed = renamed && N2.rules.size() == Nr.rules.size() && N2.boxes.size() == Nr.boxes.size() && N2.wires.size() == Nr.wires.size();
      for (std::size_t i = 0; renamed && i < Nr.rules.size(); i++) renamed = N2.rules[i].uses == Nr.rules[i].uses && N2.rules[i].wires.size() == Nr.rules[i].wires.size();
    }
    s.require("200 random records: every event in one box of the top level", covered);
    s.require("... the wires inside the patterns and those left add up to the net's", wires);
    s.require("... each pattern derives as many events as its instances hold", lengths);
    s.require("... the same patterns, uses and wires when the participants are renamed", renamed);
  }

  // ---- Growing foams ------------------------------------------------------------------------------------------------------------
  {
    // The keys of an event not yet in a foam, built from what it consumes, equal those stars() reads once it is there: on the hand
    // foam and on 200 random records, each with random states.
    std::mt19937 g(11);
    bool same = true;
    const auto check_keys = [&](const Foam& Fk) {
      const auto Sk = stars(Fk);
      for (std::size_t v = 0; v < Fk.size(); v++) {
        const Star& s = Sk[v];
        if (s.boundary) continue;
        std::vector<std::size_t> position(s.participants.size(), 0);
        std::vector<std::string> kinds, carried;
        std::vector<std::vector<std::string>> passes;
        for (std::size_t i = 0; i < s.participants.size(); i++) {
          if (!s.in[i]) continue;
          position[i] = kinds.size();
          kinds.push_back(Fk.participant[s.participants[i]]), carried.push_back(Fk.state(s.participants[i], v)), passes.emplace_back();
        }
        std::vector<std::tuple<std::size_t, std::size_t, std::string>> closes;
        for (const auto& c : s.corners) {
          const auto& pair = Fk.pair[Fk.W.lenses[c.lens].pair];
          if (c.kind == Corner::pass && s.in[c.a]) passes[position[c.a]].push_back(pair);
          if (c.kind == Corner::close) closes.emplace_back(position[c.a], position[c.b], pair);
        }
        const auto [m, ctx] = match_and_context(kinds, carried, passes, closes);
        same = same && m == s.match && ctx == s.context;
      }
    };
    auto Fh = hand();
    Fh.states.assign(6, {});
    Fh.states[0][4] = "office", Fh.states[1][5] = "retail";
    check_keys(Fh);
    for (int trial = 0; trial < 200; trial++) {
      const std::size_t P = 2 + g() % 5, V = 3 + g() % 15;
      const auto R = random_encounters(g, P, V);
      std::vector<std::string> kinds(V), pk(P), qk(R.pairs.size());
      for (auto& k : kinds) k = std::string(1, static_cast<char>('e' + g() % 2));
      for (auto& k : pk) k = std::string(1, static_cast<char>('p' + g() % 2));
      for (auto& k : qk) k = std::string(1, static_cast<char>('u' + g() % 2));
      std::vector<double> time(V);
      std::iota(time.begin(), time.end(), 0.0);
      auto Fr = make(R, kinds, time, std::vector<char>(V, 0), pk, qk);
      Fr.states.assign(P, {});
      for (std::size_t x = 0; x < P; x++)
        for (std::size_t v = 0; v < V; v++)
          if (g() % 3 == 0) Fr.states[x][v] = g() % 2 ? "A" : "B";
      check_keys(Fr);
    }
    s.require("the keys of an event not yet in a foam equal those stars() reads once it is there: the hand foam and 200 random records with states", same);

    // A planted record over 20 years: 20 buildings altered at 0.2 a year each, a further building born at year 10 against the
    // first street and altered after; three streets taken together by regimes that begin at 3 a year while none of theirs is in
    // force and end 0.02 years later; the first two streets taken by regimes of their own that begin at one a year and last half a
    // year, overlapping the others; ten buildings fronting the first street. The streets carry whether a regime is in force.
    calx::random::Stream rs(20261126, 4, 5, 6);
    struct Ev {
      double t;
      std::string kind;
      std::vector<std::size_t> who;
    };
    std::vector<Ev> evs;
    const double T = 20;
    for (std::size_t b = 3; b <= 23; b++) {
      double t = b == 23 ? 10.0 : 0.0;
      while ((t += rs.exponential() / 0.2) < T) evs.push_back({t, "alteration", {b}});
    }
    std::size_t starts = 0;
    double sum_duration = 0;
    for (double t = 0;;) {
      t += rs.exponential() / 3.0;
      if (t >= T) break;
      evs.push_back({t, "regime", {0, 1, 2}}), starts++;
      if (t + 0.02 < T) evs.push_back({t + 0.02, "regime ended", {0, 1, 2}}), sum_duration += 0.02;
      t += 0.02;
    }
    std::size_t starts2 = 0;
    for (double t = 0;;) {
      t += rs.exponential() / 1.0;
      if (t >= T) break;
      evs.push_back({t, "regime", {0, 1}}), starts2++;
      if (t + 0.5 < T) evs.push_back({t + 0.5, "regime ended", {0, 1}});
      t += 0.5;
    }
    evs.push_back({10.0, "birth", {23, 0}});
    std::sort(evs.begin(), evs.end(), [](const Ev& a, const Ev& b) { return a.t < b.t; });
    mt::Encounters E;
    E.participants = 24;
    std::vector<std::string> vk = {"start"};
    std::vector<double> times = {0};
    std::vector<std::size_t> everyone(23);
    std::iota(everyone.begin(), everyone.end(), 0);
    E.present.push_back(everyone);
    std::size_t birth_vertex = 0, alterations = 0;
    for (const auto& e : evs) {
      E.present.push_back(e.who), vk.push_back(e.kind), times.push_back(e.t);
      if (e.kind == "birth") birth_vertex = E.present.size() - 1;
      alterations += e.kind == "alteration";
    }
    everyone.push_back(23);
    E.present.push_back(everyone), vk.push_back("end"), times.push_back(T);
    const std::size_t END = E.present.size() - 1;
    for (std::size_t b = 3; b <= 12; b++) E.pairs.push_back({b, 0, {0, END}});
    E.pairs.push_back({23, 0, {birth_vertex, END}});
    std::vector<char> ends(E.present.size(), 0);
    ends.front() = ends.back() = 1;
    std::vector<std::string> pk(24, "building");
    pk[0] = pk[1] = pk[2] = "street";
    auto Fp = make(E, vk, times, ends, pk, std::vector<std::string>(E.pairs.size(), "frontage"));
    Fp.states.assign(24, {});
    for (std::size_t x = 0; x < 3; x++) {
      int open = 0;
      for (auto v : Fp.W.worldline[x]) {
        if (ends[v]) continue;
        Fp.states[x][v] = open > 0 ? "under a regime" : "";
        if (vk[v] == "regime") open++;
        if (vk[v] == "regime ended") open--;
      }
    }
    const auto Sp = stars(Fp);
    const auto Pg = learn_growth(Fp, Sp);
    double alt = 0, st = 0, dur = 0, endings = 0, born = 0, st2 = 0, dur2 = 0, endings2 = 0;
    bool alternate = true, consistent = true, after_birth = true;
    std::string first_break;
    const int H = 100;
    for (int h = 0; h < H; h++) {
      const auto Fg = grow(Fp, Sp, Pg, rs);
      consistent = consistent && replay(Fg, stars(Fg)).consistent;
      std::string expect[2] = {"regime", "regime"};
      double began[2] = {0, 0}, birth_at = -1;
      for (std::size_t v = 1; v + 1 < Fg.size(); v++) {
        const auto& k = Fg.vertex[v];
        alt += k == "alteration";
        if (k == "birth") born += 1, birth_at = Fg.time[v];
        if (k == "regime" || k == "regime ended") {
          const int site = Fg.encounters.present[v].size() == 3 ? 0 : 1;  // the three streets, or the first two
          if (alternate && k != expect[site])
            first_break = "history " + std::to_string(h) + ", site " + std::to_string(site) + " at " + std::to_string(Fg.time[v]) + ": " + k + " after " + std::to_string(began[site]) +
                          "; states carried: " + Fg.state(Fg.encounters.present[v][0], v) + "|" + Fg.state(Fg.encounters.present[v][1], v);
          alternate = alternate && k == expect[site];
          if (k == "regime") (site == 0 ? st : st2) += 1, began[site] = Fg.time[v];
          else if (site == 0) dur += Fg.time[v] - began[0], endings += 1;
          else dur2 += Fg.time[v] - began[1], endings2 += 1;
          expect[site] = k == "regime" ? "regime ended" : "regime";
        }
      }
      // The born building's alterations come after its birth (its participant is the last, when it is born).
      if (birth_at >= 0)
        for (std::size_t v = 1; v + 1 < Fg.size(); v++)
          if (Fg.vertex[v] == "alteration" && Fg.encounters.present[v].front() == Fg.participant.size() - 1) after_birth = after_birth && Fg.time[v] > birth_at;
    }
    const double a = static_cast<double>(alterations);
    s.require("a planted record grown 100 times: alterations as many as the record's on average, within four standard errors",
              std::fabs(alt / H - a) < 4 * std::sqrt(a / H), std::to_string(alt / H) + " against " + std::to_string(alterations));
    s.require("... every regime and its end alternate at each of the two sites, though they overlap: each site reads its own states", alternate, first_break);
    s.require("... regimes begun as often as in the record at each site, within 10% and 15%, lasting as long on average, within 15%",
              std::fabs(st / H - static_cast<double>(starts)) < 0.1 * static_cast<double>(starts) && std::fabs(dur / endings - 0.02) < 0.15 * 0.02 &&
                  std::fabs(st2 / H - static_cast<double>(starts2)) < 0.15 * static_cast<double>(starts2) && std::fabs(dur2 / endings2 - 0.5) < 0.15 * 0.5,
              std::to_string(st / H) + " against " + std::to_string(starts) + ", " + std::to_string(st2 / H) + " against " + std::to_string(starts2) + "; " +
                  std::to_string(dur / endings) + ", " + std::to_string(dur2 / endings2));
    s.require("... the building born at most once, in about 1 − e⁻² of the histories (0.70 to 0.97), and altered only once born",
              born / H > 0.7 && born / H < 0.97 && after_birth, std::to_string(born / H));
    s.require("... every grown foam's cells counted from its stars, χ = Σκ", consistent);
    calx::random::Stream r1(9, 1, 1, 1), r2(9, 1, 1, 1);
    const auto g1 = grow(Fp, Sp, Pg, r1), g2 = grow(Fp, Sp, Pg, r2);
    s.require("... the same stream, the same foam", g1.vertex == g2.vertex && g1.time == g2.time);

    // A state an era never saw: regimes on three streets begin at 3 a year and last 0.2 years in the first ten years only, the growth
    // learned with a break at year 10. A regime grown across year 10 stands in a state the second era never saw, and ends there at
    // the window's rate.
    {
      std::vector<Ev> ev2;
      for (double t = 0;;) {
        t += rs.exponential() / 3.0;
        if (t + 0.2 >= 10) break;
        ev2.push_back({t, "regime", {0, 1, 2}}), ev2.push_back({t + 0.2, "regime ended", {0, 1, 2}});
        t += 0.2;
      }
      mt::Encounters E2;
      E2.participants = 3;
      std::vector<std::string> vk2 = {"start"};
      std::vector<double> t2 = {0};
      E2.present.push_back({0, 1, 2});
      for (const auto& e : ev2) E2.present.push_back(e.who), vk2.push_back(e.kind), t2.push_back(e.t);
      E2.present.push_back({0, 1, 2}), vk2.push_back("end"), t2.push_back(T);
      std::vector<char> ends2(E2.present.size(), 0);
      ends2.front() = ends2.back() = 1;
      const auto F2 = make(E2, vk2, t2, ends2, std::vector<std::string>(3, "street"), {});
      const auto S2 = stars(F2);
      const auto P2 = learn_growth(F2, S2, {}, {10.0});
      int crossed = 0, left = 0, late = 0;
      for (int h = 0; h < H; h++) {
        const auto Fg = grow(F2, S2, P2, rs);
        double began = -1;
        bool open = false;
        for (std::size_t v = 1; v + 1 < Fg.size(); v++) {
          if (Fg.vertex[v] == "regime") open = true, began = Fg.time[v], late += Fg.time[v] >= 10;
          if (Fg.vertex[v] == "regime ended") {
            crossed += began < 10 && Fg.time[v] >= 10;
            open = false;
          }
        }
        left += open;
      }
      s.require("a state an era never saw: of 100 histories, a regime grown across the era's start in about a third (15 to 60), each ended "
                "at the window's rate; none begun after it, none left in force",
                crossed >= 15 && crossed <= 60 && left == 0 && late == 0,
                std::to_string(crossed) + " crossed, " + std::to_string(left) + " left in force, " + std::to_string(late) + " begun late");
    }
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
