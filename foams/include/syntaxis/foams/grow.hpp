// Growing foams from their grammar (docs/theory.md, §7.8).
//
// A foam is grown one event at a time by the event operator: at each step a site is chosen, the grammar chooses the event's rule
// among the rules that fit the site, and the event's vertex is attached there with its star. A site is where an event can happen:
//   - generic: any one participant of a kind, for an event that touches one participant passing through and opens or closes no
//     lens (an alteration of any building);
//   - concrete: the very participants and pairs of an event the record holds, and the participants and pairs it begins (a regime on
//     one set of street spaces, a building born against its neighbours); a site that begins a participant happens once.
// A site's rate in each condition (the states its participants carry, distinct and sorted) and each era of the window is the
// number of the record's events there in that condition and era over the time it spent so: rates per exposure, at the sites the
// record shows. The rule is drawn from the grammar's context-sensitive probability given the context the event would have (the
// match, the states carried into it, the lenses open beside it), renormalised over the rules whose action fits the site. The states
// are attributes the productions update: the state a participant carries out of an event is drawn from the record's own table of
// what followed that state, rule and place, and kept where the table is silent (Knuth's attributes, here learned). The process is a
// continuous-time Markov chain on the attributed network of the moment, simulated by Gillespie's direct method ("Exact stochastic
// simulation of coupled chemical reactions", J. Phys. Chem. 81, 2340–2361, 1977) on a calx stream; an era's end restarts it, the
// chain being memoryless.
//
// What it cannot do: grow a site the record never shows for more than one participant (a permit on a set of streets no permit
// used, a building born where none was); hold time beyond the rates of its eras; change a state through an event it left out.

#pragma once

#include <syntaxis/foams/foams.hpp>
#include <syntaxis/foams/grammar.hpp>
#include <syntaxis/foams/lives.hpp>

#include <calx/random/stream.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace syntaxis::foams {

// The match and context keys of an event from its consumed participants (each with its kind, the state carried in, and the kinds
// of the lenses that pass beside it) and the lenses it closes among them (their positions and the pair's kind): built as stars()
// builds them, so that an event not yet in a foam is keyed as it would be once there.
inline std::pair<std::string, std::string> match_and_context(const std::vector<std::string>& kinds, const std::vector<std::string>& states,
                                                             const std::vector<std::vector<std::string>>& passes,
                                                             const std::vector<std::tuple<std::size_t, std::size_t, std::string>>& closes) {
  mc::LabelledGraph match, context;
  match.vertex.push_back("m"), context.vertex.push_back("m");
  for (std::size_t i = 0; i < kinds.size(); i++) {
    match.vertex.push_back("x:" + kinds[i]);
    context.vertex.push_back("x:" + kinds[i] + (states[i].empty() ? "" : "|" + states[i]));
    match.edges.push_back({0, i + 1, "at"}), context.edges.push_back({0, i + 1, "at"});
  }
  for (const auto& [a, b, pair] : closes) {
    const auto label = detail::corner_label(Corner::close, pair);
    match.edges.push_back({a + 1, b + 1, label}), context.edges.push_back({a + 1, b + 1, label});
  }
  for (std::size_t i = 0; i < passes.size(); i++)
    for (const auto& pair : passes[i]) context.edges.push_back({i + 1, i + 1, detail::corner_label(Corner::pass, pair)});
  return {detail::key_any(match).code, detail::key_any(context).code};
}

struct GrowSite {
  bool generic = false;          // any one participant of `kind`, passing through, opening and closing no lens
  std::string kind;              // a generic site's participants' kind
  std::size_t example = 0;       // an event of the record at the site: its participants, roles and pairs
  std::string shape;             // its action without its kind
  bool once = false;             // it begins a participant
  std::size_t events = 0;        // the record's events there
  std::vector<std::map<std::string, double>> rate;  // era ↦ condition ↦ rate per unit of time (per participant, for a generic site)
};

struct Growth {
  double from = 0, to = 0;                           // the window
  std::vector<double> breaks;                        // the eras' boundaries inside the window, increasing
  std::vector<GrowSite> sites;
  std::map<std::string, std::vector<std::string>> rules_of;          // shape ↦ the rules with that action
  std::map<std::string, std::size_t> example;                        // rule ↦ an event of it
  std::map<std::tuple<std::string, std::string, std::string>, std::map<std::string, std::size_t>> after;  // (state in, rule, place) ↦ state out ↦ count
  std::set<std::string> left_out;                    // kinds of event the growth leaves out
  std::size_t era(double t) const { return static_cast<std::size_t>(std::upper_bound(breaks.begin(), breaks.end(), t) - breaks.begin()); }
};

namespace detail {
inline bool generic_star(const Star& s) { return s.participants.size() == 1 && s.in[0] && s.out[0] && s.opens == 0 && s.closes == 0; }
inline std::string condition(const std::vector<std::string>& states) {
  std::set<std::string> distinct(states.begin(), states.end());
  std::string out;
  for (const auto& st : distinct) out += (out.empty() ? "" : "\x1f") + st;
  return out;
}
inline std::string most(const std::map<std::string, std::size_t>& m, const std::string& fallback) {
  std::string best = fallback;
  std::size_t n = 0;
  for (const auto& [k, c] : m)
    if (c > n) best = k, n = c;
  return best;
}
}  // namespace detail

// What a foam's growth needs, read from the foam and its stars: its sites with their rates by condition and era (the eras cut at
// `breaks`), the rules each action admits, and the table of the states that followed. Events of the kinds left out that only pass are
// not grown; they are taken to change no state.
inline Growth learn_growth(const Foam& F, const std::vector<Star>& S, const std::set<std::string>& left_out_kinds = {}, std::vector<double> breaks = {}) {
  detail::require(F.size() >= 2, "a foam with its two ends");
  Growth P;
  P.from = F.time.front(), P.to = F.time.back(), P.left_out = left_out_kinds;
  std::sort(breaks.begin(), breaks.end());
  for (auto b : breaks)
    if (b > P.from && b < P.to) P.breaks.push_back(b);
  const std::size_t eras = P.breaks.size() + 1;
  const auto era_start = [&](std::size_t k) { return k == 0 ? P.from : P.breaks[k - 1]; };
  const auto era_end = [&](std::size_t k) { return k + 1 == eras ? P.to : P.breaks[k]; };
  std::vector<char> kept(F.size(), 0);
  for (std::size_t v = 0; v < F.size(); v++)
    kept[v] = !S[v].boundary && !(left_out_kinds.count(F.vertex[v]) && passes_only(S[v]));
  // The rules each action admits, and the table of what followed.
  for (std::size_t v = 0; v < F.size(); v++) {
    if (!kept[v]) continue;
    auto& rs = P.rules_of[S[v].shape];
    if (std::find(rs.begin(), rs.end(), S[v].rule) == rs.end()) rs.push_back(S[v].rule);
    P.example.emplace(S[v].rule, v);
  }
  for (auto& [shape, rs] : P.rules_of) std::sort(rs.begin(), rs.end());
  const std::size_t X = F.participant.size();
  std::vector<std::vector<std::size_t>> kept_on(X);  // each participant's kept events, in time
  for (std::size_t x = 0; x < X; x++)
    for (auto v : F.W.worldline[x])
      if (kept[v]) kept_on[x].push_back(v);
  const auto position = [&](const Star& s, std::size_t x) {
    return static_cast<std::size_t>(std::find(s.participants.begin(), s.participants.end(), x) - s.participants.begin());
  };
  for (std::size_t x = 0; x < X; x++)
    for (std::size_t k = 0; k + 1 < kept_on[x].size(); k++) {
      const auto v = kept_on[x][k], w = kept_on[x][k + 1];
      P.after[{F.state(x, v), S[v].rule, place_of(F, S[v], position(S[v], x))}][F.state(x, w)]++;
    }
  // Each participant's states in time: [its life's start, its first kept event] the state carried into that event, and so on; after
  // its last kept event, the state the table gives most often.
  struct Piece {
    double from, to;
    std::string state;
  };
  std::vector<std::vector<Piece>> pieces(X);
  std::vector<std::pair<double, double>> life(X, {0, -1});
  for (std::size_t x = 0; x < X; x++) {
    const auto& line = F.W.worldline[x];
    if (line.empty()) continue;
    life[x] = {F.time[line.front()], F.time[line.back()]};
    double t = life[x].first;
    std::string carried = kept_on[x].empty() ? F.state(x, line.back()) : F.state(x, kept_on[x].front());
    for (std::size_t k = 0; k < kept_on[x].size(); k++) {
      const auto v = kept_on[x][k];
      carried = F.state(x, v);
      pieces[x].push_back({t, F.time[v], carried});
      t = F.time[v];
      const auto it = P.after.find({carried, S[v].rule, place_of(F, S[v], position(S[v], x))});
      carried = it == P.after.end() ? carried : detail::most(it->second, carried);
    }
    pieces[x].push_back({t, life[x].second, carried});
  }
  const auto state_at = [&](std::size_t x, double t) {
    for (const auto& p : pieces[x])
      if (t <= p.to) return p.state;
    return pieces[x].empty() ? std::string() : pieces[x].back().state;
  };
  // The sites.
  std::map<std::string, std::size_t> site_of;
  std::vector<std::vector<std::size_t>> consumed_of;  // a concrete site's consumed participants
  for (std::size_t v = 0; v < F.size(); v++) {
    if (!kept[v]) continue;
    const Star& s = S[v];
    std::string key;
    std::vector<std::size_t> consumed;
    bool once = false;
    if (detail::generic_star(s)) {
      key = "generic\x1e" + F.participant[s.participants[0]] + "\x1e" + s.shape;
    } else {
      std::vector<std::size_t> begun, pairs;
      for (std::size_t i = 0; i < s.participants.size(); i++) (s.in[i] ? consumed : begun).push_back(s.participants[i]);
      for (const auto& c : s.corners)
        if (c.kind != Corner::pass) pairs.push_back(F.W.lenses[c.lens].pair);
      std::sort(consumed.begin(), consumed.end()), std::sort(begun.begin(), begun.end()), std::sort(pairs.begin(), pairs.end());
      pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
      once = !begun.empty();
      key = "concrete\x1e" + s.shape;
      for (auto x : consumed) key += "\x1e" + std::to_string(x);
      key += "\x1e|";
      for (auto x : begun) key += "\x1e" + std::to_string(x);
      key += "\x1e|";
      for (auto p : pairs) key += "\x1e" + std::to_string(p);
    }
    auto it = site_of.find(key);
    if (it == site_of.end()) {
      GrowSite g;
      g.generic = detail::generic_star(s), g.kind = g.generic ? F.participant[s.participants[0]] : "", g.example = v, g.shape = s.shape, g.once = once;
      g.rate.assign(eras, {});
      it = site_of.emplace(key, P.sites.size()).first;
      P.sites.push_back(g);
      consumed_of.push_back(consumed);
    }
    auto& g = P.sites[it->second];
    g.events++;
    std::vector<std::string> st;
    for (std::size_t i = 0; i < s.participants.size(); i++)
      if (s.in[i]) st.push_back(F.state(s.participants[i], v));
    g.rate[P.era(F.time[v])][detail::condition(st)] += 1;  // a count, divided by the exposure below
  }
  // Exposures: a generic site, every participant of its kind alive in each state; a concrete one, the time all its consumed
  // participants are alive (until its event, for a site that happens once) in each condition.
  for (std::size_t i = 0; i < P.sites.size(); i++) {
    auto& g = P.sites[i];
    std::vector<std::map<std::string, double>> exposure(eras);
    const auto add = [&](double a, double b, const std::string& c) {
      for (std::size_t k = 0; k < eras; k++) {
        const double lo = std::max(a, era_start(k)), hi = std::min(b, era_end(k));
        if (hi > lo) exposure[k][c] += hi - lo;
      }
    };
    if (g.generic) {
      for (std::size_t x = 0; x < X; x++)
        if (F.participant[x] == g.kind)
          for (const auto& p : pieces[x]) add(p.from, p.to, p.state);
    } else {
      const auto& cs = consumed_of[i];
      double lo = P.from, hi = P.to;
      for (auto x : cs) lo = std::max(lo, life[x].first), hi = std::min(hi, life[x].second);
      if (g.once) hi = std::min(hi, F.time[g.example]);
      std::set<double> cuts = {lo, hi};
      for (auto x : cs)
        for (const auto& p : pieces[x])
          if (p.to > lo && p.to < hi) cuts.insert(p.to);
      for (auto k : P.breaks)
        if (k > lo && k < hi) cuts.insert(k);
      for (auto a = cuts.begin(), b = std::next(cuts.begin()); b != cuts.end(); ++a, ++b) {
        std::vector<std::string> st;
        const double mid = (*a + *b) / 2;
        for (auto x : cs) st.push_back(state_at(x, mid));
        add(*a, *b, detail::condition(st));
      }
    }
    for (std::size_t k = 0; k < eras; k++)
      for (auto& [c, n] : g.rate[k]) {
        const auto e = exposure[k].find(c);
        n = (e == exposure[k].end() || e->second <= 0) ? 0 : n / e->second;
      }
  }
  return P;
}

// A foam grown from the grammar G and the growth P, from the record F's moment at its start (its participants then alive with the
// states they carried into their first events, its pairs then open) to its end. Participants and pairs keep their records' indices
// and kinds; those never begun are left out.
inline Foam grow(const Foam& F, const std::vector<Star>& S, const Grammar& G, const Growth& P, calx::random::Stream& rs) {
  const std::size_t X = F.participant.size(), Q = F.encounters.pairs.size();
  std::vector<char> alive(X, 0), pair_alive(Q, 0), fired(P.sites.size(), 0), begun(X, 0);
  std::vector<std::string> state(X);
  for (auto x : F.encounters.present.front()) {
    alive[x] = begun[x] = 1;
    std::size_t first = F.size();
    for (auto v : F.W.worldline[x])
      if (!S[v].boundary && !(P.left_out.count(F.vertex[v]) && passes_only(S[v]))) {
        first = v;
        break;
      }
    state[x] = first < F.size() ? F.state(x, first) : F.state(x, F.W.worldline[x].back());
  }
  for (std::size_t q = 0; q < Q; q++)
    if (!F.encounters.pairs[q].meetings.empty() && F.encounters.pairs[q].meetings.front() == 0) pair_alive[q] = 1;
  // The foam being grown.
  std::vector<std::vector<std::size_t>> present = {F.encounters.present.front()};
  std::vector<std::string> kind = {F.vertex.front()};
  std::vector<double> time = {P.from};
  std::vector<std::vector<std::size_t>> meetings(Q);
  std::vector<std::map<std::size_t, std::string>> states(X);
  for (std::size_t q = 0; q < Q; q++)
    if (pair_alive[q]) meetings[q].push_back(0);
  std::vector<std::vector<std::size_t>> pairs_of(X);
  for (std::size_t q = 0; q < Q; q++) pairs_of[F.encounters.pairs[q].a].push_back(q), pairs_of[F.encounters.pairs[q].b].push_back(q);
  double t = P.from;
  std::size_t era = 0;
  const std::size_t eras = P.breaks.size() + 1;
  const auto era_end = [&](std::size_t k) { return k + 1 == eras ? P.to : P.breaks[k]; };
  while (era < eras) {
    // The weight of every site now, and for a generic one the participants in each condition.
    std::vector<double> weight(P.sites.size(), 0);
    std::vector<std::map<std::string, std::vector<std::size_t>>> who(P.sites.size());
    double total = 0;
    for (std::size_t i = 0; i < P.sites.size(); i++) {
      const auto& g = P.sites[i];
      const auto& rates = g.rate[era];
      if (rates.empty()) continue;
      if (g.generic) {
        for (std::size_t x = 0; x < X; x++)
          if (alive[x] && F.participant[x] == g.kind) {
            const auto r = rates.find(state[x]);
            if (r != rates.end() && r->second > 0) who[i][state[x]].push_back(x), weight[i] += r->second;
          }
      } else {
        if (g.once && fired[i]) continue;
        const Star& s = S[g.example];
        bool ok = true;
        std::vector<std::string> st;
        for (std::size_t k = 0; k < s.participants.size() && ok; k++) {
          const auto x = s.participants[k];
          if (s.in[k]) ok = alive[x], st.push_back(state[x]);
          else ok = !begun[x];
        }
        if (!ok) continue;
        const auto r = rates.find(detail::condition(st));
        if (r != rates.end()) weight[i] = r->second;
      }
      total += weight[i];
    }
    const double next = total > 0 ? t + rs.exponential() / total : std::numeric_limits<double>::infinity();
    if (next >= era_end(era)) {
      t = era_end(era), era++;
      continue;
    }
    t = next;
    // The site, by weight; for a generic site, a participant in proportion to its condition's rate.
    double u = rs.uniform() * total;
    std::size_t i = 0;
    while (i + 1 < P.sites.size() && u >= weight[i]) u -= weight[i], i++;
    const auto& g = P.sites[i];
    const Star& tpl = S[g.example];
    std::vector<std::size_t> consumed, newcomers;
    std::vector<std::string> places;
    if (g.generic) {
      double w = rs.uniform() * weight[i];
      std::size_t x = 0;
      for (const auto& [c, xs] : who[i]) {
        const double r = g.rate[era].at(c);
        if (w < r * static_cast<double>(xs.size())) {
          x = xs[std::min(xs.size() - 1, static_cast<std::size_t>(w / r))];
          break;
        }
        w -= r * static_cast<double>(xs.size()), x = xs.back();
      }
      consumed = {x};
      places = {place_of(F, tpl, 0)};
    } else {
      for (std::size_t k = 0; k < tpl.participants.size(); k++) (tpl.in[k] ? consumed : newcomers).push_back(tpl.participants[k]);
    }
    // The pairs the event meets (the example's, alive), opens (those it begins) and ends (those it meets for the last time in the
    // record); the context the event would have.
    std::set<std::size_t> meets, opens, closing;
    if (!g.generic)
      for (const auto& c : tpl.corners) {
        if (c.kind == Corner::pass) continue;
        const auto q = F.W.lenses[c.lens].pair;
        (pair_alive[q] ? meets : opens).insert(q);
        if (F.encounters.pairs[q].meetings.back() == g.example) closing.insert(q);
      }
    std::vector<std::string> kinds, carried;
    std::vector<std::vector<std::string>> passes;
    for (auto x : consumed) {
      kinds.push_back(F.participant[x]), carried.push_back(state[x]), passes.emplace_back();
      for (auto q : pairs_of[x])
        if (pair_alive[q] && !meets.count(q)) passes.back().push_back(F.pair[q]);
    }
    std::vector<std::tuple<std::size_t, std::size_t, std::string>> closes;
    for (auto q : meets) {
      const auto a = static_cast<std::size_t>(std::find(consumed.begin(), consumed.end(), F.encounters.pairs[q].a) - consumed.begin());
      const auto b = static_cast<std::size_t>(std::find(consumed.begin(), consumed.end(), F.encounters.pairs[q].b) - consumed.begin());
      if (a < consumed.size() && b < consumed.size()) closes.emplace_back(a, b, F.pair[q]);
    }
    const auto [m, c] = match_and_context(kinds, carried, passes, closes);
    // The rule, among those of the site's action, by the grammar's probability in that context.
    const auto& candidates = P.rules_of.at(g.shape);
    std::vector<double> p(candidates.size());
    double sum = 0;
    for (std::size_t k = 0; k < candidates.size(); k++) sum += p[k] = p2(G, candidates[k], m, c);
    double w = rs.uniform() * sum;
    std::size_t pick = 0;
    while (pick + 1 < candidates.size() && w >= p[pick]) w -= p[pick], pick++;
    const auto& rule = candidates[pick];
    // The event attached.
    const std::size_t v = present.size();
    std::vector<std::size_t> here = consumed;
    here.insert(here.end(), newcomers.begin(), newcomers.end());
    present.push_back(here), kind.push_back(F.vertex[P.example.at(rule)]), time.push_back(t);
    for (auto x : consumed) states[x][v] = state[x];
    for (auto q : meets) meetings[q].push_back(v);
    for (auto q : opens) meetings[q].push_back(v), pair_alive[q] = 1;
    for (auto q : closing) pair_alive[q] = 0;
    for (auto x : newcomers) alive[x] = begun[x] = 1;
    // What each participant carries out: the table's draw, or what it carried in; one that ends is gone, with its pairs.
    for (std::size_t k = 0; k < here.size(); k++) {
      const auto x = here[k];
      const std::size_t at = g.generic ? 0 : static_cast<std::size_t>(std::find(tpl.participants.begin(), tpl.participants.end(), x) - tpl.participants.begin());
      if (!g.generic && tpl.in[at] && !tpl.out[at]) {
        alive[x] = 0;
        for (auto q : pairs_of[x]) pair_alive[q] = 0;
        continue;
      }
      const auto it = P.after.find({state[x], rule, g.generic ? places[0] : place_of(F, tpl, at)});
      if (it == P.after.end()) continue;
      std::size_t n = 0;
      for (const auto& [st, cnt] : it->second) n += cnt;
      auto r = rs.below(n);
      for (const auto& [st, cnt] : it->second) {
        if (r < cnt) {
          state[x] = st;
          break;
        }
        r -= cnt;
      }
    }
    if (g.once) fired[i] = 1;
  }
  // The window's end: every participant alive, every pair open.
  const std::size_t end = present.size();
  std::vector<std::size_t> last;
  for (std::size_t x = 0; x < X; x++)
    if (alive[x]) last.push_back(x);
  present.push_back(last), kind.push_back(F.vertex.back()), time.push_back(P.to);
  for (std::size_t q = 0; q < Q; q++)
    if (pair_alive[q]) meetings[q].push_back(end);
  // Participants and pairs never begun left out; the rest renumbered in their records' order.
  std::vector<std::size_t> index(X, 0);
  std::vector<std::string> participant;
  std::size_t n = 0;
  for (std::size_t x = 0; x < X; x++)
    if (begun[x]) index[x] = n++, participant.push_back(F.participant[x]);
  mt::Encounters E;
  E.participants = n;
  for (auto& pr : present) {
    std::vector<std::size_t> q;
    for (auto x : pr) q.push_back(index[x]);
    E.present.push_back(q);
  }
  std::vector<std::string> pair;
  for (std::size_t q = 0; q < Q; q++)
    if (meetings[q].size() >= 2) E.pairs.push_back({index[F.encounters.pairs[q].a], index[F.encounters.pairs[q].b], meetings[q]}), pair.push_back(F.pair[q]);
  std::vector<char> boundary(present.size(), 0);
  boundary.front() = boundary.back() = 1;
  auto out = make(std::move(E), std::move(kind), std::move(time), std::move(boundary), std::move(participant), std::move(pair));
  out.states.assign(n, {});
  for (std::size_t x = 0; x < X; x++)
    if (begun[x]) out.states[index[x]] = states[x];
  return out;
}

}  // namespace syntaxis::foams
