// The process net of a foam, and its compression by digrams: patterns of events that span several worldlines (docs/theory.md, §7.7).
//
// A foam's events are boxes and its segments wires. Each event is a box labelled by its rule, with an in-port for each participant
// it consumes (an arriving segment) and an out-port for each it leaves (a leaving segment); each segment is a wire from the out-port
// of the event it leaves to the in-port of the next event on the same worldline. A port is named by its signature: the participant's
// kind and place in its event. Events of chosen kinds that only pass (every participant through, no corner) may be left out first,
// their wires joined across them: such an event is a refinement of the worldlines it cuts (Proposition 7.1's remark). The net is the
// record's trace as a process: its wires are the dependences between events.
//
// A digram is two boxes joined by at least one wire, the first the one whose first event is earlier. Its type is the pair of labels
// and the multiset of its wires, each by the direction it runs and the signatures of the ports it joins. The type with the most
// occurrences that share no box (chosen greedily in the order of time), two at least, becomes a new label, a pattern, whose ports are
// the two boxes' ports not wired between them, each signature prefixed by the part it comes from; every chosen occurrence is replaced
// by one box of the pattern; and again, until no type occurs twice. This is Re-Pair (Larsson and Moffat 2000) on the net, as gRePair
// compresses hypergraphs (Maneth and Peternek 2018, "Grammar-based graph compression", Information Systems 76, 19–45). A pattern is
// a part of the record that recurs; its instances bind the worldlines that run through them, so that patterns spanning several
// worldlines are found from the foam itself, not from the word of one worldline (lives.hpp).
//
// A type is exact up to the symmetries of its boxes: two occurrences of a type are isomorphic as wired pairs when the ports of one
// signature are interchangeable within their box, as they are in every star with no corner between two participants and in every
// pattern made of such stars. The order of the result does not depend on the order of the record's participants.

#pragma once

#include <syntaxis/foams/foams.hpp>
#include <syntaxis/foams/lives.hpp>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace syntaxis::foams {

struct NetPort {
  std::uint32_t signature = 0;  // the participant's kind and place in its event; in a pattern, prefixed by the part it comes from
  std::size_t participant = 0;
  bool in = false;              // an arriving segment ends here (else a leaving one begins)
};

struct NetBox {
  std::uint32_t label = 0;                // a terminal (an event's rule) below Net::terminals, else a pattern
  std::vector<std::size_t> events;        // the foam's vertices it covers, increasing
  std::vector<std::size_t> participants;  // the participants of those events, increasing
  std::vector<NetPort> ports;
};

struct NetWire {
  std::size_t from = 0, from_port = 0, to = 0, to_port = 0;  // an out-port of one box to an in-port of another
};

struct NetRule {
  std::uint32_t first = 0, second = 0;  // the parts' labels, the first the one whose first event is earlier
  // The wires between the parts: the direction (0 from the first to the second, 1 back), the out-port's and the in-port's signatures.
  std::vector<std::tuple<int, std::uint32_t, std::uint32_t>> wires;
  std::size_t uses = 0, length = 0, width = 0;     // instances made; events in each; participants in each
  std::vector<std::size_t> example;                // the events of its first instance
  std::vector<std::vector<std::size_t>> spans;     // each instance's participants
};

struct Net {
  std::vector<std::string> labels;     // label ↦ a terminal's rule key, or a pattern's name
  std::vector<std::size_t> example;    // terminal label ↦ an event of that rule
  std::size_t terminals = 0;
  std::vector<std::string> signatures;
  std::vector<NetRule> rules;          // pattern label terminals + i ↦ rules[i]
  std::vector<NetBox> boxes;           // the top level of the derivation
  std::vector<NetWire> wires;          // the wires between them
  std::size_t events = 0;              // the events kept
  std::size_t wires_made = 0;          // the wires of the net before compression
  std::size_t open_ends = 0;           // segments from ∂₀ or to ∂₁, which join no two boxes
  std::size_t left_out = 0;            // events left out
  bool pattern(std::uint32_t label) const { return label >= terminals; }
  const NetRule& rule(std::uint32_t label) const { return rules.at(label - terminals); }
  // The size of the description: before, the events and wires; after, each rule's two parts and wires, the boxes and wires left.
  std::size_t size_before() const { return events + wires_made; }
  std::size_t size_after() const {
    std::size_t s = boxes.size() + wires.size();
    for (const auto& r : rules) s += 2 + r.wires.size();
    return s;
  }
};

namespace detail {
inline std::uint32_t intern(std::vector<std::string>& names, std::map<std::string, std::uint32_t>& ids, const std::string& s) {
  const auto it = ids.find(s);
  if (it != ids.end()) return it->second;
  const auto id = static_cast<std::uint32_t>(names.size());
  names.push_back(s), ids.emplace(s, id);
  return id;
}
}  // namespace detail

// The process net of a foam, events of the kinds left out that only pass taken away.
inline Net net(const Foam& F, const std::vector<Star>& S, const std::set<std::string>& left_out_kinds = {}) {
  constexpr std::size_t none = std::numeric_limits<std::size_t>::max();
  Net N;
  std::map<std::string, std::uint32_t> sig;
  std::vector<char> kept(F.size(), 0);
  for (std::size_t v = 0; v < F.size(); v++) {
    if (S[v].boundary) continue;
    if (left_out_kinds.count(F.vertex[v]) && passes_only(S[v])) {
      N.left_out++;
      continue;
    }
    kept[v] = 1;
  }
  // Terminal labels in the order of their keys, so that they do not depend on the order of the events.
  std::set<std::string> keys;
  for (std::size_t v = 0; v < F.size(); v++)
    if (kept[v]) keys.insert(S[v].rule);
  std::map<std::string, std::uint32_t> lab;
  for (const auto& k : keys) lab[k] = static_cast<std::uint32_t>(N.labels.size()), N.labels.push_back(k), N.example.push_back(none);
  N.terminals = N.labels.size();
  std::vector<std::size_t> box_of(F.size(), none);
  for (std::size_t v = 0; v < F.size(); v++) {
    if (!kept[v]) continue;
    const Star& s = S[v];
    NetBox b;
    b.label = lab.at(s.rule);
    if (N.example[b.label] == none) N.example[b.label] = v;
    b.events = {v};
    b.participants = s.participants;
    std::sort(b.participants.begin(), b.participants.end());
    for (std::size_t i = 0; i < s.participants.size(); i++) {
      const std::string place = place_of(F, s, i);
      if (s.in[i]) b.ports.push_back({detail::intern(N.signatures, sig, place + "|in"), s.participants[i], true});
      if (s.out[i]) b.ports.push_back({detail::intern(N.signatures, sig, place + "|out"), s.participants[i], false});
    }
    box_of[v] = N.boxes.size();
    N.boxes.push_back(std::move(b));
    N.events++;
  }
  const auto port = [&](std::size_t box, std::size_t x, bool in) {
    const auto& ps = N.boxes[box].ports;
    for (std::size_t k = 0; k < ps.size(); k++)
      if (ps[k].participant == x && ps[k].in == in) return k;
    detail::require(false, "a port for every segment's end");
    return none;
  };
  // A wire for each pair of consecutive kept events on a worldline.
  for (std::size_t x = 0; x < F.participant.size(); x++) {
    std::size_t prev = none;
    bool open = false;  // whether the worldline came in from the window's start
    for (auto v : F.W.worldline[x]) {
      if (!kept[v]) {
        if (S[v].boundary && prev == none) open = true;
        continue;
      }
      if (prev != none) N.wires.push_back({box_of[prev], port(box_of[prev], x, false), box_of[v], port(box_of[v], x, true)});
      else if (open) N.open_ends++;
      prev = v;
    }
    if (prev != none && F.W.worldline[x].back() != prev && S[F.W.worldline[x].back()].boundary) N.open_ends++;
  }
  N.wires_made = N.wires.size();
  return N;
}

// The net compressed by digrams until no type occurs twice (or `max_rules` patterns are made, 0 for no bound).
inline void compress(Net& N, std::size_t max_rules = 0) {
  std::vector<char> box_alive(N.boxes.size(), 1), wire_alive(N.wires.size(), 1);
  std::vector<std::vector<std::size_t>> touching(N.boxes.size());
  for (std::size_t w = 0; w < N.wires.size(); w++) touching[N.wires[w].from].push_back(w), touching[N.wires[w].to].push_back(w);
  std::map<std::string, std::uint32_t> sig;
  for (std::size_t k = 0; k < N.signatures.size(); k++) sig.emplace(N.signatures[k], static_cast<std::uint32_t>(k));
  const auto earlier = [&](std::size_t a, std::size_t b) { return N.boxes[a].events.front() < N.boxes[b].events.front(); };
  while (max_rules == 0 || N.rules.size() < max_rules) {
    // The wires between each pair of boxes, the pair ordered by their first events.
    std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> between;
    for (std::size_t w = 0; w < N.wires.size(); w++) {
      if (!wire_alive[w]) continue;
      auto a = N.wires[w].from, b = N.wires[w].to;
      if (!earlier(a, b)) std::swap(a, b);
      between[{a, b}].push_back(w);
    }
    // Each pair's type: the labels and its wires by direction and signatures.
    using Wire = std::tuple<int, std::uint32_t, std::uint32_t>;
    std::map<std::tuple<std::uint32_t, std::uint32_t, std::vector<Wire>>, std::vector<std::pair<std::size_t, std::size_t>>> types;
    for (const auto& [ab, ws] : between) {
      std::vector<Wire> key;
      for (auto w : ws) {
        const auto& W = N.wires[w];
        const int dir = W.from == ab.first ? 0 : 1;
        key.emplace_back(dir, N.boxes[W.from].ports[W.from_port].signature, N.boxes[W.to].ports[W.to_port].signature);
      }
      std::sort(key.begin(), key.end());
      types[{N.boxes[ab.first].label, N.boxes[ab.second].label, key}].push_back(ab);
    }
    // The occurrences of each type that share no box, taken in the order of time; the type with the most.
    const std::vector<std::pair<std::size_t, std::size_t>>* best = nullptr;
    const std::tuple<std::uint32_t, std::uint32_t, std::vector<Wire>>* best_key = nullptr;
    std::vector<std::pair<std::size_t, std::size_t>> chosen, trial;
    for (const auto& [key, occ] : types) {
      if (occ.size() < 2 || (best && occ.size() <= chosen.size())) continue;
      auto sorted = occ;
      std::sort(sorted.begin(), sorted.end(), [&](const auto& p, const auto& q) {
        const auto pa = N.boxes[p.first].events.front(), qa = N.boxes[q.first].events.front();
        if (pa != qa) return pa < qa;
        const auto pb = N.boxes[p.second].events.front(), qb = N.boxes[q.second].events.front();
        return pb != qb ? pb < qb : p < q;
      });
      std::set<std::size_t> used;
      trial.clear();
      for (const auto& [a, b] : sorted)
        if (!used.count(a) && !used.count(b)) trial.push_back({a, b}), used.insert(a), used.insert(b);
      if (trial.size() >= 2 && trial.size() > chosen.size()) chosen = trial, best = &occ, best_key = &key;
    }
    if (!best) break;
    // A new pattern, each chosen occurrence contracted to one box.
    const auto label = static_cast<std::uint32_t>(N.labels.size());
    N.labels.push_back("pattern " + std::to_string(N.rules.size() + 1));
    NetRule R;
    R.first = std::get<0>(*best_key), R.second = std::get<1>(*best_key), R.wires = std::get<2>(*best_key);
    for (const auto& [a, b] : chosen) {
      NetBox c;
      c.label = label;
      std::merge(N.boxes[a].events.begin(), N.boxes[a].events.end(), N.boxes[b].events.begin(), N.boxes[b].events.end(), std::back_inserter(c.events));
      std::set_union(N.boxes[a].participants.begin(), N.boxes[a].participants.end(), N.boxes[b].participants.begin(), N.boxes[b].participants.end(),
                     std::back_inserter(c.participants));
      // The ports wired between the two parts are inside; the rest are the pattern's, prefixed by their part.
      std::set<std::pair<std::size_t, std::size_t>> inside;  // (box, port)
      for (auto w : between.at({a, b})) {
        inside.insert({N.wires[w].from, N.wires[w].from_port}), inside.insert({N.wires[w].to, N.wires[w].to_port});
        wire_alive[w] = 0;
      }
      std::map<std::pair<std::size_t, std::size_t>, std::size_t> moved;  // (old box, old port) ↦ new port
      for (const auto& [part, prefix] : {std::pair{a, std::string("1·")}, std::pair{b, std::string("2·")}}) {
        for (std::size_t k = 0; k < N.boxes[part].ports.size(); k++) {
          if (inside.count({part, k})) continue;
          auto p = N.boxes[part].ports[k];
          p.signature = detail::intern(N.signatures, sig, prefix + N.signatures[p.signature]);
          moved[{part, k}] = c.ports.size();
          c.ports.push_back(p);
        }
      }
      const std::size_t ci = N.boxes.size();
      N.boxes.push_back(std::move(c));
      box_alive.push_back(1), box_alive[a] = 0, box_alive[b] = 0;
      touching.emplace_back();
      for (auto part : {a, b})
        for (auto w : touching[part]) {
          if (!wire_alive[w]) continue;
          auto& W = N.wires[w];
          if (W.from == part) W.from_port = moved.at({part, W.from_port}), W.from = ci;
          if (W.to == part) W.to_port = moved.at({part, W.to_port}), W.to = ci;
          touching[ci].push_back(w);
        }
      R.uses++;
      R.spans.push_back(N.boxes[ci].participants);
      if (R.example.empty()) R.example = N.boxes[ci].events, R.length = N.boxes[ci].events.size(), R.width = N.boxes[ci].participants.size();
    }
    N.rules.push_back(std::move(R));
  }
  // The top level: the boxes and wires left, renumbered.
  std::vector<std::size_t> index(N.boxes.size(), 0);
  std::vector<NetBox> boxes;
  for (std::size_t i = 0; i < N.boxes.size(); i++)
    if (box_alive[i]) index[i] = boxes.size(), boxes.push_back(std::move(N.boxes[i]));
  std::vector<NetWire> wires;
  for (std::size_t w = 0; w < N.wires.size(); w++)
    if (wire_alive[w]) wires.push_back({index[N.wires[w].from], N.wires[w].from_port, index[N.wires[w].to], N.wires[w].to_port});
  N.boxes = std::move(boxes), N.wires = std::move(wires);
}

// The net of a foam, compressed.
inline Net compressed(const Foam& F, const std::vector<Star>& S, const std::set<std::string>& left_out_kinds = {}, std::size_t max_rules = 0) {
  auto N = net(F, S, left_out_kinds);
  compress(N, max_rules);
  return N;
}

// The terminal labels a label derives, its first part's then its second's.
inline std::vector<std::uint32_t> derive(const Net& N, std::uint32_t label) {
  if (!N.pattern(label)) return {label};
  auto out = derive(N, N.rule(label).first);
  const auto second = derive(N, N.rule(label).second);
  out.insert(out.end(), second.begin(), second.end());
  return out;
}

// How many wires each pattern holds inside it, its parts' and its own.
inline std::size_t inner_wires(const Net& N, std::uint32_t label) {
  if (!N.pattern(label)) return 0;
  const auto& r = N.rule(label);
  return r.wires.size() + inner_wires(N, r.first) + inner_wires(N, r.second);
}

}  // namespace syntaxis::foams
