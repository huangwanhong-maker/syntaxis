// The lives of a foam's participants as words, and their grammar (docs/theory.md, §7.6).
//
// A participant's worldline passes through the events it attends; read in time, the rules of those events, each with the
// participant's role in it, form a word: the participant's life. Events on one worldline are ordered in the record's trace, so a life
// is the trace's projection on one participant, and the record's trace is determined by its projections together with the stars
// that join them. The words of all lives are compressed by Re-Pair (morphe's `repair`): each rule of the straight-line grammar is a
// pattern of consecutive events that recurs in lives, nested; its uses count how often it is spent in deriving them. Events whose kind
// is among those left out and which only pass (every participant through, no corner opened or closed: an observation, a tick of the
// clock) may be left out of the words, to read the lives of changes alone.
//
// Symbols are numbered in the order of their keys, so that a grammar does not depend on the order of the record's participants.

#pragma once

#include <morphe/formal/repair.hpp>

#include <syntaxis/foams/foams.hpp>

#include <algorithm>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace syntaxis::foams {

// A readable account of a star's action (its rule) or of the whole star: its kind; its participants by kind and role; the corners it
// opens and closes; with `full`, the passes and the states too. Deterministic, not canonical: the keys are the canonical forms.
inline std::string describe(const Foam& F, const Star& s, bool full = false) {
  std::map<std::string, int> who, corner, passing;
  for (std::size_t i = 0; i < s.participants.size(); i++) {
    const auto x = s.participants[i];
    std::string w = F.participant[x];
    if (full) {
      const auto st = F.state(x, s.vertex);
      if (!st.empty()) w += " (" + st + ")";
    }
    who[w + " " + role(s.in[i], s.out[i])]++;
  }
  for (const auto& c : s.corners) {
    const auto& label = F.pair[F.W.lenses[c.lens].pair];
    if (c.kind == Corner::pass) passing[label]++;
    else corner[std::string(c.kind == Corner::open ? "opens " : "closes ") + label]++;
  }
  const auto list = [](const std::map<std::string, int>& m) {
    std::string out;
    for (const auto& [k, n] : m) out += (out.empty() ? "" : ", ") + k + (n > 1 ? " ×" + std::to_string(n) : "");
    return out;
  };
  std::string out = F.vertex[s.vertex] + ": " + list(who);
  if (!corner.empty()) out += "; " + list(corner);
  if (full && !passing.empty()) out += "; passes " + list(passing);
  return out;
}

struct Lives {
  std::vector<std::string> symbols;                 // symbol ↦ its key (the rule, then the participant's place in it)
  std::vector<std::string> names;                   // symbol ↦ a readable account
  std::vector<std::size_t> owner;                   // word ↦ its participant
  std::vector<std::vector<std::uint32_t>> words;    // the lives, as given to Re-Pair
  morphe::formal::StraightLine grammar;
  std::vector<std::size_t> uses;                    // nonterminal ↦ how often the derivation of the words spends it
  std::size_t left_out = 0;                         // events left out of the words (identity events of the kinds asked)
};

// Whether an event only passes: every participant through and no corner opened or closed.
inline bool passes_only(const Star& s) {
  if (s.opens || s.closes) return false;
  for (std::size_t i = 0; i < s.participants.size(); i++)
    if (!(s.in[i] && s.out[i])) return false;
  return true;
}

// A participant's place in a star (the i-th present): its kind, its role and the corners it opens or closes, with the partner's kind.
inline std::string place_of(const Foam& F, const Star& s, std::size_t i) {
  std::vector<std::string> cs;
  for (const auto& c : s.corners) {
    if (c.kind == Corner::pass || (c.a != i && c.b != i)) continue;
    const std::size_t other = c.a == i ? c.b : c.a;
    cs.push_back(std::string(c.kind == Corner::open ? "opens " : "closes ") + F.pair[F.W.lenses[c.lens].pair] + " with " + F.participant[s.participants[other]]);
  }
  std::sort(cs.begin(), cs.end());
  std::string out = F.participant[s.participants[i]] + " " + role(s.in[i], s.out[i]);
  for (const auto& c : cs) out += ", " + c;
  return out;
}

// The lives of every participant, or of those of one kind (`only`, empty for all).
inline Lives lives(const Foam& F, const std::vector<Star>& S, const std::set<std::string>& left_out_kinds = {}, const std::string& only = {}) {
  Lives L;
  const auto place = [&](const Star& s, std::size_t i) { return place_of(F, s, i); };
  // Collect the symbols first, to number them in the order of their keys.
  std::vector<std::vector<std::pair<std::string, std::string>>> raw(F.participant.size());  // participant ↦ (key, name) in time
  for (std::size_t x = 0; x < F.participant.size(); x++)
    for (auto v : F.W.worldline[x]) {
      if (!only.empty() && F.participant[x] != only) break;
      const Star& s = S[v];
      if (s.boundary) continue;
      if (left_out_kinds.count(F.vertex[v]) && passes_only(s)) {
        L.left_out++;
        continue;
      }
      std::size_t i = 0;
      while (i < s.participants.size() && s.participants[i] != x) i++;
      const std::string where = place(s, i);
      raw[x].push_back({s.rule + "\x1f" + where, describe(F, s) + " [as " + where + "]"});
    }
  std::map<std::string, std::string> named;
  for (const auto& r : raw)
    for (const auto& [k, nm] : r) named.emplace(k, nm);
  std::map<std::string, std::uint32_t> id;
  for (const auto& [k, nm] : named) {
    id[k] = static_cast<std::uint32_t>(L.symbols.size());
    L.symbols.push_back(k), L.names.push_back(nm);
  }
  for (std::size_t x = 0; x < raw.size(); x++) {
    if (raw[x].empty()) continue;
    L.owner.push_back(x);
    L.words.emplace_back();
    for (const auto& [k, nm] : raw[x]) L.words.back().push_back(id.at(k));
  }
  L.grammar = morphe::formal::re_pair(L.words);
  // How often each nonterminal is spent: its occurrences in the words, then through every rule that holds it.
  const auto& G = L.grammar;
  L.uses.assign(G.rules.size(), 0);
  for (const auto& w : G.words)
    for (auto a : w)
      if (G.nonterminal(a)) L.uses[a - G.first_new]++;
  for (std::size_t k = G.rules.size(); k-- > 0;)
    for (auto a : {G.rules[k].first, G.rules[k].second})
      if (G.nonterminal(a)) L.uses[a - G.first_new] += L.uses[k];
  return L;
}

}  // namespace syntaxis::foams
