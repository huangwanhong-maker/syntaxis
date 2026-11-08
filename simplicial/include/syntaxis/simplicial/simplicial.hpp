// Grammars whose words are 2-simplices (docs/theory.md, §6).
//
// A triangle is the lowest simplex with a shape of its own: two of its angles determine it up to similarity. Its *letter* is its
// shape class, its angles sorted from the largest and taken to a resolution in degrees (a triangle and its mirror image share a
// letter, and their handedness is kept apart where it matters). A 2-complex is then a text in these letters, and two grammars read it:
//
// - **Gluing.** A complex is grown from one triangle by gluing triangles along edges, as a sentence is grown by concatenation; read
//   back, a derivation is a spanning tree of the complex's dual graph (breadth first from a root). Each triangle's edges have slots by
//   length (0 the longest, opposite the largest angle; 1; 2 the shortest). The grammar is a probabilistic tree grammar: a nonterminal is
//   a letter entered through a slot, and its production glues to each of its other slots either nothing or a triangle of some letter
//   by one of its slots and a handedness, each slot independently, with the probabilities of their counts in the derivation. Sampling
//   grows a complex by drawing these choices and placing each triangle on its parent's edge by its angles; a triangle that would
//   overlap one already placed is refused (a context condition: what may be glued depends on what is already there). The
//   letters' entropy and their conditional entropy given the parent's letter and slot measure how much the gluing is constrained.
// - **Replacement.** A triangle rewrites into a small triangulation of itself: Pachner's 1–3 move (a vertex at the centroid, three
//   triangles; Pachner 1991, "P.L. homeomorphic manifolds are equivalent by elementary shellings", Europ. J. Combin. 12), the midpoint
//   subdivision (four triangles, the medial one among them) and Sierpinski's (the three corner triangles, the medial one left out).
//   σ → K with probabilities per letter is a probabilistic context-free grammar on triangles; a derivation is a tree whose leaves are the
//   final triangles. parse_replacement recovers such a tree from the final triangles alone, bottom up, by finding the groups of
//   triangles each move makes and replacing them by their parent, until one triangle remains (the syntax analysis of a triangulation).
//
// Limits: gluing grows tree-like complexes (a triangle meets its parent along one edge and its children along others, and vertices
// are not identified across branches), so the language it generates is of simply connected patches, not of complexes with holes;
// parsing by replacement is greedy and exact geometry is assumed (a tolerance absorbs rounding only).

#pragma once

#include <physis/computational_geometry/planar.hpp>

#include <calx/random/stream.hpp>

#include <syntaxis/graphs/graphs.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace syntaxis::simplicial {

using Point = std::complex<double>;
using Triangle = std::array<std::size_t, 3>;
constexpr std::size_t none = std::numeric_limits<std::size_t>::max();

struct Complex {
  std::vector<Point> points;
  std::vector<Triangle> triangles;  // counter-clockwise
  std::size_t skipped = 0;          // rings or faces left out: not simple, or fewer than three distinct vertices
};

namespace detail {
inline double orient(Point a, Point b, Point c) { return (b.real() - a.real()) * (c.imag() - a.imag()) - (b.imag() - a.imag()) * (c.real() - a.real()); }
inline double angle_at(Point a, Point b, Point c) {  // the angle at a, in degrees
  const double d = std::abs(std::arg((b - a) / (c - a)));
  return d * 180 / std::numbers::pi;
}
}  // namespace detail

// ---- complexes from polygons and plane graphs -------------------------------------------------------------------------

// Rings triangulated by physis's constrained Delaunay triangulation of a polygon (ears, then Lawson's flips), so that a ring's
// triangles, and so its letters, are a function of its shape and not of where ear clipping began (unique up to cocircular ties,
// which keep the ears' diagonal; a rectangle's two triangles either way); points within `merge` of one another become one vertex.
inline Complex from_rings(const std::vector<std::vector<Point>>& rings, double merge = 1e-6) {
  Complex c;
  std::map<std::pair<std::int64_t, std::int64_t>, std::size_t> index;
  const auto vertex = [&](Point p) {
    const auto key = std::make_pair(std::llround(p.real() / merge), std::llround(p.imag() / merge));
    auto [it, fresh] = index.emplace(key, c.points.size());
    if (fresh) c.points.push_back(p);
    return it->second;
  };
  for (const auto& ring : rings) {
    // Each corner once: a closing point equal to the first, or a corner repeated, is one vertex.
    std::vector<std::size_t> id;
    for (const auto& p : ring) {
      const std::size_t v = vertex(p);
      if (id.empty() || id.back() != v) id.push_back(v);
    }
    while (id.size() > 1 && id.back() == id.front()) id.pop_back();
    if (id.size() < 3) { c.skipped++; continue; }
    std::vector<physis::computational_geometry::Point2> r;
    for (std::size_t v : id) r.push_back({c.points[v].real(), c.points[v].imag()});
    // Not simple if two edges that share no corner cross (a footprint drawn across itself): left out, as ears assume simplicity.
    bool simple = true;
    const std::size_t m = r.size();
    for (std::size_t i = 0; i < m && simple; i++)
      for (std::size_t j = i + 2; j < m && simple; j++)
        if (!(i == 0 && j == m - 1) && physis::computational_geometry::segments_cross(r[i], r[(i + 1) % m], r[j], r[(j + 1) % m])) simple = false;
    if (!simple) { c.skipped++; continue; }
    std::vector<std::array<std::size_t, 3>> T;
    try {
      T = physis::computational_geometry::delaunay_ring(r);
    } catch (const std::invalid_argument&) {  // a ring that touches or crosses itself: no ear
      c.skipped++;
      continue;
    }
    for (const auto& t : T) {
      Triangle tr{id[t[0]], id[t[1]], id[t[2]]};
      if (tr[0] == tr[1] || tr[1] == tr[2] || tr[0] == tr[2]) continue;
      if (detail::orient(c.points[tr[0]], c.points[tr[1]], c.points[tr[2]]) <= 0) continue;  // degenerate
      c.triangles.push_back(tr);
    }
  }
  return c;
}

// The bounded simple faces of a plane geometric graph (physis's plane faces), each triangulated as from_rings's: the graph's 2-complex.
inline Complex from_faces(const graphs::Graph& g, double max_area = std::numeric_limits<double>::infinity()) {
  std::vector<physis::computational_geometry::Point2> P;
  for (const auto& n : g.nodes) P.push_back({n.at.real(), n.at.imag()});
  std::vector<physis::computational_geometry::PlaneEdge> E;
  for (const auto& [u, v] : g.edges()) E.push_back({u, v, {}});
  Complex c;
  for (const auto& n : g.nodes) c.points.push_back(n.at);
  for (const auto& f : physis::computational_geometry::plane_faces(P, E)) {
    if (!(f.area > 0) || f.area > max_area) continue;  // the outer face, or a face larger than a block
    if (!f.simple || f.vertices.size() < 3) { c.skipped++; continue; }
    std::vector<physis::computational_geometry::Point2> ring;
    for (std::size_t v : f.vertices) ring.push_back(P[v]);
    for (const auto& t : physis::computational_geometry::delaunay_ring(ring)) {
      Triangle tr{f.vertices[t[0]], f.vertices[t[1]], f.vertices[t[2]]};
      if (detail::orient(c.points[tr[0]], c.points[tr[1]], c.points[tr[2]]) > 0) c.triangles.push_back(tr);
    }
  }
  return c;
}

// ---- letters ------------------------------------------------------------------------------------------------------------

struct Letter {
  int largest = 0, middle = 0, smallest = 0;  // the angles, in units of the resolution
  bool operator<(const Letter& o) const { return std::tie(largest, middle, smallest) < std::tie(o.largest, o.middle, o.smallest); }
  bool operator==(const Letter& o) const { return largest == o.largest && middle == o.middle && smallest == o.smallest; }
};

// The angles of a triangle at its vertices (degrees), in the triangle's order.
inline std::array<double, 3> angles(const Complex& c, const Triangle& t) {
  const Point a = c.points[t[0]], b = c.points[t[1]], d = c.points[t[2]];
  return {detail::angle_at(a, b, d), detail::angle_at(b, d, a), detail::angle_at(d, a, b)};
}

// A triangle's letter at a resolution in degrees: the largest and middle angles rounded, the smallest what makes 180.
inline Letter letter(const Complex& c, const Triangle& t, double degrees) {
  auto a = angles(c, t);
  std::sort(a.begin(), a.end(), std::greater<>());
  Letter l;
  const int total = static_cast<int>(std::lround(180.0 / degrees));
  l.largest = static_cast<int>(std::lround(a[0] / degrees));
  l.middle = static_cast<int>(std::lround(a[1] / degrees));
  l.smallest = total - l.largest - l.middle;
  return l;
}

inline std::string name(const Letter& l, double degrees) {
  const auto f = [&](int k) { const double d = k * degrees; return std::abs(d - std::round(d)) < 1e-9 ? std::to_string(static_cast<int>(std::lround(d))) : std::to_string(d); };
  return f(l.largest) + "·" + f(l.middle) + "·" + f(l.smallest);
}

// ---- the dual graph and gluing derivations -----------------------------------------------------------------------------

// For each triangle, the triangle across each of its edges (edge k joins vertices k and k + 1), or none.
inline std::vector<std::array<std::size_t, 3>> dual(const Complex& c) {
  std::map<std::pair<std::size_t, std::size_t>, std::vector<std::pair<std::size_t, int>>> by_edge;
  for (std::size_t t = 0; t < c.triangles.size(); t++)
    for (int k = 0; k < 3; k++) by_edge[std::minmax(c.triangles[t][k], c.triangles[t][(k + 1) % 3])].push_back({t, k});
  std::vector<std::array<std::size_t, 3>> out(c.triangles.size(), {none, none, none});
  for (const auto& [e, list] : by_edge)
    if (list.size() == 2) {
      out[list[0].first][list[0].second] = list[1].first;
      out[list[1].first][list[1].second] = list[0].first;
    }
  return out;
}

// The slot of edge k of a triangle: 0 the longest, 1 the middle, 2 the shortest (ties by position).
inline int slot(const Complex& c, const Triangle& t, int k) {
  double len[3];
  for (int j = 0; j < 3; j++) len[j] = std::abs(c.points[t[(j + 1) % 3]] - c.points[t[j]]);
  int rank = 0;
  for (int j = 0; j < 3; j++) if (len[j] > len[k] + 1e-12 || (std::fabs(len[j] - len[k]) <= 1e-12 && j < k)) rank++;
  return rank;
}

// The handedness of a triangle glued by edge k: 1 when the larger of the two angles at the edge's ends is at its first end.
inline int hand(const Complex& c, const Triangle& t, int k) {
  const auto a = angles(c, t);
  return a[k] >= a[(k + 1) % 3] ? 1 : 0;
}

struct Gluing {
  std::size_t parent = none, child = none;
  int parent_slot = 0, child_slot = 0, child_hand = 0;
};

struct Derivation {
  std::vector<std::size_t> roots;   // one per connected component, the first triangle reached
  std::vector<Gluing> gluings;      // in the order they were made
};

inline Derivation derive(const Complex& c) {
  const auto d = dual(c);
  Derivation D;
  std::vector<char> seen(c.triangles.size(), 0);
  for (std::size_t start = 0; start < c.triangles.size(); start++) {
    if (seen[start]) continue;
    seen[start] = 1;
    D.roots.push_back(start);
    std::vector<std::size_t> queue{start};
    for (std::size_t head = 0; head < queue.size(); head++) {
      const std::size_t t = queue[head];
      for (int k = 0; k < 3; k++) {
        const std::size_t u = d[t][k];
        if (u == none || seen[u]) continue;
        seen[u] = 1;
        int ku = 0;
        for (int j = 0; j < 3; j++) if (d[u][j] == t) ku = j;
        D.gluings.push_back({t, u, slot(c, c.triangles[t], k), slot(c, c.triangles[u], ku), hand(c, c.triangles[u], ku)});
        queue.push_back(u);
      }
    }
  }
  return D;
}

// ---- the gluing grammar -------------------------------------------------------------------------------------------------

struct Glued {  // what a slot holds: nothing, or a triangle of a letter by one of its slots and a handedness
  bool empty = true;
  Letter letter;
  int slot = 0, hand = 0;
  bool operator<(const Glued& o) const { return std::tie(empty, letter, slot, hand) < std::tie(o.empty, o.letter, o.slot, o.hand); }
};

struct GluingGrammar {
  double degrees = 15;
  std::map<Letter, double> letters;                                      // the letters' counts
  std::map<std::pair<Letter, int>, std::map<Glued, double>> productions;  // (letter, slot) → what is glued there, counted
  std::map<Letter, double> roots;
  // Bits per triangle: the letters' entropy; that of the glued triangles' letters; and that of a glued letter given its parent's
  // letter and slot (its own slot and handedness summed out). glued_entropy − conditional_entropy is the information the neighbour
  // carries about a letter, a mutual information, never negative.
  double entropy = 0, glued_entropy = 0, conditional_entropy = 0;
};

inline GluingGrammar learn_gluing(const Complex& c, double degrees = 15) {
  GluingGrammar G;
  G.degrees = degrees;
  const Derivation D = derive(c);
  std::vector<Letter> L(c.triangles.size());
  for (std::size_t t = 0; t < c.triangles.size(); t++) { L[t] = letter(c, c.triangles[t], degrees); G.letters[L[t]] += 1; }
  for (std::size_t r : D.roots) G.roots[L[r]] += 1;
  // Each triangle's slots: the child glued there, or nothing (its entry slot is not counted).
  std::vector<std::array<Glued, 3>> held(c.triangles.size());
  std::vector<int> entry(c.triangles.size(), -1);
  for (const auto& g : D.gluings) {
    held[g.parent][static_cast<std::size_t>(g.parent_slot)] = {false, L[g.child], g.child_slot, g.child_hand};
    entry[g.child] = g.child_slot;
  }
  for (std::size_t t = 0; t < c.triangles.size(); t++)
    for (int s = 0; s < 3; s++)
      if (s != entry[t]) G.productions[{L[t], s}][held[t][static_cast<std::size_t>(s)]] += 1;
  double n = 0;
  for (const auto& [l, k] : G.letters) n += k;
  for (const auto& [l, k] : G.letters) G.entropy -= k / n * std::log2(k / n);
  double glued = 0, h = 0;
  std::map<Letter, double> children;
  for (const auto& [ctx, dist] : G.productions) {
    std::map<Letter, double> by;  // the glued letter alone
    double m = 0;
    for (const auto& [x, k] : dist)
      if (!x.empty) by[x.letter] += k, children[x.letter] += k, m += k;
    for (const auto& [l, k] : by) h -= k * std::log2(k / m);
    glued += m;
  }
  G.conditional_entropy = glued > 0 ? h / glued : 0;
  for (const auto& [l, k] : children) G.glued_entropy -= k / glued * std::log2(k / glued);
  return G;
}

namespace detail {
// The angles a letter stands for when a triangle is placed: its classes' values, except that an angle taken to 0 (below half the
// resolution) stands for the middle of its class, a quarter of the resolution, taken from the largest; placed at 0, a sliver's two
// long sides would be parallel and never meet.
inline std::array<double, 3> placement_angles(const Letter& l, double degrees) {
  double A = l.largest * degrees, B = l.middle * degrees, C = l.smallest * degrees;
  const double least = degrees / 4;
  if (B < least) A -= least - B, B = least;
  if (C < least) A -= least - C, C = least;
  return {A, B, C};
}
// The angles a letter gives at the two ends of its slot-s edge, the larger first when hand = 1.
inline std::pair<double, double> end_angles(const Letter& l, double degrees, int s, int hand) {
  const auto [A, B, C] = placement_angles(l, degrees);
  // The edge of slot s is opposite the angle of rank s; the other two are at its ends.
  double x, y;
  if (s == 0) x = B, y = C;
  else if (s == 1) x = A, y = C;
  else x = A, y = B;
  return hand ? std::make_pair(x, y) : std::make_pair(y, x);
}
inline bool overlaps(const Complex& c, const std::array<Point, 3>& t) {
  const Point centre = (t[0] + t[1] + t[2]) / 3.0;
  const auto inside = [](const std::array<Point, 3>& tri, Point p) {
    return orient(tri[0], tri[1], p) > 1e-9 && orient(tri[1], tri[2], p) > 1e-9 && orient(tri[2], tri[0], p) > 1e-9;
  };
  for (const auto& tr : c.triangles) {
    const std::array<Point, 3> q{c.points[tr[0]], c.points[tr[1]], c.points[tr[2]]};
    if (inside(q, centre) || inside(t, (q[0] + q[1] + q[2]) / 3.0)) return true;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        if (shape::segments_cross(t[i], t[(i + 1) % 3], q[j], q[(j + 1) % 3])) return true;
  }
  return false;
}
inline double total(const std::map<Glued, double>& d) {
  double s = 0;
  for (const auto& [x, k] : d) s += k;
  return s;
}
}  // namespace detail

struct Grown {
  Complex complex;
  std::vector<Letter> letters;
  std::size_t refused = 0;  // triangles refused because they would overlap
  std::size_t draws = 1;    // the draws made, conditioned on size
  bool reached = true;      // whether one of them had at least the size asked
};

// The mean number of triangles glued to a triangle in the grammar's derivations: (T − roots)/T for a complex of T triangles, below
// one. The grammar's growth is a subcritical branching process (consistent, as Chi and Geman (1998) prove of every
// maximum-likelihood PCFG): its patches are finite with probability one, and small, unless conditioned on their size.
inline double mean_offspring(const GluingGrammar& G) {
  double children = 0, triangles = 0;
  for (const auto& [l, k] : G.letters) triangles += k;
  for (const auto& [ctx, dist] : G.productions)
    for (const auto& [x, k] : dist) if (!x.empty) children += k;
  return triangles > 0 ? children / triangles : 0;
}

namespace detail {
inline Grown grow_once(const GluingGrammar& G, calx::random::Stream& s, std::size_t max_triangles);
}

// A complex grown from a root triangle of the most common root letter, its base on [0, 1], to at most max_triangles. Conditioned on
// its size: up to `attempts` draws, the first of at least `at_least` triangles kept, else the largest (a draw from the grammar's law
// given a large patch, by rejection).
inline Grown grow(const GluingGrammar& G, calx::random::Stream& s, std::size_t max_triangles = 400, std::size_t at_least = 0, std::size_t attempts = 1) {
  Grown best;
  std::size_t a = 0;
  bool enough = false;
  for (; a < std::max<std::size_t>(1, attempts) && !enough; a++) {
    Grown g = detail::grow_once(G, s, max_triangles);
    enough = g.complex.triangles.size() >= at_least;
    if (g.complex.triangles.size() > best.complex.triangles.size() || a == 0) best = std::move(g);
  }
  best.draws = a;
  best.reached = enough;
  return best;
}

inline Grown detail::grow_once(const GluingGrammar& G, calx::random::Stream& s, std::size_t max_triangles) {
  Grown out;
  if (G.roots.empty()) return out;
  Letter root = G.roots.begin()->first;
  for (const auto& [l, k] : G.roots) if (k > G.roots.at(root)) root = l;
  // The root: its slot-0 edge on [0, 1], the triangle above it.
  const auto [a0, a1] = detail::end_angles(root, G.degrees, 0, 1);
  const auto apex = [&](Point u, Point v, double au, double av) {
    const double t = std::sin(av * std::numbers::pi / 180) / std::sin((au + av) * std::numbers::pi / 180);
    return u + (v - u) * std::polar(t, au * std::numbers::pi / 180);
  };
  Complex& c = out.complex;
  c.points = {Point{0, 0}, Point{1, 0}, apex(Point{0, 0}, Point{1, 0}, a0, a1)};
  c.triangles.push_back({0, 1, 2});
  out.letters.push_back(root);
  struct Open { std::size_t triangle; int edge; };
  std::vector<Open> queue;
  const auto open_slots = [&](std::size_t t, int entry_edge) {
    for (int k = 0; k < 3; k++) if (k != entry_edge) queue.push_back({t, k});
  };
  open_slots(0, -1);
  for (std::size_t head = 0; head < queue.size() && c.triangles.size() < max_triangles; head++) {
    const auto [t, k] = queue[head];
    const Triangle tri = c.triangles[t];
    const int sl = slot(c, tri, k);
    const auto it = G.productions.find({out.letters[t], sl});
    if (it == G.productions.end()) continue;
    // Draw what is glued here.
    double u = s.uniform() * detail::total(it->second);
    Glued choice;
    for (const auto& [x, n] : it->second) { if (u < n) { choice = x; break; } u -= n; choice = x; }
    if (choice.empty) continue;
    // The shared edge, from tri[k+1] to tri[k] for the new triangle (counter-clockwise on the other side).
    const Point p = c.points[tri[(k + 1) % 3]], q = c.points[tri[k]];
    const auto [ap, aq] = detail::end_angles(choice.letter, G.degrees, choice.slot, choice.hand);
    const Point w = apex(p, q, ap, aq);
    const std::array<Point, 3> nt{p, q, w};
    if (detail::orient(p, q, w) <= 0 || detail::overlaps(c, nt)) { out.refused++; continue; }
    c.points.push_back(w);
    c.triangles.push_back({tri[(k + 1) % 3], tri[k], c.points.size() - 1});
    out.letters.push_back(choice.letter);
    open_slots(c.triangles.size() - 1, 0);
  }
  return out;
}

// ---- replacement grammars -----------------------------------------------------------------------------------------------

enum class Move { pachner13 = 0, midpoint4 = 1, sierpinski3 = 2 };
inline const char* move_name(Move m) { return m == Move::pachner13 ? "Pachner 1–3" : m == Move::midpoint4 ? "midpoint" : "Sierpinski"; }

struct TreeNode {
  std::array<Point, 3> corners;
  int move = -1;                     // the move applied, or −1 for a leaf (a final triangle)
  std::vector<std::size_t> children;
  std::size_t depth = 0;
};

struct Replacement {
  std::vector<TreeNode> tree;        // tree[0] the root
  std::vector<std::array<Point, 3>> leaves;
  std::array<std::size_t, 3> moves{0, 0, 0};
};

inline std::vector<std::array<Point, 3>> subdivide(Move m, const std::array<Point, 3>& t) {
  const Point a = t[0], b = t[1], c = t[2], ab = (a + b) / 2.0, bc = (b + c) / 2.0, ca = (c + a) / 2.0;
  if (m == Move::pachner13) {
    const Point g = (a + b + c) / 3.0;
    return {{a, b, g}, {b, c, g}, {c, a, g}};
  }
  if (m == Move::midpoint4) return {{a, ab, ca}, {ab, b, bc}, {ca, bc, c}, {ab, bc, ca}};
  return {{a, ab, ca}, {ab, b, bc}, {ca, bc, c}};
}

// Derivation by a PCFG on triangles: each triangle stops with probability `stop[depth]`-free weight, else rewrites by a move drawn from
// `weights` (stop, Pachner 1–3, midpoint, Sierpinski); at max_depth every triangle stops.
inline Replacement replace(const std::array<Point, 3>& start, const std::array<double, 4>& weights, std::size_t max_depth, calx::random::Stream& s) {
  Replacement R;
  R.tree.push_back({start, -1, {}, 0});
  for (std::size_t i = 0; i < R.tree.size(); i++) {
    if (R.tree[i].depth >= max_depth) continue;
    const double total = weights[0] + weights[1] + weights[2] + weights[3];
    double u = s.uniform() * total;
    int pick = 0;
    for (int k = 0; k < 4; k++) { if (u < weights[static_cast<std::size_t>(k)]) { pick = k; break; } u -= weights[static_cast<std::size_t>(k)]; pick = k; }
    if (pick == 0) continue;
    const Move m = static_cast<Move>(pick - 1);
    R.tree[i].move = static_cast<int>(m);
    R.moves[static_cast<std::size_t>(m)]++;
    for (const auto& t : subdivide(m, R.tree[i].corners)) {
      R.tree[i].children.push_back(R.tree.size());
      R.tree.push_back({t, -1, {}, R.tree[i].depth + 1});
    }
  }
  for (const auto& n : R.tree) if (n.move < 0) R.leaves.push_back(n.corners);
  return R;
}

struct ReplacementParse {
  std::vector<TreeNode> nodes;                 // the leaves first, then each parent found
  std::vector<std::size_t> roots;              // what remains
  std::array<std::size_t, 3> moves{0, 0, 0};
  std::vector<std::vector<std::size_t>> steps; // per pass, the parents found
  bool derived = false;                        // one triangle remains
};

// The syntax analysis of a triangulation by the three moves, bottom up.
inline ReplacementParse parse_replacement(const std::vector<std::array<Point, 3>>& leaves, double tolerance = 1e-7) {
  ReplacementParse P;
  for (const auto& t : leaves) P.nodes.push_back({t, -1, {}, 0});
  std::vector<std::size_t> alive(leaves.size());
  for (std::size_t i = 0; i < alive.size(); i++) alive[i] = i;
  double scale = 0;
  for (const auto& t : leaves) for (const auto& p : t) scale = std::max(scale, std::abs(p));
  const double q = tolerance * std::max(1.0, scale);
  const auto key = [&](Point p) { return std::make_pair(std::llround(p.real() / q), std::llround(p.imag() / q)); };
  const auto tkey = [&](const std::array<Point, 3>& t) {
    std::array<std::pair<std::int64_t, std::int64_t>, 3> k{key(t[0]), key(t[1]), key(t[2])};
    std::sort(k.begin(), k.end());
    return k;
  };
  for (int pass = 0; pass < 64 && alive.size() > 1; pass++) {
    std::map<std::array<std::pair<std::int64_t, std::int64_t>, 3>, std::size_t> by_key;
    for (std::size_t i : alive) by_key[tkey(P.nodes[i].corners)] = i;
    const auto find = [&](const std::array<Point, 3>& t) { const auto it = by_key.find(tkey(t)); return it == by_key.end() ? none : it->second; };
    std::set<std::size_t> used;
    std::vector<std::size_t> made;
    for (std::size_t i : alive) {
      if (used.count(i)) continue;
      const auto& t = P.nodes[i].corners;
      bool done = false;
      // Each move, each vertex of t taken as the parent's corner a (or, for Pachner, as the centroid).
      for (int r = 0; r < 3 && !done; r++) {
        const Point a = t[static_cast<std::size_t>(r)], m1 = t[static_cast<std::size_t>((r + 1) % 3)], m2 = t[static_cast<std::size_t>((r + 2) % 3)];
        // Midpoint and Sierpinski: t = (a, ab, ca), so b = 2 ab − a, c = 2 ca − a.
        const std::array<Point, 3> parent{a, 2.0 * m1 - a, 2.0 * m2 - a};
        for (Move mv : {Move::midpoint4, Move::sierpinski3}) {
          const auto kids = subdivide(mv, parent);
          std::vector<std::size_t> ids;
          for (const auto& k : kids) { const std::size_t id = find(k); if (id == none || used.count(id)) break; ids.push_back(id); }
          if (ids.size() != kids.size()) continue;
          // Sierpinski's three are also midpoint's first three: take midpoint when its medial triangle is there.
          if (mv == Move::sierpinski3) {
            const std::size_t medial = find(subdivide(Move::midpoint4, parent)[3]);
            if (medial != none && !used.count(medial)) continue;
          }
          for (std::size_t id : ids) used.insert(id);
          P.nodes.push_back({parent, static_cast<int>(mv), ids, 0});
          made.push_back(P.nodes.size() - 1);
          P.moves[static_cast<std::size_t>(mv)]++;
          done = true;
          break;
        }
        if (done) break;
        // Pachner 1–3: t = (a, b, g) with g the centroid of (a, b, c), so c = 3g − a − b.
        const Point g = m2, b = m1;
        const Point c = 3.0 * g - a - b;
        const std::array<Point, 3> par{a, b, c};
        if (detail::orient(a, b, c) <= 0) continue;
        const auto kids = subdivide(Move::pachner13, par);
        std::vector<std::size_t> ids;
        for (const auto& k : kids) { const std::size_t id = find(k); if (id == none || used.count(id)) break; ids.push_back(id); }
        if (ids.size() == 3) {
          for (std::size_t id : ids) used.insert(id);
          P.nodes.push_back({par, static_cast<int>(Move::pachner13), ids, 0});
          made.push_back(P.nodes.size() - 1);
          P.moves[static_cast<std::size_t>(Move::pachner13)]++;
          done = true;
        }
      }
    }
    if (made.empty()) break;
    P.steps.push_back(made);
    std::vector<std::size_t> next;
    for (std::size_t i : alive) if (!used.count(i)) next.push_back(i);
    next.insert(next.end(), made.begin(), made.end());
    alive = std::move(next);
  }
  P.roots = alive;
  P.derived = alive.size() == 1;
  return P;
}

}  // namespace syntaxis::simplicial
