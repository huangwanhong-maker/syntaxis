// Shapes in the plane up to similarity: the geometry a geometric grammar compares (docs/theory.md, §2).
//
// A point of the plane is a complex number z. An orientation-preserving similarity is z ↦ a z + b with a ≠ 0 (a rotation by arg a, a
// scaling by |a|, a translation by b); they form a group, and two configurations are the same shape when one is carried onto the other
// by such a map. Reflections are not similarities here: a shape and its mirror image are two shapes, as in Fiser et al. (2016), §3.1.
//
// - align: the similarity carrying points P onto Q with least squares Σ |a p_i + b − q_i|², in closed form (the complex form of
//   Procrustes analysis; Umeyama 1991, IEEE PAMI 13(4), for the general case): with P̃, Q̃ centred, a = Σ conj(p̃) q̃ / Σ |p̃|², and the
//   relative residual √(Σ|a p_i + b − q_i|² / Σ |q̃|²).
// - canonical: the key of a labelled configuration with edges, the same for any two configurations carried onto each other by a
//   similarity, to within a tolerance. For every ordered pair (i, j) of distinct points, the frame z ↦ (z − p_i)/(p_j − p_i) sends p_i to
//   0 and p_j to 1; the points are quantised in it to a grid of side `tolerance` (relative to |p_j − p_i|), sorted, and written with their
//   labels and edges; the key is the least such writing over all pairs. It also gives the canonical order of the points and the frame,
//   so that two configurations with one key correspond point by point.
// - hull_area: the area of the convex hull (Andrew's monotone chain); segments_cross, physis's.
//
// Limits: quantisation puts two nearly equal shapes in different cells when a coordinate falls near a cell's edge; a tolerance is a
// resolution, not a guarantee. The cost of a key is O(n³ log n) for n points, meant for the small patterns of a grammar (n ≤ 10).

#pragma once

#include <physis/computational_geometry/planar.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace syntaxis::shape {

using Point = std::complex<double>;

struct Similarity {
  Point a{1, 0}, b{0, 0};  // z ↦ a z + b
  Point operator()(Point z) const { return a * z + b; }
  Similarity then(const Similarity& s) const { return {s.a * a, s.a * b + s.b}; }  // first this, then s
  Similarity inverse() const { return {1.0 / a, -b / a}; }
  double scale() const { return std::abs(a); }
  double rotation() const { return std::arg(a); }
};

struct Alignment {
  Similarity map;        // carries P onto Q
  double residual = 0;   // relative RMS misfit
};

inline Alignment align(std::span<const Point> P, std::span<const Point> Q) {
  if (P.size() != Q.size() || P.size() < 2) throw std::invalid_argument("syntaxis::shape::align: two configurations of the same size, at least two points");
  Point cp = 0, cq = 0;
  for (std::size_t i = 0; i < P.size(); i++) { cp += P[i]; cq += Q[i]; }
  cp /= static_cast<double>(P.size()), cq /= static_cast<double>(Q.size());
  Point num = 0;
  double den = 0, spread = 0;
  for (std::size_t i = 0; i < P.size(); i++) {
    num += std::conj(P[i] - cp) * (Q[i] - cq);
    den += std::norm(P[i] - cp);
    spread += std::norm(Q[i] - cq);
  }
  if (!(den > 0)) throw std::invalid_argument("syntaxis::shape::align: the first configuration is a single point");
  Alignment r;
  r.map.a = num / den;
  r.map.b = cq - r.map.a * cp;
  double ss = 0;
  for (std::size_t i = 0; i < P.size(); i++) ss += std::norm(r.map(P[i]) - Q[i]);
  r.residual = spread > 0 ? std::sqrt(ss / spread) : std::sqrt(ss);
  return r;
}

struct Canonical {
  std::string key;
  std::vector<std::size_t> order;  // order[k] = the point given k-th place
  Similarity frame;                // the canonical frame: carries the configuration to its canonical coordinates
  std::vector<Point> coordinates;  // the points in canonical order, in the canonical frame (unquantised)
};

// The key of points with labels and edges (index pairs) up to orientation-preserving similarity, at a relative tolerance.
inline Canonical canonical(std::span<const Point> p, std::span<const int> label, std::span<const std::pair<std::size_t, std::size_t>> edges, double tolerance) {
  const std::size_t n = p.size();
  if (n < 2 || label.size() != n) throw std::invalid_argument("syntaxis::shape::canonical: at least two points, each with a label");
  if (!(tolerance > 0)) throw std::invalid_argument("syntaxis::shape::canonical: a tolerance > 0");
  Canonical best;
  bool have = false;
  std::vector<std::size_t> idx(n), rank(n);
  std::vector<std::int64_t> qx(n), qy(n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) {
      if (i == j) continue;
      const Point d = p[j] - p[i];
      if (std::abs(d) < 1e-12) continue;
      const Similarity f{1.0 / d, -p[i] / d};
      for (std::size_t k = 0; k < n; k++) {
        const Point c = f(p[k]);
        qx[k] = std::llround(c.real() / tolerance);
        qy[k] = std::llround(c.imag() / tolerance);
        idx[k] = k;
      }
      std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
        if (qx[a] != qx[b]) return qx[a] < qx[b];
        if (qy[a] != qy[b]) return qy[a] < qy[b];
        return label[a] < label[b];
      });
      for (std::size_t k = 0; k < n; k++) rank[idx[k]] = k;
      std::vector<std::pair<std::size_t, std::size_t>> e;
      e.reserve(edges.size());
      for (const auto& [u, v] : edges) e.push_back(std::minmax(rank[u], rank[v]));
      std::sort(e.begin(), e.end());
      std::string key;
      key.reserve(n * 20 + e.size() * 8);
      for (std::size_t k = 0; k < n; k++) {
        key += std::to_string(qx[idx[k]]);
        key += ',';
        key += std::to_string(qy[idx[k]]);
        key += ':';
        key += std::to_string(label[idx[k]]);
        key += ';';
      }
      key += '|';
      for (const auto& [u, v] : e) {
        key += std::to_string(u);
        key += '-';
        key += std::to_string(v);
        key += ';';
      }
      if (!have || key < best.key) {
        have = true;
        best.key = std::move(key);
        best.order = idx;
        best.frame = f;
      }
    }
  if (!have) throw std::invalid_argument("syntaxis::shape::canonical: all points coincide");
  best.coordinates.resize(n);
  for (std::size_t k = 0; k < n; k++) best.coordinates[k] = best.frame(p[best.order[k]]);
  return best;
}

// The area of the convex hull.
inline double hull_area(std::vector<Point> p) {
  if (p.size() < 3) return 0;
  std::sort(p.begin(), p.end(), [](Point a, Point b) { return a.real() < b.real() || (a.real() == b.real() && a.imag() < b.imag()); });
  const auto cross = [](Point o, Point a, Point b) { return (a.real() - o.real()) * (b.imag() - o.imag()) - (a.imag() - o.imag()) * (b.real() - o.real()); };
  std::vector<Point> h(2 * p.size());
  std::size_t k = 0;
  for (std::size_t i = 0; i < p.size(); i++) {
    while (k >= 2 && cross(h[k - 2], h[k - 1], p[i]) <= 0) k--;
    h[k++] = p[i];
  }
  for (std::size_t i = p.size() - 1, t = k + 1; i > 0; i--) {
    while (k >= t && cross(h[k - 2], h[k - 1], p[i - 1]) <= 0) k--;
    h[k++] = p[i - 1];
  }
  double a = 0;
  for (std::size_t i = 0; i + 1 < k; i++) a += h[i].real() * h[i + 1].imag() - h[i + 1].real() * h[i].imag();
  return std::fabs(a) / 2;
}

// Whether the open segments ab and cd cross at a point interior to both (sharing an end is not crossing): physis's test.
inline bool segments_cross(Point a, Point b, Point c, Point d) {
  using physis::computational_geometry::Point2;
  return physis::computational_geometry::segments_cross(Point2{a.real(), a.imag()}, Point2{b.real(), b.imag()}, Point2{c.real(), c.imag()}, Point2{d.real(), d.imag()});
}

}  // namespace syntaxis::shape
