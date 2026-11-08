// Checks of syntaxis/shape: alignment and canonical keys up to similarity.
//
// - align recovers a similarity applied to random points, with zero residual; a configuration not similar has a residual.
// - canonical gives one key to a configuration and to its image under any similarity with its points permuted (200 random cases),
//   and the image's canonical coordinates agree with the original's; a mirror image of an asymmetric configuration has another key;
//   so does a configuration moved beyond the tolerance; labels and edges count in the key.
// - The convex hull's area: the unit square, a triangle, collinear points.

#include <syntaxis/shape/shape.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

using namespace syntaxis::shape;

int main() {
  calx::verification::Suite s;
  calx::random::Stream st(3, 1, 0, 0);

  {
    std::vector<Point> P, Q;
    const Similarity f{std::polar(2.5, 0.7), Point{3, -1}};
    for (int i = 0; i < 7; i++) { P.push_back({st.normal(), st.normal()}); Q.push_back(f(P.back())); }
    const Alignment a = align(P, Q);
    s.check("align: recovers a similarity applied to seven points", std::abs(a.map.a - f.a) + std::abs(a.map.b - f.b), 1e-12);
    s.check("align: with zero residual", a.residual, 1e-12);
    Q[3] += Point{0.5, 0};
    s.require("align: a configuration not similar has a residual", align(P, Q).residual > 1e-3);
  }
  {
    int same = 0, coords = 0;
    for (int r = 0; r < 200; r++) {
      std::vector<Point> p;
      std::vector<int> label;
      for (int i = 0; i < 6; i++) { p.push_back({st.uniform(), st.uniform()}); label.push_back(i % 2); }
      std::vector<std::pair<std::size_t, std::size_t>> e{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}, {0, 3}};
      const Canonical c = canonical(p, label, e, 1e-4);
      // The image under a random similarity, its points permuted.
      const Similarity f{std::polar(st.uniform(0.1, 10), st.uniform(-3.14, 3.14)), Point{st.normal() * 100, st.normal() * 100}};
      std::vector<std::size_t> perm(6);
      std::iota(perm.begin(), perm.end(), 0);
      for (std::size_t i = 5; i > 0; i--) std::swap(perm[i], perm[static_cast<std::size_t>(st.below(i + 1))]);
      std::vector<std::size_t> where(6);
      for (std::size_t i = 0; i < 6; i++) where[perm[i]] = i;
      std::vector<Point> q(6);
      std::vector<int> ql(6);
      for (std::size_t i = 0; i < 6; i++) { q[where[i]] = f(p[i]); ql[where[i]] = label[i]; }
      std::vector<std::pair<std::size_t, std::size_t>> qe;
      for (const auto& [u, v] : e) qe.push_back({where[u], where[v]});
      const Canonical d = canonical(q, ql, qe, 1e-4);
      same += c.key == d.key;
      double worst = 0;
      for (std::size_t k = 0; k < 6; k++) worst = std::max(worst, std::abs(c.coordinates[k] - d.coordinates[k]));
      coords += worst < 1e-9;
    }
    s.require("canonical: one key for a configuration and its image under any similarity, points permuted (200 cases)", same == 200, std::to_string(same));
    s.require("canonical: and the same canonical coordinates", coords == 200, std::to_string(coords));
  }
  {
    const std::vector<Point> p{{0, 0}, {3, 0}, {1, 2}, {0.2, 0.7}};
    const std::vector<int> label(4, 0);
    const std::vector<std::pair<std::size_t, std::size_t>> e{{0, 1}, {1, 2}, {2, 0}, {0, 3}};
    std::vector<Point> mirror;
    for (const auto& z : p) mirror.push_back(std::conj(z));
    const std::string k = canonical(p, label, e, 0.01).key;
    s.require("canonical: a mirror image is another shape", canonical(mirror, label, e, 0.01).key != k);
    std::vector<Point> moved = p;
    moved[2] += Point{0.3, 0};
    s.require("canonical: a point moved beyond the tolerance changes the key", canonical(moved, label, e, 0.01).key != k);
    const std::vector<int> other{0, 0, 1, 0};
    s.require("canonical: a label changes the key", canonical(p, other, e, 0.01).key != k);
    const std::vector<std::pair<std::size_t, std::size_t>> fewer{{0, 1}, {1, 2}, {2, 0}};
    s.require("canonical: an edge changes the key", canonical(p, label, fewer, 0.01).key != k);
  }
  {
    s.check("hull: the unit square's area", hull_area({{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}}) - 1, 1e-15);
    s.check("hull: a triangle's", hull_area({{0, 0}, {4, 0}, {0, 3}}) - 6, 1e-15);
    s.check("hull: collinear points have none", hull_area({{0, 0}, {1, 1}, {2, 2}}), 1e-15);
    s.require("segments_cross: crossing, and touching at an end", segments_cross({0, 0}, {2, 2}, {0, 2}, {2, 0}) && !segments_cross({0, 0}, {1, 1}, {1, 1}, {2, 0}));
  }
  return s.done();
}
