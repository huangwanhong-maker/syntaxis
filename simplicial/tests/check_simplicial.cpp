// Checks of syntaxis/simplicial: 2-simplices as words.
//
// - Rings triangulated: a rectangle gives two triangles of the letter 90·45·45 (at 15°); an L-shaped ring, four triangles whose
//   areas sum to its own; a closed ring is read once, and rings of fewer than three corners or crossing themselves are left out.
// - A grid's faces, triangulated, form a complex whose dual is spanned by its derivation (T − components gluings).
// - A complex of one letter has a gluing grammar of entropy zero and conditional entropy zero, and a mean offspring of (T − 1)/T,
//   below one (a subcritical branching process); grown and conditioned on its size, every triangle has that letter and none
//   overlaps another.
// - Replacement: the leaves of a random derivation number 1 + 2·(Pachner) + 3·(midpoint) + 2·(Sierpinski), and the syntax analysis
//   of the leaves alone recovers a single root and the same count of each move, for twenty random derivations.

#include <syntaxis/simplicial/simplicial.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <complex>
#include <numbers>
#include <string>
#include <vector>

using namespace syntaxis;
using namespace syntaxis::simplicial;

int main() {
  calx::verification::Suite s;

  {
    const Complex r = from_rings({{{0, 0}, {4, 0}, {4, 4}, {0, 4}}});
    bool letters = r.triangles.size() == 2;
    for (const auto& t : r.triangles) letters = letters && name(letter(r, t, 15), 15) == "90·45·45";
    s.require("from_rings: a square is two triangles of the letter 90·45·45", letters);
    const Complex l = from_rings({{{0, 0}, {4, 0}, {4, 1}, {1, 1}, {1, 4}, {0, 4}}});
    double area = 0;
    for (const auto& t : l.triangles) area += detail::orient(l.points[t[0]], l.points[t[1]], l.points[t[2]]) / 2;
    s.require("from_rings: an L is four triangles", l.triangles.size() == 4);
    s.check("from_rings: whose areas sum to the L's (7)", area - 7, 1e-12);
    const Complex k = from_rings({{{0, 0}, {4, 0}, {4, 4}, {0, 4}, {0, 0}}, {{0, 0}, {1, 0}, {0, 0}}, {{10, 0}, {12, 2}, {12, 0}, {10, 2}}});
    std::vector<Point> star;
    for (int q = 0; q < 24; q++) star.push_back(std::polar(q % 2 ? 1.0 : 2.5 + 0.1 * (q % 5), 2 * std::numbers::pi * q / 24));
    const GluingGrammar S = learn_gluing(from_rings({star}), 15);
    s.require("gluing: on a star of several letters, the neighbour's information H(glued) − H(glued | neighbour) is not negative",
              S.letters.size() > 2 && S.glued_entropy - S.conditional_entropy >= -1e-12,
              std::to_string(S.letters.size()) + " letters, " + std::to_string(S.glued_entropy) + " − " + std::to_string(S.conditional_entropy) + " bits");
    {  // A sliver, 2 by 0.05: its letter 90·90·0 at 15°; grown, its triangles keep a finite size (no side parallel to another).
      const GluingGrammar Q = learn_gluing(from_rings({{{0, 0}, {2, 0}, {2, 0.05}, {0, 0.05}}, {{2, 0}, {4, 0}, {4, 0.05}, {2, 0.05}}}), 15);
      calx::random::Stream sq(3, 7, 0, 0);
      const Grown gq = grow(Q, sq, 50, 3, 50);
      double far = 0;
      for (const auto& p : gq.complex.points) far = std::max(far, std::abs(p));
      s.require("grow: a sliver's letter (90·90·0) is placed with its smallest angle a quarter of the resolution, its patch bounded",
                Q.letters.count(Letter{6, 6, 0}) == 1 && far < 100, name(Q.letters.begin()->first, 15) + ", farthest point " + std::to_string(far));
    }
    s.require("from_rings: a closed ring (its first point repeated) is still two triangles; a ring of two corners and a bow tie are left out",
              k.triangles.size() == 2 && k.skipped == 2, std::to_string(k.triangles.size()) + " triangles, " + std::to_string(k.skipped) + " skipped");
  }
  {
    const Complex c = from_faces(graphs::grid(5, 6));
    const Derivation D = derive(c);
    s.require("from_faces: a 5 × 6 grid's 20 cells are 40 triangles", c.triangles.size() == 40, std::to_string(c.triangles.size()));
    s.require("derive: the gluings span the dual (T − components)", D.gluings.size() + D.roots.size() == c.triangles.size() && D.roots.size() == 1);
    const GluingGrammar G = learn_gluing(c, 15);
    s.require("gluing: one letter (90·45·45)", G.letters.size() == 1 && name(G.letters.begin()->first, 15) == "90·45·45");
    s.check("gluing: the letters' entropy is zero", G.entropy, 1e-12);
    s.check("gluing: and the conditional entropy", G.conditional_entropy, 1e-12);
    s.check("gluing: the mean offspring is (T − 1)/T, subcritical", mean_offspring(G) - 39.0 / 40, 1e-12);
    calx::random::Stream st(9, 5, 0, 0);
    const Grown g = grow(G, st, 300, 30, 200);
    bool same = g.complex.triangles.size() > 1;
    for (const auto& t : g.complex.triangles) same = same && name(letter(g.complex, t, 15), 15) == "90·45·45";
    s.require("grow: every triangle grown has the letter", same, std::to_string(g.complex.triangles.size()) + " triangles, " + std::to_string(g.refused) + " refused");
    bool apart = true;
    const auto& P = g.complex.points;
    for (std::size_t i = 0; i < g.complex.triangles.size() && apart; i++)
      for (std::size_t j = i + 1; j < g.complex.triangles.size() && apart; j++) {
        const auto& a = g.complex.triangles[i];
        const auto& b = g.complex.triangles[j];
        for (int x = 0; x < 3 && apart; x++)
          for (int y = 0; y < 3 && apart; y++)
            if (shape::segments_cross(P[a[x]], P[a[(x + 1) % 3]], P[b[y]], P[b[(y + 1) % 3]])) apart = false;
      }
    s.require("grow: no two triangles overlap (no edges cross)", apart);
  }
  {
    calx::random::Stream st(11, 6, 0, 0);
    int ok_count = 0, ok_parse = 0;
    for (int r = 0; r < 20; r++) {
      const std::array<Point, 3> start{Point{0, 0}, Point{1, 0}, Point{0.3, 0.8}};
      const Replacement R = replace(start, {1.0, 1.0, 1.0, 1.0}, 4, st);
      ok_count += R.leaves.size() == 1 + 2 * R.moves[0] + 3 * R.moves[1] + 2 * R.moves[2];
      const ReplacementParse P = parse_replacement(R.leaves);
      ok_parse += P.derived && P.moves == R.moves;
    }
    s.require("replacement: the leaves number 1 + 2·Pachner + 3·midpoint + 2·Sierpinski (twenty derivations)", ok_count == 20, std::to_string(ok_count));
    s.require("parse_replacement: one root and the same moves recovered from the leaves alone (twenty derivations)", ok_parse == 20, std::to_string(ok_parse));
  }
  return s.done();
}
