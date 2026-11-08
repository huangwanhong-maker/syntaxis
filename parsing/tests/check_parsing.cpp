// Checks of syntaxis/parsing: the syntax analysis of geometric graphs.
//
// - A graph parsed by the grammar learned from it is derived, by the same applications of each rule as learning made; the bridged
//   triangles to their top triangle, every node covered.
// - So is any rotated, scaled and moved copy of it: parsing is invariant under similarity.
// - The Sierpinski grammar parses a grid only in part (not derived, covered less than wholly).
// - The parse forest's leaves are the input's nodes; its steps replay in order.

#include <syntaxis/induction/induction.hpp>
#include <syntaxis/parsing/parsing.hpp>

#include <calx/verification/check.hpp>

#include <string>
#include <vector>

using namespace syntaxis;

int main() {
  calx::verification::Suite s;
  const graphs::Graph S = graphs::sierpinski(4);
  const induction::Learned L = induction::learn(S);
  const parsing::Parse self = parsing::parse(L.grammar, S);
  bool same_uses = self.uses.size() == L.grammar.rules.size();
  for (std::size_t k = 0; same_uses && k < self.uses.size(); k++) same_uses = self.uses[k] == L.grammar.rules[k].uses;
  s.require("sierpinski(4) by its own grammar: derived, by the same applications of each rule as learning made", self.derived && same_uses,
            "remaining " + std::to_string(self.remaining.size()) + ", covered " + std::to_string(self.covered));
  const graphs::Graph B = graphs::bridged_triangles(4);
  induction::Options three;
  three.max_order = 3;
  const parsing::Parse bridged = parsing::parse(induction::learn(B, three).grammar, B);
  s.require("bridged triangles by their own grammar: derived to the top triangle, every node covered", bridged.derived && bridged.covered == 1.0 && bridged.remaining.size() == 3);
  const graphs::Graph T = graphs::transformed(S, {std::polar(0.37, -1.3), {12, 3}});
  const parsing::Parse moved = parsing::parse(L.grammar, T);
  s.require("a rotated, scaled and moved copy: derived (invariance under similarity)", moved.derived, "remaining " + std::to_string(moved.remaining.size()));
  const parsing::Parse grid = parsing::parse(L.grammar, graphs::grid(6, 6));
  s.require("a grid by the Sierpinski grammar: not derived, covered in part", !grid.derived && grid.covered < 1.0, "covered " + std::to_string(grid.covered));
  std::size_t leaves = 0;
  for (const auto& in : self.instances) leaves += in.children.empty();
  s.require("the parse forest's leaves are the input's nodes", leaves == S.size());
  bool ordered = true;
  for (std::size_t k = 1; k < self.steps.size(); k++) ordered = ordered && self.steps[k].before.size() < self.steps[k - 1].before.size();
  s.require("the steps replay in order, each on a smaller graph", ordered && !self.steps.empty(), std::to_string(self.steps.size()) + " steps");
  return s.done();
}
