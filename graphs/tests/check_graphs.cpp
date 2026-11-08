// Checks of syntaxis/graphs: geometric graphs, induced subgraphs, the examples.
//
// - The Sierpinski graph of depth d has (3^(d+1) + 3)/2 nodes and 3^(d+1) edges; the grid of n × m has n·m nodes and
//   n(m − 1) + m(n − 1) edges; neither crosses itself.
// - Induced subgraphs: their edges, their connectedness; a key is unchanged by a similarity of the whole graph.

#include <syntaxis/graphs/graphs.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace syntaxis::graphs;

int main() {
  calx::verification::Suite s;
  for (int d = 0; d <= 5; d++) {
    const Graph g = sierpinski(d);
    const double p = std::pow(3.0, d + 1);
    s.require("sierpinski(" + std::to_string(d) + "): (3^(d+1) + 3)/2 nodes and 3^(d+1) edges", g.size() == static_cast<std::size_t>((p + 3) / 2) && g.edge_count() == static_cast<std::size_t>(p),
              std::to_string(g.size()) + " nodes, " + std::to_string(g.edge_count()) + " edges");
  }
  s.require("sierpinski(3): no edges cross", crossings(sierpinski(3)) == 0);
  const Graph G = grid(4, 5);
  s.require("grid(4, 5): 20 nodes and 31 edges, none crossing", G.size() == 20 && G.edge_count() == 31 && crossings(G) == 0);
  const std::vector<std::size_t> square{0, 1, 5, 6}, apart{0, 2, 7};
  s.require("induced edges: a cell of the grid has four", induced_edges(G, square).size() == 4);
  s.require("connected: a cell is, three scattered nodes are not", connected(G, square) && !connected(G, apart));
  Graph X;
  X.add({0, 0}), X.add({2, 2}), X.add({0, 2}), X.add({2, 0});
  X.connect(0, 1), X.connect(2, 3);
  s.require("crossings: two diagonals cross once", crossings(X) == 1);
  const Graph T = transformed(G, {std::polar(3.0, 1.1), {5, 7}});
  s.require("key_of: a cell's key is unchanged by a similarity of the graph", key_of(G, square, 0.01).key == key_of(T, square, 0.01).key);
  return s.done();
}
