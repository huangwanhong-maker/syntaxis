# syntaxis

*σύνταξις, an arranging together; whence "syntax".* syntaxis holds grammars whose words are geometric: graphs in the plane whose patterns repeat up to rotation, scale and position, and 2-complexes whose letters are triangles. It learns such grammars from structures by compression, makes them probabilistic and context-sensitive, derives new structures from them, and parses structures by them, recording every derivation and parse. It is a foundation library of [Noetarium, the Science and Humanity Infrastructure](https://github.com/huangwanhong-maker/Noetarium-The-Science-and-Humanity-Infrastructure), and makes sense without it. Status: **draft** (2026-11-08).

## 1. Concepts

A string grammar rewrites letters standing in a line. Here the letters stand in the plane ([docs/theory.md](docs/theory.md)):

- **Geometric graph grammars.** A rule replaces one node by a pattern of nodes and edges, placed by a similarity $z \mapsto az + b$. A pattern is the same wherever it occurs up to such a map. Two points have no shape, since the similarities are transitive on pairs, so a pattern has three nodes at least (§2).
- **Learning as compression.** What repeats in one graph becomes a rule, chosen by the description it saves, until nothing repeats. The grammar found is straight-line: its language is the input alone, which its derivation rebuilds exactly (§3, after Fiser et al. 2016).
- **The hierarchy.** Rules of the same shape are merged into categories, which makes them recursive on self-similar input. Weighted by their use, they form a PCFG; conditioned on their parent's category, and refusing to cross what is drawn, a PCSG (§4).
- **Syntax analysis.** A structure is parsed by a grammar into a forest that says which rules derive which parts. The parse is invariant under similarity: a turned and scaled copy parses as the original (§5).
- **2-simplices as words.** A triangle's shape class is a letter. A 2-complex is a text, read by a gluing grammar (which letters are glued to which edges, with its entropies) or by a replacement grammar (a triangle rewritten by Pachner's, the midpoint's or Sierpinski's move), whose derivations the analyser recovers from the final triangles alone (§6).

A method belongs here when it is a grammar over geometric or topological objects, or its learning, derivation or analysis. What it is a grammar *of*, such as streets, buildings or molecules, belongs to the domains ([decision 0002](development/decisions/0002-scope-boundaries-and-dependencies.md)).

## 2. Modules

| Module | Holds | Status |
|---|---|---|
| [shape](shape/README.md) | Similarities of the plane, Procrustes alignment, canonical keys of labelled point sets up to similarity | draft |
| [graphs](graphs/README.md) | Geometric graphs; induced patterns and their keys; crossings; examples (Sierpinski's gasket, bridged triangles, grids) | draft |
| [induction](induction/README.md) | Learning a straight-line geometric graph grammar by compression: isogroups, selection by description length, contraction, the derivation | draft |
| [derivation](derivation/README.md) | Exact replay of a learned derivation; expansion of any symbol with the relaxed embedding | draft |
| [parsing](parsing/README.md) | The syntax analyser: a graph parsed by a grammar into a forest, with coverage and whether it is derived | draft |
| [stochastic](stochastic/README.md) | Categories by shape, the PCFG by relative frequency, the PCSG by parent annotation and a non-crossing condition; sampling and its statistics | draft |
| [simplicial](simplicial/README.md) | 2-complexes from polygons and plane graphs; letters; the gluing grammar with its entropies and conditioned growth; the replacement grammar and its analyser | draft |
| hyperedge | Hyperedge-replacement grammars, which share nodes between copies (Sierpinski's gasket), and their parsing | planned |
| exact | An exact parser for small inputs by backtracking over the greedy analyser's choices | planned |
| solid | Grammars on 3-complexes whose words are tetrahedra | planned |

## 3. Checks

`node development/tools/run-cpp.mjs foundations/syntaxis foundations/syntaxis` from the infrastructure's root builds and runs the 70 checks. Each holds against an answer known exactly: canonical keys invariant under similarity; self-similar structures compressed to an axiom of three nodes; derivations replayed to rounding; a turned, scaled and moved copy parsed wholly by the grammar learned from the original; the leaf count of a replacement derivation; and twenty derivations of triangles recovered by the analyser from their leaves alone.

## 4. Dependencies

calx (random streams), morphe (as the siblings' order requires) and physis (segments and their crossings, plane faces of plane graphs, the constrained Delaunay triangulation of polygons). Header-only C++20.

## 5. Licence

MIT; see [LICENSE](LICENSE) and [CONTRIBUTORS.txt](CONTRIBUTORS.txt).
