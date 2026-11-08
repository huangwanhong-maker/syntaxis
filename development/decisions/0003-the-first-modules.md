# 0003. The First Modules

- **Status:** Accepted, provisionally
- **Date:** 2026-11-08
- **Decided by:** Claude, after decision 0002

## 1. Context

The first stage should carry G³ to each step the user named: the deterministic grammar learned by compression, its probabilistic and context-sensitive forms, the syntax analyser, and grammars whose words are 2-simplices. Each module needs an answer it can be checked against exactly.

## 2. Decision

Seven modules, in this order:

| Module | Content | Checked against |
|---|---|---|
| `shape` | Similarities of the plane, Procrustes alignment, canonical keys | keys invariant under similarity; alignment exact on similar sets |
| `graphs` | Geometric graphs, induced patterns, crossings; Sierpinski's gasket, bridged triangles, grids | their counts of nodes, edges and faces |
| `induction` | A straight-line grammar learned by compression | bridged triangles compressed to an axiom of three nodes by rules used 81, 27, 9, 3 times |
| `derivation` | Exact replay; expansion with the relaxed embedding | the input rebuilt, edges exactly, positions to rounding |
| `parsing` | The syntax analyser | a turned, scaled and moved copy covered and derived with the same uses |
| `stochastic` | Categories by shape, PCFG, PCSG (parent annotation, no crossing); sampling and its statistics | one recursive category on self-similar input; samples of many sizes; the PCSG crossing no more than the PCFG |
| `simplicial` | Complexes from polygons and plane faces, letters, the gluing grammar and its growth, the replacement grammar and its analyser | rectangles' and L's triangles and areas; entropies zero on one letter; subcritical offspring (T − 1)/T; leaf counts; twenty derivations recovered |

Each is a header-only target `syntaxis::<module>` with a check program `check_syntaxis_<module>`.

## 3. Alternatives Considered

- **Hyperedge replacement in the first stage.** Sierpinski's gasket showed its need (the theory note, Remark 3.4). It follows once node replacement and its analyser are in use, so that the two share their patterns and keys.
- **Grammars on 3-complexes.** They wait for a domain to need them (architecture's districts as 3-complexes are the likely first).

## 4. Consequences

- 70 checks in the first stage.
- Next: hyperedge replacement, an exact parser for small inputs, keys continuous in the shape space, molecular graphs for a chemistry domain, and 3-complexes.
