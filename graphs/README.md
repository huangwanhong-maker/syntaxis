# Geometric Graphs

Status: **draft** (2026-11-08). Module `syntaxis::graphs`.

## 1. Concepts

The words of a geometric graph grammar: simple graphs whose nodes carry a position and a symbol (0 a terminal, s > 0 the nonterminal of rule s) and, for a nonterminal, the frame of its rule. The theory, with proofs, is in [docs/theory.md](../docs/theory.md), §2.

## 2. Content

- `graphs.hpp`: `Node`, `Graph` (`add`, `connect`, `joined`, `edges`); `induced_edges`, `key_of` (the key of an induced pattern), `connected`, `crossings`; examples with known structure: `sierpinski(depth)`, `bridged_triangles(depth, ratio)` (copies joined by an edge rather than a shared corner), `grid(n, m)`, `transformed`.

## 3. Checks

The gasket's and the grid's counts of nodes and edges and their planarity; induced edges and connectivity of a cell; two diagonals crossing once; a cell's key unchanged by a similarity.

## 4. Dependencies

`syntaxis::shape`.

## 5. References

Fiser et al. (2016), §3.1.
