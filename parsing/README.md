# The Syntax Analyser

Status: **draft** (2026-11-08). Module `syntaxis::parsing`.

## 1. Concepts

A geometric graph parsed by a given grammar: the rules applied in order of creation to every occurrence of their keys, greedily and disjointly, in passes to a fixpoint. The parse forest says which rule derives which part, how much is covered, and whether what remains is the axiom. It is invariant under similarity; it is not a complete membership test (membership is NP-complete already for regular DNLC grammars). The theory, with proofs, is in [docs/theory.md](../docs/theory.md), §5.

## 2. Content

- `parsing.hpp`: `Options`, `Parse` (instances, steps, roots, remaining, uses, covered, derived), `parse`.

## 3. Checks

The gasket and the bridged triangles derived by their own grammars with the same uses as learning; a turned, scaled and moved copy derived; a grid not derived by the gasket's grammar, covered in part; the forest's leaves the input's nodes; the steps replaying in order.

## 4. Dependencies

`syntaxis::induction`, `syntaxis::derivation`.

## 5. References

Fiser et al. (2016); Aalbersberg, Ehrenfeucht and Rozenberg (1986), *Discrete Appl. Math.* 13.
