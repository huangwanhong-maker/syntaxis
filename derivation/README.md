# Derivation

Status: **draft** (2026-11-08). Module `syntaxis::derivation`.

## 1. Concepts

A grammar's productions applied: to rebuild what was learned, exactly, or to expand any symbol into a structure. The theory, with proofs, is in [docs/theory.md](../docs/theory.md), §3–4.

## 2. Content

- `derivation.hpp`: `replay` (the learned derivation run backwards, each instance's swallowed edges and frame restoring its occurrence) and `expand` (any symbol expanded with the relaxed embedding, an edge to an expanded node given to its pattern's nearest node).

## 3. Checks

The input's edges replayed exactly and its positions to 10⁻⁹ of its extent, on the gasket, bridged triangles and a grid; expansion making as many nodes as the input; the relaxed embedding keeping most edges, as stated.

## 4. Dependencies

`syntaxis::induction`.

## 5. References

Fiser et al. (2016), §3.3.
