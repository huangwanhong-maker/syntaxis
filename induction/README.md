# Induction

Status: **draft** (2026-11-08). Module `syntaxis::induction`.

## 1. Concepts

A grammar learned from one geometric graph by compression: what repeats up to similarity becomes a rule, chosen by the description it saves, until nothing repeats. The grammar is straight-line, the graph form of Sequitur and Re-Pair, and its derivation is recorded for exact replay. The theory, with proofs, is in [docs/theory.md](../docs/theory.md), §3.

## 2. Content

- `induction.hpp`: `Rule`, `Grammar`, `Instance`, `Step`, `Options` (tolerance, orders, pruning); `occurrences` (every connected node set of a key, exhaustively, shared with parsing); `contract`; `learn`, returning the grammar, the derivation forest with the edges each instance swallowed, and the steps.

## 3. Checks

Bridged triangles of depth 4 (243 nodes) learned as four rules used 81, 27, 9 and 3 times with an axiom of three nodes; the gasket compressed below a fifth of its nodes and a third of its description, the same rules learned from a turned and scaled copy; a grid compressed, random points scarcely; the derivation's leaves the input's nodes.

## 4. Dependencies

`syntaxis::graphs`.

## 5. References

Fiser, Benes, Garcia Galicia, Abdul-Massih, Aliaga and Krs (2016), "Learning geometric graph grammars", *SCCG '16*, §4–6; Nevill-Manning and Witten (1997); Larsson and Moffat (2000); Charikar et al. (2005).
