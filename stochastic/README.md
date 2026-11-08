# Probabilistic and Context-Sensitive Grammars

Status: **draft** (2026-11-08). Module `syntaxis::stochastic`.

## 1. Concepts

A learned grammar made generative: its rules merged into categories by shape (symbols ignored), which makes them recursive on self-similar input; a PCFG by the relative frequency of their uses; a PCSG conditioned on the parent's category, smoothed, refusing productions that cross what is drawn. The theory, with proofs, is in [docs/theory.md](../docs/theory.md), §4.

## 2. Content

- `stochastic.hpp`: `shape_key`, `Category`, `Model`, `model`, `conditional`; `Kind` (`pcfg`, `pcsg`), `SampleOptions`, `Sample`, `sample`; `Statistics`, `statistics`, `degree_distance`.

## 3. Checks

One recursive category from the bridged triangles' four rules, its samples of several sizes; probabilities summing to one as shares of uses; the gasket's sampled edges its base halved a whole number of times; the PCSG crossing no more than the PCFG over twenty seeds; an unseen context falling back to the PCFG.

## 4. Dependencies

`syntaxis::induction`, `syntaxis::derivation`, calx `random`.

## 5. References

Stolcke and Omohundro (1994); Chi and Geman (1998), *Comput. Linguist.* 24; Johnson (1998), *Comput. Linguist.* 24.
