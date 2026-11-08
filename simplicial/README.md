# Words of Triangles

Status: **draft** (2026-11-08). Module `syntaxis::simplicial`.

## 1. Concepts

Grammars whose words are 2-simplices. A triangle's letter is its shape class, its angles to a resolution. A 2-complex is read by a gluing grammar (which letters are glued to which edges of which letters) or rewritten by a replacement grammar (Pachner's 1–3 move, the midpoint subdivision, Sierpinski's corners), whose derivations the analyser recovers from the final triangles. The theory, with proofs, is in [docs/theory.md](../docs/theory.md), §6.

## 2. Content

- `simplicial.hpp`: `Complex`; `from_rings` and `from_faces` (by physis's constrained Delaunay triangulation; rings not simple left out and counted); `Letter`, `angles`, `letter`, `name`; `dual`, `derive`; `GluingGrammar` with its entropies (letters, glued letters, glued given the neighbour), `learn_gluing`, `mean_offspring`, `grow` (conditioned on size by rejection, saying how many draws it took); `Move`, `subdivide`, `replace`; `ReplacementParse`, `parse_replacement`.

## 3. Checks

A square's and an L's triangles and areas; a closed ring read once, a bow tie left out; a grid's derivation spanning its dual; entropies zero on one letter and the neighbour's information never negative; mean offspring (T − 1)/T; growth with one letter and no overlaps; a sliver placed with a finite smallest angle; leaf counts 1 + 2p + 3q + 2s; twenty derivations recovered from their leaves.

## 4. Dependencies

`syntaxis::graphs`, physis `computational_geometry`, calx `random`.

## 5. References

Lawson (1977); Lee and Lin (1986), *Discrete Comput. Geom.* 1; Pachner (1991), *Europ. J. Combin.* 12; Chi and Geman (1998); Miller (1955); Paninski (2003), *Neural Comput.* 15.
