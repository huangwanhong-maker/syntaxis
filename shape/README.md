# Shapes

Status: **draft** (2026-11-08). Module `syntaxis::shape`.

## 1. Concepts

Points in the plane compared up to the orientation-preserving similarities z ↦ az + b. Two points have no shape (the similarities are transitive on pairs); shape begins at three. The theory, with proofs, is in [docs/theory.md](../docs/theory.md), §2.

## 2. Content

- `shape.hpp`: `Point` (a complex number), `Similarity` (`then`, `inverse`, `scale`, `rotation`); `align`, the least-squares similarity between two labelled point lists (complex Procrustes) and its residual; `canonical`, the key of labelled points with edges up to similarity at a tolerance, with the canonical order and frame; `hull_area`; `segments_cross`, physis's.

## 3. Checks

A similarity applied to seven points recovered with zero residual; one key for 200 configurations and their images under random similarities, points permuted; a mirror image, a point moved beyond the tolerance, a label or an edge each change the key; hull areas; crossing and touching segments.

## 4. Dependencies

physis `computational_geometry`.

## 5. References

Kendall (1984), *Bull. London Math. Soc.* 16; Fiser et al. (2016), *SCCG '16*, §3.1.
