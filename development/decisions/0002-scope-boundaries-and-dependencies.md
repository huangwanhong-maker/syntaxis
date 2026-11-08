# 0002. Scope, Boundaries and Dependencies

- **Status:** Accepted
- **Date:** 2026-11-08
- **Decided by:**
  - **The request.** The user asked to begin with G³, their geometric graph grammar (a preliminary two-dimensional version in Python, presented at https://serendip.ngo/g3/), and to develop its deterministic context-free grammar into PCFGs and (P)CSGs, and into grammars whose words are 2-simplices: "try all, and visualize them". They asked for a syntax analyser that shows its process and results and saves its parse graphs and trees, for the code to be migrated to the project's style rather than reused, and, seeing that such grammars "can also parse other structures, like chemistry", for "a geometry and topology grammars domain". They agreed to a foundation named syntaxis.
  - **The criterion.** Claude, after the infrastructure's ADR 0083.

## 1. Context

A grammar of streets, a grammar of building plans and a grammar of molecules share their mathematics. Each needs:
- patterns that are the same up to a group of motions;
- productions that replace a symbol by a pattern placed by a motion;
- learning from examples;
- derivation;
- parsing.

None of this names an object in nature, so it does not belong in urban science, and morphe's string grammars and double-pushout rewriting do not compare patterns by their geometry. By the infrastructure's rule of reuse, a method that several domains need is written once in a foundation.

## 2. Decision

### 2.1 Scope

syntaxis holds grammars whose words are geometric or topological objects, and the methods to learn, derive and parse by them:
- geometric graph grammars by node replacement, their patterns compared up to the similarities of the plane;
- their induction from one structure as compression;
- their probabilistic and context-sensitive forms, and sampling;
- syntax analysis into derivation forests;
- grammars whose words are 2-simplices, by gluing and by replacement, with their analysers;
- later, hyperedge replacement, exact parsing, and grammars on 3-complexes.

**The criterion.** A method belongs here when it is a grammar over geometric or topological objects, or a method to learn, derive or parse by one. Grammars on strings and graph rewriting without geometry are morphe's. Computational geometry for its own sake is physis's. What a grammar is of, such as streets, buildings or molecules, is a domain's.

### 2.2 Dependencies

calx (random streams), morphe (as the siblings' order requires, and later its formal languages) and physis (segments, plane faces, the constrained Delaunay triangulation, which physis gained for this in its decision 0024). syntaxis stands beside tyche, peira, nomos, eikon and heuresis, and depends on none of them.

## 3. Alternatives Considered

- **A domain of geometric and topological grammars.** Rejected: it would name no object in nature, and a domain is one science's objects. The user's word "domain" is met by the `syntax` modules of urban science and architecture, which hold what is theirs (a street network as a graph, a footprint as a polygon) and use syntaxis for the rest.
- **Modules in morphe.** Rejected: morphe's formal languages are exact and combinatorial, while these grammars compare by geometry within a tolerance and are learned statistically; they would stretch morphe's scope and its invariants.
- **Modules in urban science.** Rejected: chemistry, architecture and the life sciences will need the same grammars.

## 4. Consequences

- The infrastructure gains a ninth foundation, `foundations/syntaxis`, a submodule with its own repository.
- urban science and architecture gain `syntax` modules on it, and the Urban Science Workbench a plug-in, *Geometric syntax*.
