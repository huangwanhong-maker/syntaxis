# CLAUDE.md

This file guides Claude Code in the syntaxis repository.

syntaxis is one of the foundation libraries of [Noetarium, the Science and Humanity Infrastructure](https://github.com/huangwanhong-maker/Noetarium-The-Science-and-Humanity-Infrastructure), which includes it as the git submodule `foundations/syntaxis`. It is a project in its own right and must make sense without the infrastructure.

## 1. Purpose and Scope

syntaxis implements grammars whose words are geometric or topological objects, and their learning, derivation and analysis:

- geometric graph grammars by node replacement, up to similarity;
- their induction from one structure as compression;
- their probabilistic and context-sensitive forms, and sampling from them;
- syntax analysis: a structure parsed by a grammar into a derivation forest;
- grammars whose words are 2-simplices: gluing and replacement, with their analysers;
- later, hyperedge replacement, exact parsing, and grammars on 3-complexes.

It builds on three sibling libraries:

- **calx**, for random streams;
- **morphe**, as the siblings' order requires (and later, for its formal languages and graph rewriting);
- **physis**, for segments, plane faces of plane graphs and the constrained Delaunay triangulation of polygons.

Out of scope:

| Topic | Where it belongs |
|---|---|
| Grammars on strings and their parsers (CFGs, Earley, BNF), double-pushout graph rewriting, regular complexes and their rewriting | morphe's `formal` and `topology` |
| Computational geometry for its own sake: triangulations, plane faces, spherical Delaunay | physis's `computational_geometry` |
| Search for the best grammar or design over a large space | heuresis |
| What a grammar is *of*: streets, buildings, molecules, organisms | the domains of the infrastructure |
| Drawing derivations and parse forests | the infrastructure's `visualization/` and its works |

The boundary is [decision 0002](development/decisions/0002-scope-boundaries-and-dependencies.md): a grammar whose words are geometric or topological, or a method to learn, derive or parse by one, is syntaxis's.

## 2. Invariants

1. **Dependencies.** syntaxis depends only on calx, morphe and physis among project libraries, and never on tyche, nomos, peira, eikon, heuresis or the infrastructure's domains, edges or applications.
2. **No particular objects.** Nothing in syntaxis names an object in nature; its examples are mathematical (gaskets, grids).
3. **Pure.** No I/O, no global state. Every draw comes from a calx stream given by the caller.
4. **Invariance stated.** Each method says under which group its results are invariant (similarities of the plane, here) and where rounding breaks that invariance.
5. **Exact where it can be.** A learned derivation replays its input exactly; a parse records what it covered and whether it derived the input. Nothing is reported as derived that was not.
6. **Truth known.** Checks are held to structures whose answers are known exactly (self-similar families, random derivations whose trees are kept). Reference values from other software come from running it (a script in `tests/fixtures/`), never from memory.
7. **Reuse.** What a sibling holds is used through its public headers, never written again here.
8. **Cited sources.** Every implemented method cites its source; equation or section numbers are given only where certain.

## 3. Design Philosophy

syntaxis shares its governing principle with its sibling repositories: **modular, extensible, revisable**.

- **Status vocabulary:**
  - `planned`: named, no content yet;
  - `seed`: an idea, possibly wrong;
  - `draft`: being worked out;
  - `stable`: tested and relied upon;
  - `deprecated`: superseded, with a pointer to its replacement.
- **One structure.** A grammar is a set of productions on words of some geometry. Learning, derivation and analysis are written once for that structure, and each kind of word (graph, 2-complex) supplies its patterns and its keys.

## 4. Repository Layout

```
syntaxis/
├── CMakeLists.txt   the project: one header-only target syntaxis::<module> per module (syntaxis_module); the siblings added once
├── shape/ · graphs/ · induction/ · derivation/ · parsing/ · stochastic/ · simplicial/
├── docs/            the theory note
└── development/     decisions · ideas · release notes
```

Each module directory holds:
- `README.md`: concepts first, then the content, then the checks;
- `include/syntaxis/<module>/`;
- `tests/check_<module>.cpp`, a check program on calx's `verification`, registered with CTest by `syntaxis_module` in the top-level `CMakeLists.txt`, with reference values in `tests/fixtures/` (`SYNTAXIS_FIXTURES`) where there are any.

The modules and their order are in [decision 0003](development/decisions/0003-the-first-modules.md).

## 5. Working Procedures

### 5.1 Decisions and Ideas

- Decisions are numbered records in `development/decisions/`. Records are never rewritten; a later record supersedes an earlier one.
- Ideas are files named `development/ideas/YYYY-MM-DD-slug.md`, with status `seed`, `promoted` or `abandoned`.
- **Release notes.** At each commit, add to `development/release-notes/YYYY-MM-DD.md` (one file per day, one short section per commit) a concise summary in plain prose of what is new, what one can now do and what works better. It is not a list of changes (the git log keeps those), and it is not written between commits. This convention is shared with the sibling repositories.

### 5.2 Change Discipline

- Before restructuring or changing an invariant (§2), propose the change and wait for the user's agreement.
- Keep each change to one concept.

### 5.3 Shared Conventions

calx, morphe, physis, tyche, peira, nomos, eikon, heuresis, syntaxis and the infrastructure keep these conventions identical:

- writing style;
- status vocabulary;
- decision, idea and release-note conventions.

Change them in every repository together.

### 5.4 Working Inside the Infrastructure

When you edit syntaxis inside the infrastructure's `foundations/syntaxis`:

1. commit in syntaxis first;
2. then commit the updated pointer in the infrastructure.

When a change needs a sibling to change as well, commit them in dependency order: calx and morphe first, then physis, then syntaxis.

## 6. Conventions

- **The plane.** Points are complex numbers, x + iy, in the caller's units; angles in letters are in degrees, counter-clockwise positive.
- **Similarities.** z ↦ az + b with a ≠ 0, orientation-preserving; a mirror image is another shape.
- **Triangles.** Vertex triples counter-clockwise; slots by edge length, 0 the longest.
- **Symbols.** 0 is the terminal (a vertex of the input); s > 0 is the nonterminal of rule s, numbered in order of creation.
- **Files.** A grammar is written as `noetarium.grammar/1` and a derivation or parse as `noetarium.derivation/1` by the works that save them.
- **Writing style.**
  - Headings are declarative, numbered and academic, never questions.
  - Concepts come before structure.
  - Label asides with *Rationale.*, *Remark.* or *Caveat.*
  - British spelling; English throughout.

## 7. Technology

- **Language: C++.** This follows ADR 0003 of the infrastructure.
- **Provisional:** C++20, header-only targets `syntaxis::<module>`; CMake 3.20 or later. Checks are programs on calx's `verification`, registered with CTest. On the development machine the compiler is GCC in WSL.
- **Dependency resolution.** Inside the infrastructure, syntaxis uses the targets of `foundations/calx`, `foundations/morphe` and `foundations/physis`, adding each that no larger project has added. Never nest submodules.
