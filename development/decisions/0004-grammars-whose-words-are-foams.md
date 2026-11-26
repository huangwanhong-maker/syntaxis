# 0004. Grammars Whose Words Are Foams

- **Status:** Accepted, provisionally (2026-11-26).
- **Date:** 2026-11-26
- **Decided by:** the infrastructure's user, who asked whether a 2-complex generative grammar already exists and, if so, to learn one and use it to analyse how the structure of a city's events in spacetime is generated (the infrastructure's ADR 0168).

## 1. Context

syntaxis's `simplicial` holds grammars whose words are triangles in the plane, with letters by shape: a gluing grammar that grows tree-like patches, and a replacement grammar with its analyser. A city's record of events makes a 2-complex of another kind. The infrastructure's event complex, morphe's worldline complex of a record of encounters, is abstract and typed, has holes, and is generated in time, one event at a time. The infrastructure's design of event spacetime reads each event as a rule schema, whose event operator attaches the event's vertex with its link (its Definition 4.8). It does not learn those rules from a record.

## 2. Decision

A new module `foams`, header-only, on morphe's `topology`, `combinatorics` and `formal`:

- a foam, the labelled worldline complex, and the star of each vertex with its roles and corners;
- four canonical keys of each star (star, rule, match, context) and its curvature κ, with χ = Σκ exact;
- the Foata steps of the record's trace;
- the stochastic grammar at three levels by a hierarchical Dirichlet backoff, its escape in closed form and its concentrations by leave-one-out;
- judgement in bits per event, on held-out events and prequentially (each event priced by the grammar of the events before it, refitted at given times), plug-in entropies, rates per exposure, and the partial order of grammars by dominance;
- the lives of the participants as words compressed by Re-Pair (morphe's decision 0022).

## 3. Alternatives Considered

- **Extending `simplicial`.** Its words are geometric triangles with letters by angle. A foam's letters are combinatorial types of stars, and its cells are lenses whose boundaries are two worldlines. One structure for both would hide the difference.
- **Hyperedge replacement over the whole foam.** It is the general form, and the author's draft of grammars of records takes it. Its induction needs subgraph matching on typed complexes and an occurrence condition. The star grammar is exact, local and checkable now, and Re-Pair over lives gives its hierarchy along worldlines. Nonterminals across worldlines remain planned.
- **A split by date as the only test.** A record whose sources begin at different times confounds it: the rules of a source that begins after the date are new to a grammar learned before it, and the test measures the source's arrival. The prequential judgement charges each new rule once and then learns it.
- **Probabilities by relative frequency alone.** A record's later events bring rules and contexts its earlier ones never had. The backoff gives every event a probability, with its parameters chosen by leave-one-out, not by hand.

## 4. Consequences

- `check_syntaxis_foams`, 25 checks (the module's README); syntaxis 95.
- The infrastructure's urban science reads a city's event complex as a foam and learns its grammar (its ADR 0168).
- A state given at a vertex is the one that vertex's context sees; a caller who wants a context to hold only what came before gives the state carried into each event (`Foam::states`).
- Not built: nonterminals across worldlines, a generator by the event operator, exposure for matches of several participants.
