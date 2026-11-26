# foams

*Status: draft (2026-11-26; [decision 0004](../development/decisions/0004-grammars-whose-words-are-foams.md)).* Grammars whose words are foams: typed 2-complexes swept in time, read as the gluing of their vertices' stars ([docs/theory.md](../docs/theory.md), §7).

## 1. Concepts

- **Foam.** morphe's worldline complex of a record of encounters, with a kind for each vertex, participant and pair, a time for each vertex, and states of the participants from vertices on. The moments at the ends of a record (∂₀, ∂₁) may be marked.
- **Star.** A vertex's kind, the participants present with their roles (through, begins, ends, instant), and the corners of lenses that open, close or pass there. Its curvature κ = 1 − ends/2 + (opens + closes)/2 sums to the Euler characteristic (Proposition 7.1).
- **Keys.** Four canonical forms of a star (morphe's `canonical_form`): the whole star, the rule (the action), the match (what it consumes) and the context (the match with states and passes).
- **Production and sweep.** A star rewrites the network of the moment: its match is its left side, what it leaves its right side, the participants that go on its interface. A foam is the sweep of its stars in time; events on disjoint participants commute, so its order is a trace with Foata steps.
- **Grammar.** Rules counted with their matches and contexts, read at three levels by a hierarchical Dirichlet backoff whose escape and concentrations are chosen by leave-one-out. Level 1 is context-free (a rule given what it consumes), level 2 context-sensitive (given also the states and the lenses open around it). A grammar is judged on held-out events, or prequentially (each event priced by the grammar of every event before it), in bits per event; grammars are compared as a partial order by dominance.
- **Lives.** Each participant's worldline as a word of the rules it takes part in, with its place in each, compressed by Re-Pair (morphe's `repair`).

## 2. Content

| Header | Content |
|---|---|
| `foams.hpp` | `Foam`, `make` (labels on a record of encounters; the worldline complex built), `Foam::state`; `Corner`, `Star`, `stars` (roles, corners, κ, the four keys, exact for every star without corners between participants and otherwise where the canonical form's search is exact), `role`; `Replay`, `replay` (cells counted from the stars, χ both ways); `Steps`, `steps` (the Foata steps of the trace) |
| `grammar.hpp` | `Grammar`, `learn` (counts; γ in closed form, α₁ and α₂ by leave-one-out), `observe`, `p0`, `p1`, `p2`, `probability`; `HeldOut`, `held_out` (bits per event at each level, what the training saw, the counts each level reads); `Prequential`, `prequential` (each event in time priced by the grammar of those before it, refitted at given times); `Entropies`, `entropies`; `single_match`, `empty_match`, `Exposure`, `exposure`, `Rate`, `rates`; `Dominance`, `dominance` |
| `lives.hpp` | `describe` (a readable account of a star's rule or of the whole star); `passes_only`; `Lives`, `lives` (the words, their Re-Pair grammar, each rule's uses, the identity events of chosen kinds left out) |

## 3. Checks

`check_syntaxis_foams` (25):
- **A foam by hand.** Three buildings, a street and two junctions, with a birth, an alteration, two events on the street and censuses: 13 vertices, 21 segments and 10 lenses. Among its 11 events: 8 kinds of star, 7 rules, 3 matches and 5 contexts. κ = 0 at every event and 1 at each end, so Σκ = χ = 2. Its Foata steps have widths 3, 4, 2, 1, 1. Also held: its entropies by hand, the census at 5 per 61 building-years, its lives and their Re-Pair grammar with and without the censuses, and a readable account of its birth.
- **Worldlines pinched.** Three participants passing through one event: κ = −2 there and Σκ = χ = −3.
- **200 random records.** The cells counted from the stars equal the foam's. Every key is unchanged under renaming of the participants. The Foata steps equal morphe's causal order. Probabilities sum to one at every match and context. γ in closed form is no worse than any value of a fine grid.
- **The prequential code by hand.** Three events of kinds a, a, b on one participant, priced one by one as the grammar learns them with γ = α₁ = 1: 1 + log₂ 3 bits at level 0 and log₂(4/3) + log₂ 9 at level 1.
- **A planted context.** Events whose kind is set by their participant's state: level 2 predicts held-out events in under 0.5 bits, level 1 needs a bit. With the kind drawn apart from the state, level 2 costs at most 0.05 bits more. Prequentially, refitted every 200 events, level 2's code is a third of level 1's or less when the context is planted, and within 2% of it when it is not.
- **Dominance and refusals.** The dominance order of three candidates on two measures; refusals.

## 4. Dependencies

morphe's `topology` (worldline complexes), `combinatorics` (canonical forms) and `formal` (Re-Pair; traces in the checks); calx's random streams in the checks.

## 5. References

- Cartier, P. and Foata, D. (1969), *Problèmes combinatoires de commutation et réarrangements*, Lecture Notes in Mathematics 85, Springer.
- Dawid, A. P. (1984), Statistical theory: the prequential approach, *Journal of the Royal Statistical Society A* 147, 278–292.
- Diekert, V. and Rozenberg, G. (eds) (1995), *The Book of Traces*, World Scientific.
- Grünwald, P. D. (2007), *The Minimum Description Length Principle*, MIT Press.
- Larsson, N. J. and Moffat, A. (2000), Off-line dictionary-based compression, *Proceedings of the IEEE* 88(11), 1722–1732.
- MacKay, D. J. C. and Peto, L. C. B. (1995), A hierarchical Dirichlet language model, *Natural Language Engineering* 1, 289–308.
- Mazurkiewicz, A. (1977), Concurrent program schemes and their interpretations, DAIMI PB-78, Aarhus University.
- McKay, B. D. and Piperno, A. (2014), Practical graph isomorphism, II, *Journal of Symbolic Computation* 60, 94–112.
- Huang, W. (2026), An inquiry into the formation of evidence-based partial orders in the governance of human–AI evolutionary production, Knowledge Commons, doi:10.17613/73180-mhn03.
- Huang, W. (2026), A discussion of archival representation and generative inheritance in self-evolving AI inquiry, Knowledge Commons, doi:10.17613/9g6zv-df566.

## 6. Open Questions

- Nonterminals spanning several worldlines, found by compressing the foam itself (hyperedge replacement), beyond Re-Pair along single lives.
- A generator that grows foams from the grammar by the event operator, drawing matches from the network of the moment.
- Exposure for matches of several participants (the number of their occurrences in the moment over time).
