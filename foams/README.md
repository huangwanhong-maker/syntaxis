# foams

*Status: draft (2026-11-26; [decision 0004](../development/decisions/0004-grammars-whose-words-are-foams.md), [decision 0005](../development/decisions/0005-patterns-across-worldlines-and-growing-foams.md), [decision 0006](../development/decisions/0006-sites-as-automata.md), [decision 0007](../development/decisions/0007-a-sites-state-is-its-last-production.md)).* Grammars whose words are foams: typed 2-complexes swept in time, read as the gluing of their vertices' stars ([docs/theory.md](../docs/theory.md), §7).

## 1. Concepts

- **Foam.** morphe's worldline complex of a record of encounters, with a kind for each vertex, participant and pair, a time for each vertex, and states of the participants from vertices on. The moments at the ends of a record (∂₀, ∂₁) may be marked.
- **Star.** A vertex's kind, the participants present with their roles (through, begins, ends, instant), and the corners of lenses that open, close or pass there. Its curvature κ = 1 − ends/2 + (opens + closes)/2 sums to the Euler characteristic (Proposition 7.1).
- **Keys.** Four canonical forms of a star (morphe's `canonical_form`): the whole star, the rule (the action), the match (what it consumes) and the context (the match with states and passes).
- **Production and sweep.** A star rewrites the network of the moment: its match is its left side, what it leaves its right side, the participants that go on its interface. A foam is the sweep of its stars in time; events on disjoint participants commute, so its order is a trace with Foata steps.
- **Grammar.** Rules counted with their matches and contexts, read at three levels by a hierarchical Dirichlet backoff whose escape and concentrations are chosen by leave-one-out. Level 1 is context-free (a rule given what it consumes), level 2 context-sensitive (given also the states and the lenses open around it). A grammar is judged on held-out events, or prequentially (each event priced by the grammar of every event before it), in bits per event; grammars are compared as a partial order by dominance.
- **Lives.** Each participant's worldline as a word of the rules it takes part in, with its place in each, compressed by Re-Pair (morphe's `repair`).
- **Patterns across worldlines.** The foam read as a process net (events as boxes, segments as wires) and compressed by digrams, Re-Pair on the net after gRePair: a pattern is two boxes joined by wires that recurs, nested, and its instances bind the worldlines that run through them (§7.6).
- **Growing.** A foam grown one event at a time by the event operator. Each site (any participant of a kind, or the participants of an event the record holds) is an automaton whose state is its last production. A site is chosen by its rate in that state and its era, and the rule is drawn from the context-sensitive grammar among the rules that followed that production at a site of that action in the record; the participants' states follow the record's table of what followed (§7.7).

## 2. Content

| Header | Content |
|---|---|
| `foams.hpp` | `Foam`, `make` (labels on a record of encounters; the worldline complex built), `Foam::state`; `Corner`, `Star`, `stars` (roles, corners, κ, the four keys, exact for every star without corners between participants and otherwise where the canonical form's search is exact), `role`; `Replay`, `replay` (cells counted from the stars, χ both ways); `Steps`, `steps` (the Foata steps of the trace) |
| `grammar.hpp` | `Grammar`, `learn` (counts; γ in closed form, α₁ and α₂ by leave-one-out), `observe`, `p0`, `p1`, `p2`, `probability`; `HeldOut`, `held_out` (bits per event at each level, what the training saw, the counts each level reads); `Prequential`, `prequential` (each event in time priced by the grammar of those before it, refitted at given times); `Entropies`, `entropies`; `single_match`, `empty_match`, `Exposure`, `exposure`, `Rate`, `rates`; `Dominance`, `dominance` |
| `lives.hpp` | `describe` (a readable account of a star's rule or of the whole star); `place_of`; `passes_only`; `Lives`, `lives` (the words, their Re-Pair grammar, each rule's uses, the identity events of chosen kinds left out) |
| `nets.hpp` | `Net`, `NetBox`, `NetPort`, `NetWire`, `NetRule`; `net` (the process net, events of chosen kinds that only pass left out), `compress`, `compressed` (digrams replaced by patterns until none recurs), `derive`, `inner_wires`; the sizes before and after |
| `grow.hpp` | `match_and_context` (the keys of an event not yet in a foam); `GrowSite`, `Growth`, `learn_growth` (sites as automata with their rates by last production and era, the rules each action admits and those that followed each production, the table of the states that followed and how far it is from a function (`unsure`), the grammar); `grow` (a foam grown on a calx stream) |

## 3. Checks

`check_syntaxis_foams` (41):
- **A foam by hand.** Three buildings, a street and two junctions, with a birth, an alteration, two events on the street and censuses: 13 vertices, 21 segments and 10 lenses. Among its 11 events: 8 kinds of star, 7 rules, 3 matches and 5 contexts. κ = 0 at every event and 1 at each end, so Σκ = χ = 2. Its Foata steps have widths 3, 4, 2, 1, 1. Also held: its entropies by hand, the census at 5 per 61 building-years, its lives and their Re-Pair grammar with and without the censuses, and a readable account of its birth.
- **Worldlines pinched.** Three participants passing through one event: κ = −2 there and Σκ = χ = −3.
- **200 random records.** The cells counted from the stars equal the foam's. Every key is unchanged under renaming of the participants. The Foata steps equal morphe's causal order. Probabilities sum to one at every match and context. γ in closed form is no worse than any value of a fine grid.
- **The prequential code by hand.** Three events of kinds a, a, b on one participant, priced one by one as the grammar learns them with γ = α₁ = 1: 1 + log₂ 3 bits at level 0 and log₂(4/3) + log₂ 9 at level 1.
- **A planted context.** Events whose kind is set by their participant's state: level 2 predicts held-out events in under 0.5 bits, level 1 needs a bit. With the kind drawn apart from the state, level 2 costs at most 0.05 bits more. Prequentially, refitted every 200 events, level 2's code is a third of level 1's or less when the context is planted, and within 2% of it when it is not.
- **A process net by hand.** Four episodes of a regime and its end on three streets, censuses of a building between: 8 events and 21 wires; the episode found four times, then two episodes in a row twice; two boxes and three wires left, the description 15 against 29.
- **200 random nets.** Every event in one box of the top level; the wires inside the patterns and those left add up to the net's; each pattern derives as many events as its instances hold; the same patterns when the participants are renamed.
- **Keys before the event.** The match and context of an event built from what it consumes equal those `stars` reads once it is there, on the hand foam and 200 random records with states.
- **A planted record grown.** Twenty buildings altered at 0.2 a year, a building born at year 10, three streets taken by regimes at 3 a year lasting 0.02 years, and two of them by regimes of their own at one a year lasting half a year, overlapping the others: grown 100 times, the alterations and each site's regimes as many as the record's on average, every regime and its end alternating at each site though they overlap, regimes lasting as long, the birth in about 1 − e⁻² of the histories and its building altered only once born; every grown foam's χ = Σκ; the same stream, the same foam. Its table of what followed, read as a function: with whether a regime is in force, some states out are not the table's commonest; with how many are in force, none.
- **A state an era never saw.** Regimes on three streets in the first ten of twenty years, the growth learned with a break at year 10: of 100 grown histories about a third carry a regime across the break, and each ends there at the window's rate; none is begun after the break, none left in force.
- **Dominance and refusals.** The dominance order of three candidates on two measures; refusals.

## 4. Dependencies

morphe's `topology` (worldline complexes), `combinatorics` (canonical forms) and `formal` (Re-Pair; traces in the checks); calx's random streams in the checks.

## 5. References

- Cartier, P. and Foata, D. (1969), *Problèmes combinatoires de commutation et réarrangements*, Lecture Notes in Mathematics 85, Springer.
- Dawid, A. P. (1984), Statistical theory: the prequential approach, *Journal of the Royal Statistical Society A* 147, 278–292.
- Diekert, V. and Rozenberg, G. (eds) (1995), *The Book of Traces*, World Scientific.
- Gillespie, D. T. (1977), Exact stochastic simulation of coupled chemical reactions, *Journal of Physical Chemistry* 81, 2340–2361.
- Grünwald, P. D. (2007), *The Minimum Description Length Principle*, MIT Press.
- Knuth, D. E. (1968), Semantics of context-free languages, *Mathematical Systems Theory* 2, 127–145.
- Maneth, S. and Peternek, F. (2018), Grammar-based graph compression, *Information Systems* 76, 19–45.
- Larsson, N. J. and Moffat, A. (2000), Off-line dictionary-based compression, *Proceedings of the IEEE* 88(11), 1722–1732.
- MacKay, D. J. C. and Peto, L. C. B. (1995), A hierarchical Dirichlet language model, *Natural Language Engineering* 1, 289–308.
- Mazurkiewicz, A. (1977), Concurrent program schemes and their interpretations, DAIMI PB-78, Aarhus University.
- McKay, B. D. and Piperno, A. (2014), Practical graph isomorphism, II, *Journal of Symbolic Computation* 60, 94–112.
- Huang, W. (2026), An inquiry into the formation of evidence-based partial orders in the governance of human–AI evolutionary production, Knowledge Commons, doi:10.17613/73180-mhn03.
- Huang, W. (2026), A discussion of archival representation and generative inheritance in self-evolving AI inquiry, Knowledge Commons, doi:10.17613/9g6zv-df566.

## 6. Open Questions

- Patterns of more than two boxes chosen at once, and the smallest net grammar against the greedy one.
- Sites a record never shows for several participants (a regime on a new set of streets), grown from the network of the moment's adjacencies.
- Exposure for matches of several participants (the number of their occurrences in the moment over time).
