# 0006. Sites as Automata

- **Status:** Accepted, provisionally (2026-11-26). Amends decision 0005, §2.2; §2.1 to §2.3 amended by [decision 0007](0007-a-sites-state-is-its-last-production.md).
- **Date:** 2026-11-26
- **Decided by:** the growth of decision 0005 checked against a city's record, which showed it wrong where sites overlap (the infrastructure's ADR 0171).

## 1. Context

Decision 0005 grew a foam with each site's condition read from the states its participants carry, and drew rules from the grammar the record was read with. Melbourne's record of 555 permits, many of them on streets other permits also held, was grown with about a thousand permits begun and half again as many ended. A site whose streets lay under another permit's regime looked as if under its own, ended at its own rate what it never began, freed the streets and began again. The planted check had a single regime site and could not show it.

## 2. Decision

1. **Each site is an automaton.** Its condition is the set of states its own events have left on its participants, each from the state the participant carried into its first event, updated after each of its events by the table's commonest outcome. Rates per exposure are counted on that condition.
2. **The grammar the growth draws from** is learned by the backoff of decision 0004 with each event's context keyed by its own site's states.
3. **A production's precondition.** The rule is drawn among those whose action fits the site and that the record shows at a site of that action in that condition, weighed by the grammar.
4. **A participant's state in the grown foam** is a state some site holding it has left, none while no site has.

## 3. Alternatives Considered

- **A count of the regimes in force as the state.** It pairs ends with beginnings in number, not in kind: a site could still end another's regime.
- **Episodes as productions** (a regime and its end grown together). It would replay the record's patterns and leave the grammar unused.

## 4. Consequences

- The planted check gains a second regime site overlapping the first: each alternates, and each grows as many regimes as the record holds, lasting as long. `check_syntaxis_foams`, 39 checks; syntaxis 109.
- Not built: sites that see each other (a permit refused while another holds the street), which the record would have to show.
