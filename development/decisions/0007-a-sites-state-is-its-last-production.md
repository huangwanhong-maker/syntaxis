# 0007. A Site's State Is Its Last Production

- **Status:** Accepted, provisionally (2026-11-26). Amends decision 0006, §2.1 to §2.3.
- **Date:** 2026-11-26
- **Decided by:** decision 0006's growth checked against a city's record, which showed it still wrong where sites overlap (the infrastructure's ADR 0171).

## 1. Context

Decision 0006 made each site an automaton whose condition is the state its own events leave, by the record's table of what followed. That table is learned from the states the participants share. Where a street is often held by other permits, the commonest state after a site's own end is still "under a regime", so the site ended twice: Melbourne's record grown so held 555 beginnings and 610 ends.

## 2. Decision

1. **A site's state is its last production**: the rule of its last event, or none before its first (a generic site's, on each participant). Rates per exposure are counted by last production and era.
2. **The rule** is drawn from the record's grammar, its contexts keyed by the states the participants carry, among the rules of the site's action that followed that production at a site of that action in the record.
3. **The participants' states** follow the record's table of what followed; they set the contexts and no site's state.

## 3. Alternatives Considered

- **A table of what followed each site's own events.** Learned from the shared states, it leaks overlaps; learned from a site's own successions, it is the last production.
- **A longer memory** (the last two productions). The record's sites hold too few events to learn it.

## 4. Consequences

- The planted check's two overlapping sites alternate and grow as many regimes as the record holds, lasting as long. `check_syntaxis_foams`, 39 checks; syntaxis 109.
- Not built: a site's memory beyond its last production.
