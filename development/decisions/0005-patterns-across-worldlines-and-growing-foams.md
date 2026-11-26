# 0005. Patterns across Worldlines and Growing Foams

- **Status:** Accepted, provisionally (2026-11-26).
- **Date:** 2026-11-26
- **Decided by:** the infrastructure's user, who asked for the three things decision 0004 left unbuilt: patterns spanning several worldlines found from the complex itself, a generator that grows a foam from its grammar, and a second city's record read the same way (the infrastructure's ADR 0171).

## 1. Context

Decision 0004's lives compress each participant's worldline alone, so a pattern that binds several worldlines shows only as runs that several lives share. Its grammar could judge a record and could not grow one.

## 2. Decision

1. **`nets.hpp`.** The foam as a process net (events as boxes with ports named by signatures, segments as wires), compressed by digrams as Re-Pair compresses words and gRePair hypergraphs: the type with the most occurrences that share no box, two at least, becomes a pattern, until none recurs.
2. **`grow.hpp`.** A foam grown by the event operator: generic sites (any participant of a kind, for an event on one participant that opens or closes no lens) and concrete ones (an event's participants, pairs and newcomers); rates per exposure by condition and era; the rule drawn from the context-sensitive grammar among the rules of the site's shape; the states updated by the record's table of what followed; Gillespie's method on a calx stream.
3. **`foams.hpp`.** A star's shape, the rule's key with the event's kind left out. **`lives.hpp`.** `place_of`, a participant's place in a star, shared by the lives, the net and the growth.

## 3. Alternatives Considered

- **Hyperedge replacement with induction by subgraph matching.** It finds patterns of any size at once and needs matching on typed complexes. Digrams nest to any size, each step exact and cheap.
- **Patterns as the generator's productions.** Growing episodes (a regime and its end) as units would pair them without states, but would replay the record's patterns. The growth uses the grammar's rules and the states as attributes, and the patterns stay a check.
- **A generator drawing every match from the network of the moment.** A set of street spaces taken together is a fact of the record: drawn from all subsets of the moment's participants, nearly every match would be one no record shows. Sites of several participants are those the record shows.

## 4. Consequences

- `check_syntaxis_foams`, 39 checks; syntaxis 109.
- The infrastructure's urban science finds the patterns of a city's changes and grows its histories from its grammar (its ADR 0171).
- Not built: patterns of more than two boxes chosen at once; sites of several participants that no record shows.
