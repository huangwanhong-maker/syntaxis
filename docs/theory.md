# The Theory of syntaxis: Grammars Whose Words Are Geometric

*Status: draft (2026-11-08). The results the code relies on, with their hypotheses and proofs or citations, and the limits found while checking them. Each section names the module that implements it.*

## 1. Grammars over Geometric Objects

A formal grammar rewrites words. A string is the one-dimensional case: its letters stand in a line, and a production replaces a letter by a word glued in at the letter's two ends. syntaxis holds grammars whose words are geometric objects in the plane:

- **geometric graphs**, whose letters are labelled points joined by edges, rewritten by node replacement (§2–5);
- **2-complexes**, whose letters are triangles up to similarity, glued along edges or rewritten into triangulations of themselves (§6).

Geometry enters twice. A pattern is the same wherever it occurs *up to a similarity*: a rule found at one place applies at any position, angle and scale. And geometry constrains what may be derived: two edges may not cross, and a triangle may not overlap another. The second is a context condition that the productions alone cannot state.

The grammars come in the order of the classical hierarchy:

| Grammar | Its productions | How they are found | Module |
|---|---|---|---|
| straight-line, context-free | one per nonterminal; its language is the input alone | compression of one graph | `induction` |
| probabilistic context-free (PCFG) | rules merged into categories by shape, each with a probability | relative frequency over the derivation | `stochastic` |
| probabilistic, context-sensitive (PCSG) | a production's probability depends on its parent's category, and it is refused if it crosses what is drawn | parent annotation, smoothed | `stochastic` |
| gluing grammar on 2-simplices | a triangle's letter entered by one edge glues letters to its other edges | relative frequency over a spanning tree of the dual graph | `simplicial` |
| replacement grammar on 2-simplices | a triangle rewrites into three or four by a move | stated, with probabilities | `simplicial` |

Each grammar comes with a syntax analyser: given a structure, it recovers the derivation that a grammar gives it (§5, §6.4).

## 2. Shapes up to Similarity, and Geometric Graphs

*Modules `shape` and `graphs`.*

**Definition 2.1 (similarities).** Identify the plane with $\mathbb{C}$. The orientation-preserving similarities are the maps $g(z) = az + b$ with $a \in \mathbb{C}^\times$, $b \in \mathbb{C}$; they form the group $\mathrm{Sim}^+(2) \cong \mathbb{C} \rtimes \mathbb{C}^\times$ (`shape::Similarity`, with `then` and `inverse`). Reflections are excluded, so a shape and its mirror image are two shapes, as in Fiser et al. (2016), §3.1.

**Proposition 2.2 (pairs have no shape).** $\mathrm{Sim}^+(2)$ acts simply transitively on ordered pairs of distinct points.

*Proof.* For $p \ne q$ the unique similarity with $p \mapsto 0$, $q \mapsto 1$ is $z \mapsto (z - p)/(q - p)$. ∎

So every edge is similar to every other edge, and shape begins at three points. The shapes of $k \ge 3$ labelled points form Kendall's shape space $\Sigma^k_2 \cong \mathbb{CP}^{k-2}$ (Kendall 1984, *Bull. London Math. Soc.* 16). Learning therefore chooses patterns of three nodes or more (`induction::Options::min_order`).

**Definition 2.3 (canonical key).** Let $X = (x_1, \dots, x_k)$ be points with labels $\lambda_i$ and edges $E$ among them, and let $\tau > 0$ be a tolerance. For each ordered pair $(i, j)$, $i \ne j$, map $x_i \mapsto 0$ and $x_j \mapsto 1$. Then round each image to the lattice $\tau\mathbb{Z}[i]$, and sort the points by their rounded coordinates and labels. This gives a word $w_{ij}$ in rounded coordinates, labels and edges. The *key* of $X$ is $\min_{(i,j)} w_{ij}$ (`shape::canonical`). The minimising pair fixes the *canonical order* of the nodes and the *frame*, the similarity from canonical coordinates to $X$.

**Proposition 2.4 (invariance).** For every $g \in \mathrm{Sim}^+(2)$, $\mathrm{key}(gX) = \mathrm{key}(X)$, except where an image coordinate lies within rounding error of a cell boundary of $\tau\mathbb{Z}[i]$.

*Proof.* The map normalising the pair $(i, j)$ of $gX$ is the one for $X$ composed with $g^{-1}$, so the normalised images are equal and so are the words. ∎

*Caveat.* Rounding is not continuous. Two occurrences that differ by less than $\tau$ can still fall on two sides of a cell boundary and receive different keys. Learning is therefore sensitive to $\tau$. On the bridged triangles (`graphs::bridged_triangles`), the hierarchy is found exactly at $\tau = 0.06$ but broken at $0.08$. Keys in a quotient that is continuous, such as Procrustes distances clustered with a margin, are planned.

**Definition 2.5 (geometric graph).** A *geometric graph* over terminals $\Sigma$ and nonterminals $N$ is $G = (V, E, x, \lambda, \phi)$:
- $(V, E)$ is a simple undirected graph;
- $x : V \to \mathbb{C}$ gives the positions;
- $\lambda : V \to \Sigma \sqcup N$ gives the symbols;
- $\phi$ gives each nonterminal node a frame in $\mathrm{Sim}^+(2)$ (`graphs::Graph`, `graphs::Node`).

A *pattern* is the subgraph induced on a connected node set, compared with other patterns by its key (`graphs::key_of`).

## 3. Induction as Compression

*Module `induction`.* After Fiser, Benes, Garcia Galicia, Abdul-Massih, Aliaga and Krs (2016), "Learning geometric graph grammars", *SCCG '16*, §4–6.

**Definition 3.1 (rule and contraction).** A *rule* $A \to P_A$ gives a nonterminal $A$ a pattern $P_A$ in its canonical frame: nodes with symbols, frames for its nonterminal nodes, and edges. To *contract* a node-disjoint set of occurrences $O_1, \dots, O_m$ of $P_A$ in $G$ (`induction::contract`):
- replace each $O_t$ by one node labelled $A$, at the centroid of $O_t$, whose frame maps $P_A$'s canonical frame onto $O_t$;
- give every edge from $O_t$ to the rest of the graph to the new node;
- record the edges inside $O_t$ as *swallowed* by the new instance.

This is node replacement in the sense of NLC grammars (Janssens and Rozenberg 1980, *Inform. Control* 47). The connection relation is not stated by the rule: it is recorded for each instance, which is what makes the derivation exact (Proposition 3.3).

**Definition 3.2 (description length).** A graph costs $|V| + |E|$, and a grammar costs its axiom plus $\sum_A (|V(P_A)| + |E(P_A)|)$. Replacing $m$ occurrences of a pattern of $k$ nodes and $e$ edges saves $m(k - 1 + e)$ and costs a rule of $k + e$. The *gain* is $m(k - 1 + e) - (k + e)$.

Learning repeats the following until no gain is positive (`induction::learn`):
- **Detection.** Grow connected node sets from every edge, one neighbour at a time, keeping the most frequent keys at each size (Fiser et al.'s pruning by frequency).
- **Selection.** Choose the key of greatest gain, its occurrences re-enumerated exhaustively (`induction::occurrences`) and taken greedily node-disjoint.
- **Contraction.** Contract them as in Definition 3.1.

Counting edges in the gain makes a triangle worth more than a path through the same three nodes. That is the role of Fiser et al.'s convex-hull heuristic, which prefers bifurcations to chains.

**Proposition 3.3 (straight-line grammar, exact replay).**
- (i) Each nonterminal of a learned grammar has exactly one rule, and the nonterminals are ordered by creation, each rule using only earlier ones. The grammar is therefore acyclic and deterministic, and its language is one graph.
- (ii) Expanding the axiom, with the recorded frames placing each rule and the recorded swallowed edges as the connection relation, rebuilds the input: the edges exactly, and the positions to rounding (`derivation::replay`).

*Proof.* (i) holds by construction. For (ii), each contraction is a bijection between the occurrence's nodes and the rule's nodes in canonical order. The frame carries the canonical coordinates back onto the occurrence's positions. The edges that left the occurrence are those now incident to the new node, and they return to the endpoints recorded with the instance. Induct on the steps in reverse. ∎

So induction is grammar-based compression, the graph form of Sequitur and Re-Pair on strings (Nevill-Manning and Witten 1997, *J. Artif. Intell. Res.* 7; Larsson and Moffat 2000, *Proc. IEEE* 88). Finding the smallest grammar is NP-hard already for strings (Charikar et al. 2005, *IEEE Trans. Inf. Theory* 51), and the greedy search finds a good grammar, not the smallest.

**Remark 3.4 (what node replacement cannot share).** Sierpinski's gasket of depth $d$ is three copies of depth $d - 1$, each pair sharing one corner. A node-replacement step contracts node-*disjoint* occurrences, so a shared corner can belong to only one copy. The gasket therefore keeps an axiom of many nodes: under a fifth of the input at depth 4, but not $O(1)$. The same family with the copies joined by an edge instead of a shared corner is compressed to an axiom of three nodes by four rules used 81, 27, 9 and 3 times (`graphs::bridged_triangles`, with patterns of up to three nodes). Hyperedge-replacement grammars glue a rule's right-hand side at shared external nodes (Habel and Kreowski 1987, *STACS 87*; Drewes, Kreowski and Habel 1997, in Rozenberg (ed.), *Handbook of Graph Grammars*, vol. 1). The gasket is a hyperedge-replacement language, and such grammars are planned.

## 4. Derivation, and Grammars Made Probabilistic

*Modules `derivation` and `stochastic`.*

**Definition 4.1 (expansion).** To expand a nonterminal node $v$ labelled $A$, place $P_A$ by $v$'s frame. Each edge $\{u, v\}$ is given to the node of $P_A$ nearest to $u$: this is the *relaxed embedding* of Fiser et al. (2016), §3.3. It always rewrites, and for a straight-line grammar it differs from Proposition 3.3's replay only where the recorded connection is not the nearest node (`derivation::expand`).

**Definition 4.2 (categories).** Two rules fall in one *category* when their patterns have the same coarse key: the canonical key at a coarser tolerance, with node symbols ignored (`stochastic::shape_key`). This is state merging in the sense of Stolcke and Omohundro (1994), "Inducing probabilistic grammars by Bayesian model merging", *ICGI-94*, with shape as the merging criterion. A category's rules are its alternative productions. A rule whose pattern contains a node of its own category makes the category recursive. On self-similar input, the rules of the successive levels merge into one recursive category, and the grammar generates structures of every size.

**Definition 4.3 (PCFG).** The probability of a production within its category is its relative frequency of use in the derivation (`stochastic::model`). This is the maximum-likelihood estimate given the derivation, and it is consistent, assigning probability one to finite derivations (Chi and Geman 1998, *Comput. Linguist.* 24).

**Definition 4.4 (PCSG).** Two conditions make the grammar sensitive to context (`stochastic::Kind::pcsg`):
- **Parent annotation.** The probability of a production depends on the category of the node's parent in the derivation as well as its own (Johnson 1998, *Comput. Linguist.* 24). The counts are smoothed by adding a constant $\alpha$ to each.
- **Geometric context.** A production whose new edges would cross an edge already drawn is refused, and another is drawn (`SampleOptions::tries` attempts).

The second condition depends on the whole sentential form, not on a bounded neighbourhood. It is a context condition in the sense of graph grammars with application conditions. What it does to the language is not characterised here; the checks compare sampled graphs with the input by their degree distributions (total variation) and their crossings, and the PCSG's samples cross no more than the PCFG's.

## 5. Syntax Analysis of a Geometric Graph

*Module `parsing`.*

**Definition 5.1 (the analyser).** Given a grammar $\mathcal{G}$ and a geometric graph $H$, take the rules in order of creation. For each rule, find every occurrence of its key among the connected node sets of its size whose symbols match (`induction::occurrences`, the routine learning uses). Contract a greedy node-disjoint choice of them, and repeat passes to a fixpoint (`parsing::parse`). The result is a parse forest:
- every node made, with its children and the input edges it swallowed;
- the uses of each rule;
- the share of the input covered;
- whether what remains has the axiom's key (*derived*).

**Proposition 5.2 (invariance).** If $H = g \cdot X$ for $g \in \mathrm{Sim}^+(2)$, the analyses of $H$ and of $X$ make the same choices and the same uses, up to the rounding of Proposition 2.4.

*Proof.* Every comparison the analyser makes is by keys, which are invariant (Proposition 2.4). The greedy order is by node index, which $g$ does not change. ∎

The check parses a copy of a learned structure that has been turned by 37°, scaled by 2.3 and moved: it is covered wholly, derived, with the rules used as often as in learning.

*Caveat.* The analyser is deterministic and greedy. It is not a decision procedure for membership: a graph of the language can be missed if a greedy choice takes the wrong occurrences of an ambiguous pattern. Membership is NP-complete already for regular DNLC grammars (Aalbersberg, Ehrenfeucht and Rozenberg 1986, *Discrete Appl. Math.* 13), so a complete parser is exponential in the worst case. An exact parser for small inputs by backtracking over choices, beside the greedy one, is planned.

## 6. Words of Dimension Two

*Module `simplicial`.*

### 6.1 Letters

**Definition 6.1 (letter).** The triangles up to orientation-preserving similarity and reflection are their angle triples $\alpha \ge \beta \ge \gamma > 0$ with $\alpha + \beta + \gamma = \pi$, a two-dimensional open simplex. At a resolution $\delta$ (in degrees, dividing 180), the *letter* of a triangle rounds $\alpha$ and $\beta$ to multiples of $\delta$ and takes $\gamma$ as what makes 180 (`simplicial::letter`, written for instance `90·45·45`). Handedness, which reflection forgets, is kept where a triangle is glued.

**Definition 6.2 (the complex of a polygon).** A polygon is cut into triangles by its *constrained Delaunay triangulation* (physis's `delaunay_ring`, its decision 0024): ears are clipped, and diagonals are then flipped until each is locally Delaunay (Lawson 1977). This triangulation maximises the smallest angle among all triangulations of the polygon, and it is unique where no four vertices are cocircular (Lee and Lin 1986, *Discrete Comput. Geom.* 1). A polygon's word therefore depends on its geometry alone, not on where ear clipping began. Rings that cross themselves, or have fewer than three distinct corners, are left out and counted (`Complex::skipped`).

### 6.2 The Gluing Grammar

**Definition 6.3.** Let $K$ be a 2-complex and $K^*$ its dual graph, with triangles as nodes and shared edges as arcs. A *derivation* of $K$ is a breadth-first spanning forest of $K^*$ (`simplicial::derive`). Each triangle's edges have *slots* by length: 0 is the longest, opposite the largest angle, then 1, then 2.

The gluing grammar's nonterminals are pairs (letter, entry slot). A production assigns to each other slot either nothing or a triple (letter, slot, handedness): a triangle glued there by one of its own slots, on one side or the other. Slots are chosen independently, each with the relative frequency of its choice in the derivation (`simplicial::learn_gluing`). Growth draws these choices and places each triangle on its parent's edge by its letter's angles (`simplicial::grow`). A triangle that would overlap one already placed is refused: a geometric context condition, as in Definition 4.4.

**Definition 6.4 (entropies).** Let three quantities be per triangle:
- $H(L)$ is the entropy of the letters;
- $H(L_c)$ is the entropy of the glued triangles' letters;
- $H(L_c \mid L_p, s)$ is the entropy of a glued letter given its parent's letter and the slot it is glued to, with its own slot and handedness summed out.

Then $I = H(L_c) - H(L_c \mid L_p, s) \ge 0$ is the information a neighbour carries about a letter, a mutual information (`GluingGrammar::entropy`, `glued_entropy`, `conditional_entropy`).

*Caveat.* These are plug-in estimates from counts, biased low when contexts are many and each is seen a few times (Miller 1955; Paninski 2003, *Neural Comput.* 15). A star of 22 triangles in 6 letters has its conditional entropy estimated at 0, every context being seen once.

**Proposition 6.5 (subcritical growth).** In a derivation of a complex of $T$ triangles in $r$ components, the mean number of triangles glued to a triangle is $(T - r)/T < 1$ (`simplicial::mean_offspring`).

*Proof.* Every triangle except the $r$ roots is glued to exactly one parent. ∎

Read as a branching process, the grammar's patches are finite with probability one: relative-frequency estimates give a consistent grammar (Chi and Geman 1998). For one type, a Galton–Watson process with mean $m = (T - 1)/T$ has expected total progeny $1/(1 - m) = T$. Most draws are nevertheless small, since the law of the total progeny of a near-critical process is heavy-tailed, and refusals for overlap cut branches further. A patch of a given size is therefore drawn by rejection: up to a stated number of draws, keeping the first that reaches the size asked, or else the largest (`grow`'s `at_least`, `attempts`, and `Grown::draws`, `Grown::reached`). That is a draw from the grammar's law conditioned on its size, which the output states.

*Remark.* A letter whose smallest angle rounds to 0 stands for a sliver. Placed with an angle of 0, its two long sides would be parallel. Growth gives such an angle the middle of its class, $\delta/4$, taken from the largest (`simplicial::placement_angles`).

### 6.3 The Replacement Grammar

**Definition 6.6 (moves).** A triangle $\sigma$ rewrites by one of three moves (`simplicial::subdivide`), or stops:
- **Pachner 1–3:** a vertex at the centroid, giving three triangles;
- **midpoint:** the edges' midpoints joined, giving four triangles similar to $\sigma$;
- **Sierpinski:** the three corner triangles of the midpoint move, the middle removed.

With probabilities $(w_\text{stop}, w_{13}, w_\text{mid}, w_\text{S})$ to a maximum depth, this is a PCFG whose words are triangles (`simplicial::replace`). Pachner's move and the midpoint subdivision are subdivisions of $\sigma$: their triangles cover it. Sierpinski's is not: it removes a quarter of the area and opens a hole, and iterated it tends to the gasket.

**Proposition 6.7 (leaf count).** A derivation with $p$ Pachner moves, $q$ midpoint moves and $s$ Sierpinski moves has $1 + 2p + 3q + 2s$ leaves.

*Proof.* Each move replaces one leaf by three, four and three leaves respectively. ∎

### 6.4 Syntax Analysis of Triangles

The analyser of the replacement grammar works from the leaves alone, bottom up (`simplicial::parse_replacement`). In each pass, it supposes each current triangle $t$ to be a child of a move, with one of its vertices $a$ as the parent's corner. That determines the parent:
- for a corner of the midpoint or Sierpinski move, $t = (a, m_1, m_2)$ with $m_1, m_2$ the midpoints of the parent's edges at $a$, so the parent is $(a, 2m_1 - a, 2m_2 - a)$;
- for a child of Pachner's move, $t = (a, b, g)$ with $g$ the parent's centroid, so the parent is $(a, b, 3g - a - b)$.

The parent is accepted when every child of the move is among the current triangles, compared by their vertices to a tolerance. It then replaces them. Where Sierpinski's three corners are found and the midpoint move's medial triangle is present too, the midpoint is taken. Every hypothesis determines its parent and its move, so an accepted group has one reading. The analysis ends with one root when the leaves came from one derivation. The check recovers the root and the count of each move from the leaves of twenty random derivations.

## 7. Placement in the Infrastructure

The grammars name no object in nature, which is why they are a foundation. The infrastructure's domains supply the objects:
- urban science, a street network as a geometric graph and its blocks as faces;
- architecture, footprints as polygons and plans as outlines.

Chemistry's molecular graphs and the cell complexes of districts are next (decision 0003). Hyperedge replacement (Remark 3.4), grammars on 3-complexes whose words are tetrahedra, and an exact parser (§5) are planned there.
