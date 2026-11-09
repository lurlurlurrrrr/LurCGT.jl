# Part 2: Demazure enumeration and string lowering schedules

Status: revised design proposal following completion of part 1. Originally split
from the source review at LurCGT commit `d1f4176`.
No representation code or cache data is changed by these documents.

This is one of three implementation plans:

1. [Unified crystal representation](crystal_basis_irrep_plan.md).
2. [Demazure enumeration and string lowering schedules](crystal_basis_string_plan.md).
3. [Numerical representations, symmetry integration, and caches](crystal_basis_integration_plan.md).

## Scope and dependencies

Enumerate vertices through Demazure subsets, then use string coordinates for a
fixed reduced expression of the longest Weyl element to label crystal vertices,
order columns within each weight, and specify reproducible lowering schedules.
Start with existing SU crystals, including a weight with multiplicity greater
than one, before applying this machinery to the other models from part 1.

Consume `CrystalVertex{S}` values and their weight, raising/lowering, and
epsilon/phi operations from [part 1](crystal_basis_irrep_plan.md). Fundamental
vertices already have tensor-word embeddings; string extraction operates through
these vertex methods without requiring a whole-crystal type or precomputed edge
tables. New spin or exceptional support is not a prerequisite. Retain `Irep` and
its sparse weight-block storage. Changes to ordering and numerical schedules
belong here, after the vertex replacement, and require part 3's exact-arithmetic
and cache activation gates before production use.

## Current lowering-sequence assignment

`get_irep_frommw` in `src/clebsch.jl` currently couples crystal discovery to
numerical vector construction. Starting at the highest vertex, it processes a
queue of newly discovered vertices and tries colors in ascending order.
`det_getvector` checks whether the remaining monochromatic lowering string
contains any undiscovered vertex. If so, the traversal walks that string,
applying `Sl[i]` and `divcfac` at every step and saving each new vertex/vector.

The first discovery selects the stored numerical vector; no explicit lowering
schedule is retained. When a walk passes through an already discovered vertex,
it continues with its current numerical vector rather than switching to the
stored vector at that vertex. Columns within a weight are subsequently filled
by dictionary iteration. Thus both path choice and column order depend on
traversal/storage conventions rather than explicit string coordinates.

Extending monochromatic strings resembles Demazure closure, but this adaptive
queue is not a prescribed sequence of Demazure subsets. Replace that scheduling
policy while retaining the existing vertex operators and numerical route as a
comparison baseline during development.

## Demazure enumeration

Fix the written reduced expression ``w_0=s_{i_1}\cdots s_{i_L}``, where L is the
number of positive roots. Use this convention throughout enumeration, coordinate
extraction, and numerical evaluation; operator products act rightmost first.

For a subset X of the normal highest-weight crystal, define the set closure
from §13.2, equation (13.4):

```math
\mathfrak D_i(X)
 =\{\widetilde f_i^a b : b\in X,\ a\ge0\}\setminus\{0\}.
```

Starting with ``X_0=\{u_\lambda\}``, compute

```math
X_k=\mathfrak D_{i_{L-k+1}}(X_{k-1}),\qquad k=1,\ldots,L.
```

These are the Demazure subsets for the successive suffixes of the written
reduced word. Theorem 13.5 guarantees that the final set is the whole crystal,
that the subset for a Weyl element is independent of its reduced expression,
and that each intermediate subset intersects any color string in nothing,
its highest vertex alone, or the entire string.

Implement each stage by extending color strings from the previous stage and
merging vertices by tensor-word equality. Include the previous set (exponent
zero); processing only the vertices newly added at the preceding stage is not
sufficient. Keep numerical vectors out of enumeration. Never mutate a word
already published as a dictionary or set key. Initially use straightforward
closure and deduplication; use the string-intersection property to avoid
redundant walks only after validating that optimization.

Demazure closure enumerates vertices but does not select a unique path to each
vertex. Do not use first discovery to assign the final schedule. Extract the
string coordinates below after enumeration. This is the set operator of §13.2,
not the signed linear crystal Demazure operator of §13.3; the latter is not
needed for vertex enumeration.

## Reduced words and string patterns

Fix one reduced expression ``w_0=s_{i_1}\cdots s_{i_L}``, where L is the number
of positive roots. Given a vertex b of a normal highest-weight crystal, set

```math
b_0=b,\qquad a_k=\varepsilon_{i_k}(b_{k-1}),\qquad
b_k=\widetilde e_{i_k}^{a_k}b_{k-1}.
```

Then ``b_L=u_\lambda`` and

```math
b=\widetilde f_{i_1}^{a_1}\cdots
  \widetilde f_{i_L}^{a_L}u_\lambda.
```

Thus the chronological lowering schedule is **k=L,L−1,…,1**, with node
``i_k`` applied ``a_k`` times. This distinction between written operator order
and execution order must appear in the implementation and tests.

Section 11.1 defines these coordinates and proves the type A result; the
all-normal-crystals justification is Corollary 13.8, not Theorem 11.6 alone.
The string tuple identifies a vertex only relative to the chosen word and
highest-weight component. It is not a word-independent canonical label.

Suggested implementation:

1. Generate a deterministic reduced word using exact root reflections and a
   fixed smallest-node choice among length-increasing simple reflections.
   Check that its inversion count/length equals the number of positive roots
   and that the resulting element maps every positive root to a negative root.
   Specify whether reflections are multiplied on the left or right, and return
   the word in the written-product convention above. Persist the chosen word
   as part of the basis convention.
2. Enumerate vertices through the Demazure stages above using the reversible
   vertex operators from part 1. Compare the resulting vertex set with the
   existing traversal during development. Require inverse-edge and color-string
   checks before extracting coordinates; no stored graph is required.
3. Compute the above string tuple for each vertex. Verify termination at the
   unique highest vertex, reconstruction, and uniqueness of tuples.
4. Group vertices by weight, retain the established weight order where valid,
   and sort by lexicographic string tuple within each weight. Dictionary
   iteration and graph-discovery order must not select basis columns.
5. Evaluate numerical vectors from the extracted schedules in the predetermined
   weight/column slots, subject to the numerical gates below. Enumeration paths
   must not determine numerical vectors.
6. Store exponents in a concrete integer matrix or compressed nonzero runs.
   Do not expand long powers into repeated operator indices. Once enumeration
   is complete, string extraction and independent vector evaluation can be
   parallelized with predetermined output slots.

Use existing structures and small functions with concrete containers; no new
whole-crystal wrapper or type hierarchy is needed. Document each new named
function, including ordering and mutation semantics. Demazure stages depend on
earlier stages; parallel work within a stage must use independent candidate
buffers followed by a merge, rather than concurrent mutation of shared sets.

Do not enumerate an arbitrary rectangular box of exponents and assume it is
valid. Strings satisfy additional constraints. Initially derive them from the
verified crystal; string-cone/polytope enumeration is a separate optimization.
Do not confuse these strings with Lusztig/PBW data: their coordinate maps differ.

### A crystal path is not a numerical representation matrix

Kashiwara operators ``\widetilde f_i`` are combinatorial maps. Ordinary lowering
generators ``F_i`` can produce linear combinations and have normalization
coefficients. A crystal graph alone does not determine the existing `Sl` and
`innerprod` data. The explicit coefficient 2 and zero-weight metric in `G2.jl`
already illustrate why general edge adjacency is insufficient.

Use the string schedule to propose numerical vectors
``F_{i_1}^{(a_1)}\cdots F_{i_L}^{(a_L)}v_\lambda``, with divided powers
``F_i^{(a)}=F_i^a/a!`` and an explicit integral rescaling convention. This does
not assert that these are global/canonical basis vectors. The combinatorial
reconstruction theorem by itself does not prove independence after replacing
Kashiwara operators by classical generators. Before making this the production
basis, supply an appropriate monomial-basis/triangularity justification and
verify exact rank and expected weight multiplicity for each constructed block.
A failed rank check must stop construction, not cache an incomplete irrep.
Retain the existing numerical route during this development gate.

Preserve Gram matrices and solve for lowering coefficients exactly. Primitive
rescaling by gcd must be deterministic, including the sign, and recorded when
reusing partial paths. Such rescaling can change integrality of the operator
matrices: `Irep.Sl` currently has `Int` entries even when the metric type is
`BigInt`. Audit denominator clearing and bounded integer conversion explicitly;
do not hide rational or overflowing coefficients in the existing storage.

## Acceptance and handoff

Follow the implementation sequence above, starting with an A2 repeated-weight
example. Check reduced-word length/action, string round trips and tuple
uniqueness. Use an existing relevant test file and small regression cases.
Also check intermediate Demazure string intersections and final-set equality
with the current traversal. For small examples, compare the sets from distinct
reduced expressions of the same Weyl element. Coordinates and numerical bases
may differ between words; only the Demazure sets must agree. Include a
non-palindromic reduced word to expose accidental reversal of execution order.
For a fixed word, different discovery orders must produce identical ordered data.
Compare dimensions, exact ranks, generator relations, and metrics with current
results before changing numerical construction; compare different bases through
intertwining transformations. Hand the validated word, ordering, rescaling rules,
and schedules to part 3 for production activation and cache versioning.

Use `@code_warntype` or JET at string extraction and payload barriers, and measure
allocations/runtime against existing benchmarks where applicable. Apply the
[shared implementation policy](crystal_basis_integration_plan.md#shared-implementation-policy).

## References

Bump and Schilling, *Crystal Bases: Representations and Combinatorics* (2017),
§11.1 (printed pp. 144 onward), §13.2 equation (13.4) and Theorem 13.5
(printed p. 174), and Corollary 13.8 (printed p. 176). The supplied source is
`/home/lurlurlur/Downloads/Daniel Bump, Anne Schilling - Crystal Bases_ Representations and Combinatorics-World Scientific (2017).pdf`. See
[part 1's provenance](crystal_basis_irrep_plan.md#references-and-provenance)
for the user-supplied PDF and original review details.
