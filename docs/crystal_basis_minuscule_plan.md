# Part 3a: Minuscule irreps and fundamental construction recipes

Status: implementation and validation specification. See the [part 3 overview](crystal_basis_integration_plan.md)
for dependencies and shared requirements.

The [recipe-struct cleanup](crystal_basis_recipe_structs_plan.md) is implemented:
per-symmetry construction and validation catalogues drive recursive requests.

## Scope

Replace the alternating-product irrep backend with direct minuscule seeds and
selected tensor-product components for SU(N), Sp(2N), SO(2N), E6, and E7.
Virtual crystals and adjoint seed backends belong to later parts; SO(2N+1),
F4, G2, and E8 are outside this implementation scope.

Retain `Irep`, `CrystalVertex`, exact arithmetic, and part 2's Demazure
traversal. Implement seed construction and generalize `get_irep_frommw` to
select a requested component of two known irreps, then dispatch fundamental
construction through fixed recipes. A fundamental need not be the maximal
component of its source product.

## First example: SO(8)

| Target q-label | Dimension | Hard-coded construction |
| --- | ---: | --- |
| `(1, 0, 0, 0)` | 8 | Minuscule vector seed |
| `(0, 1, 0, 0)` | 28 | Select this component of `(1, 0, 0, 0)` tensor `(1, 0, 0, 0)` |
| `(0, 0, 1, 0)` | 8 | Minuscule half-spin seed |
| `(0, 0, 0, 1)` | 8 | Minuscule half-spin seed of the other chirality |

The largest component of the vector square has label `(2, 0, 0, 0)` and dimension
35. The target 28-dimensional component therefore needs its own highest-weight
vector inside the ambient product weight space. Passing a different `mw` while
keeping the scalar initial state `RT[1]` does not perform that selection.
The vector square decomposes as 35 + 28 + 1.

## Conventions and label adapters

Use the book's node numbering: the E6/E7 chain is `1–3–4–5–6[–7]`,
with node 2 attached to 4; the D fork has terminal nodes `r-1,r`.

Record the Cartan matrix with the explicit convention
``A_{ij}=\langle\alpha_i^\vee,\alpha_j\rangle``. In fundamental-weight
coordinates, reflection in node i sends a weight x to
``x-x_i A_{:,i}``; lowering subtracts ``A_{:,i}``.
Provide explicit conversions between these coordinates and existing z-labels.
Do not infer root lengths from the current storage coordinates. In particular,
SO uses doubled, reversely arranged orthogonal coordinates; preserve its external
q-label meaning and establish the node permutation by checking each `getdz`.

Spin modules represent ``\mathfrak{so}(N)`` and ``\mathrm{Spin}(N)``;
they do not generally descend to ``\mathrm{SO}(N)``. For this extension, `SO{N}` must support the requested
Lie-algebra/Spin representations, including the spin labels in the example.
Document that meaning and update existing parity restrictions consistently. Handle C1=A1, D3=A3, and D2=A1⊕A1 separately; SO(2) is Abelian.

## Minuscule seeds

For a minuscule highest weight, enumerate its Weyl orbit. Each weight occurs once,
and ``\widetilde f_i(\mu)=\mu-\alpha_i`` exactly when
``\langle\alpha_i^\vee,\mu\rangle=1``. Raising uses pairing −1.
This gives a uniform construction for all minuscule seeds, including spin.
In the minuscule case, unit-coefficient edge maps also realize Chevalley
raising/lowering generators with diagonal coroot action; adapt any existing
normalization explicitly. This special numerical shortcut does not extend to
arbitrary crystals. See the [Sage minuscule representation construction](https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/quantum_groups/representations.html)
and [Geck's minuscule construction](https://arxiv.org/abs/1603.07179).

### Integer symbols and vertex storage

Keep `CrystalVertex{S}.word::Vector{Int}` for every family. Assign each vertex
of a minuscule seed one integer symbol `a`; its crystal vertex is
`CrystalVertex{S}([a])`. Weights, sign tuples, and orbit coordinates are metadata
associated with that symbol, not elements of the stored word. Do not introduce
a new vertex type or a vector of structured symbols.

Use deterministic, disjoint symbol ranges for distinct minuscule seeds of the
same symmetry, preserving the existing vector alphabet. A symbol must identify
both its seed and its orbit vertex without extra per-word metadata, even when
different seeds contain the same weight. Assign ranges independently of request
order; check integer overflow. Within each seed use a fixed orbit traversal
(highest weight first, lowering colors in ascending order, first discovery),
or the explicit compatible encoding specified below.

Provide symbol-to-weight and single-symbol raising/lowering access through the
existing crystal interface, using concrete integer tables or family-specific
integer decoding. Local epsilon and phi are 0 or 1 for minuscule symbols.
Missing edges follow the existing interface's no-action convention.

Tensor products concatenate these integer words in the prescribed factor
order: `[a]` tensor `[b]` gives `[a, b]`. Selected components retain their
words in the original minuscule alphabet; do not replace a nonminuscule
component vertex with a new atomic symbol. This keeps the existing word
signature rule applicable to staged fundamental constructions.

### D-family spin words

The notation `s(...)` below denotes an encoded integer, not a stored object.
Keep the existing vector alphabet compatible. Encode each spin symbol
``s(\sigma_1,\ldots,\sigma_r)`` as one integer in a disjoint range, using
sign bits with overflow checks. Its weight is
``\frac12\sum_j\sigma_j e_j``. The two sign-parity classes give the two
half-spin modules, each of dimension ``2^{r-1}``:

| Fundamental node | Highest word |
| --- | --- |
| r−1 | `[s(+,...,+,-)]` |
| r | `[s(+,...,+)]` |

For i<r, lowering changes signs `(+,-)` at positions i,i+1 to `(-,+)`.
At node r it changes the final pair `(+,+)` to `(-,-)`. Otherwise lowering
is undefined; raising reverses these changes. Local epsilon/phi are 0 or 1.
The stored SO z-label is ``(\sigma_r,\ldots,\sigma_1)``; validate transitions
against `getdz` and use the node adapter to identify the terminal q-labels.
No separate chirality identifier is needed.

## Fundamental recipes

Write ``C_\lambda(X)`` for one connected component of crystal X whose highest
vertex has weight λ. It never means that all of X is irreducible. Specify tensor
factor order and the signature convention; the book and other software may use
opposite conventions. Check a small tensor example before copying highest tuples.

### Construction order

Write `F_k` for the fundamental irrep with unit Dynkin label at node k.
`select(F_k; X, Y)` means construct that target component of the ordered product
`X tensor Y`, using the preceding-multiplet algorithm below. All table entries
use deterministic copy 1; determine the highest crystal vertex in the fixed
local tensor convention. This selection must not require a preexisting target.

| Family | Direct minuscule seeds | Binary construction of remaining fundamentals |
| --- | --- | --- |
| A_r, `SU{r+1}` | Every `F_k`, 1≤k≤r | None: construct each requested fundamental directly from its Weyl orbit |
| C_r, `Sp{2r}` | `F_1` | For k=2,…,r: `select(F_k; F_{k-1}, F_1)` |
| D_r, `SO{2r}`, r≥4 | `F_1`, `F_{r-1}`, `F_r` | For k=2,…,r−2: `select(F_k; F_{k-1}, F_1)` |
| E6 | `F_1` (27), `F_6` (27) | `F_2 ← F_1 tensor F_6`; `F_3 ← F_1 tensor F_1`; `F_5 ← F_6 tensor F_6`; `F_4 ← F_3 tensor F_1` |
| E7 | `F_7` (56) | `F_1 ← F_7 tensor F_7`; `F_6 ← F_7 tensor F_7`; `F_2 ← F_1 tensor F_7`; `F_5 ← F_6 tensor F_7`; `F_3 ← F_1 tensor F_1`; `F_4 ← F_5 tensor F_7` |

Each arrow selects only the named target, not the entire product. In particular,
the E7 recipe for `F_3` uses the adjoint square instead of constructing a full
fourfold minuscule product. Validate the occurrence of every target and the
node-label adapters before enabling these recipes.

Independent dimension targets:

- A_r: ``dim F_k=\binom{r+1}{k}``.
- C_r: ``dim F_k=\binom{2r}{k}-\binom{2r}{k-2}``.
- D_r: ``dim F_k=\binom{2r}{k}`` for k≤r−2; both half-spin dimensions are
  ``2^{r-1}``.
- E6, nodes 1,…,6: 27, 78, 351, 2,925, 351, 27.
- E7, nodes 1,…,7: 133, 912, 8,645, 365,750, 27,664, 1,539, 56.

A binomial coefficient with negative lower index is zero. These dimension
formulas do not prescribe an exterior-product implementation.

For C1 use the A1 minuscule construction. For D3 use the A3 adapter and its
three minuscule fundamentals; for D2 use the two A1 factors. Record the explicit
node permutation in each adapter and verify it against `getdz`. SO(2) remains
on its Abelian path.

### Hard-coded recipes and recursive dependencies

Each symmetry provides a `Vector{FundamentalRecipe}` indexed by fundamental
node. Concrete minuscule, adjoint, and product recipe subtypes specify:

- target fundamental q-label and expected dimension;
- seed kind, or two ordered source labels for component extraction;
- selected component/copy and its highest crystal vertex or a fixed rule to
  obtain it in the library's tensor convention.

Use the binary dependencies above directly; do not search tensor powers for
a construction at runtime. Where the book gives a highest tuple, translate
it to the local factor order and verify its weight and vanishing raising actions.

Start work on a new symmetry with `initialize_symmetry(S, RT; verbose)`: call
`isvalidsymm`, reject invalid symmetries before storage access, and save only
the trivial irrep. This initializer does not enumerate fundamentals.

The bulk `get_fundamental_ireps` entry point is removed. On a
`getNsave_irep` cache miss, execute only the requested fundamental's hard-coded
recipe and its dependencies. Keep trivial initialization separate from this
per-target construction; do not trigger bulk generation on a cache miss.

Resolve dependencies recursively and cache each completed irrep. Distinguish
actual unit-label fundamentals from convenient doubled-spin bootstrap modules.
A fundamental cache miss must execute its own recipe rather than regenerate
all fundamentals and hope the requested one was saved. Keep recipes queryable
without enumeration; the exceptional dimension tables are not an instruction
to allocate every fundamental eagerly.

## Generalized tensor-component construction

Modify `get_irep_frommw` to accept a product representation and the maximal
weight of the target irrep, returning the selected `Irep`. With alternating
products removed, both `inner_prod1` and `inner_prod2` are concretely typed
`Dict{NTuple{NZ, Int}, Matrix{RT}}`; product `weight_info` is likewise a
`Dict{NTuple{NZ, Int}, WeightInfo{NZ}}`, with no `nothing` branch.
The maximal Cartan component remains the special case used by dominant-weight
recursion.

### 1. Build the ambient product

Reuse `tensor_prod_reps` and its sparse lowering blocks, weight ranges, and
`WeightInfo`. Preserve ordered input factors. Compute the target weight using
`qlab2mwz`; do not choose it from the largest ambient weight.

### 2. Select the maximal-weight state through preceding multiplets

Follow the three-leg CG construction: process product highest weights in
descending order, construct all larger irreps and their copies, and retain
their embedded multiplets and exact metric data. At the target weight, project
candidate product-basis vectors orthogonally to those already-found states,
as in `get_mwstate`. This yields the target maximal-weight state without a
direct raising-kernel solve.

Use deterministic candidate ordering, exact projection, and sign conventions.
For repeated targets, project out preceding copies as well and select the
requested copy deterministically. Reject an absent target or an out-of-range
copy. A crystal highest tuple selects a combinatorial component for traversal;
it is not generally the numerical highest-weight pure tensor, nor does it
provide a canonical matching to numerical copies.

Construct the preceding irreps from the same ambient product using the basis
builder below. Do not load an unknown output irrep to generate its multiplet:
that would make fundamental bootstrap circular.

### 3. Generate the selected component without an output irrep

Initialize the target component's highest crystal vertex together with the
full ambient numerical highest vector. Replace the assumptions that the ambient
highest-weight space is one-dimensional and that the initial vector is `RT[1]`.
The initial weight-basis column must have the actual ambient sector dimension.

Run part 2's Demazure traversal from that crystal vertex. Apply the ambient
lowering matrices to the selected numerical vector, retain first discovery,
and stop a string at an already stored vertex. Use the same insertion order
for numerical weight-basis columns. No string-coordinate extraction or separate
rank check is required.

`generate_every_CGT` currently calls `getNsave_irep` for the output, and
`get_multiplet`/`get_wspace!` use that output's lowering matrices. Do not use
those paths to bootstrap the unknown irrep: that would be circular. Reuse
independent lower-level operations and the generalized Demazure basis builder.

### 4. Assemble the new irrep

Reuse `get_irep_frommw`'s Gram-matrix construction and exact coefficient solves
to express ambient lowering in the generated basis. Normalize the completed target basis uniformly: if its maximal-weight vector
has squared norm h in the product metric, every basis vector receives the
same factor ``1/\sqrt{h}``. Thus every Gram block is divided by h, inverse
Gram blocks are multiplied by h, and lowering matrices are unchanged. Do not
normalize weight spaces or basis columns independently. Represent this common
normalization through exact Gram-block division rather than introducing
floating-point square roots into the basis construction. If this division or
construction of integer basis/generator data requires a nonintegral rational
value, raise an error; do not truncate or promote the storage to rationals.
The existing rational factors representing inverse metrics remain supported.

Construct inverse metrics consistently with the normalized Gram blocks.
Preserve primitive rescaling, sign conventions, integrality, and overflow
handling; let existing numerical failures propagate. Save only the successfully
constructed, normalized `Irep`.

Return an `Irep`; a full CG tensor is not required for bootstrap.

## Integration sequence

1. Implement and integrate direct minuscule seeds, including the SO(8) vector and both half-spin seeds.
2. Generalize the selected-component builder and extract SO(8) `(0, 1, 0, 0)`
   from two defining irreps without preloading the target.
3. Encode the in-scope classical and E6/E7 fundamental recipes with acyclic
   dependencies, including staged binary forms of the tensor powers.
4. Replace the blanket exterior-power bootstrap with per-target recipe dispatch.
   Remove the alternating-product irrep backend after migrating its in-scope
   callers. Audit out-of-scope callers before deleting shared code. Existing
   realizations may serve as migration comparisons, not retained backends.
5. Keep general dominant-weight recursion through the maximal Cartan component
   of known fundamentals. Update label adapters, duality, and cache conventions.

## Canonical consistency of incidental irreps

Validate every nonfundamental irrep generated while obtaining fundamentals
from product representations, including every multiplicity copy. Here the
validation scope is all such nonfundamental results, not irreps ordered by
dimension or enumeration of arbitrary tensor-product cycles.

The canonical path uses the fixed fundamental recipes and the fixed
dominant-weight split/maximal-component construction for nonfundamental labels.
Factor order and split selection are deterministic.

The highest-weight crystal and fixed Demazure traversal determine basis order.
Use the same operator sequence, frozen source counts, first-discovery rule,
and divided-power factors in every embedding. Do not order traversal by
ambient integer symbol values or dictionary iteration. Uniform highest-weight
normalization then gives the same intrinsic basis convention. No separate
basis recipe or change-of-basis stage is needed.

Validation proceeds in two phases:

1. Construct the fundamentals needed for the validation run and save them in
   the database. Record every incidental nonfundamental `Irep`, including all
   copies, with its source pair and target label. Do not save these incidental
   candidates in the canonical irrep cache or database.
2. After the fundamentals are saved, obtain each recorded target through its
   canonical path using those saved fundamentals, and compare every candidate
   against that reference. Reuse a canonical reference for repeated targets.
   Reference construction does not recursively launch another validation run.

Require the same `Irep` type and equality of every representation field: `Sl`, `Sz`,
`innerprod`, `inv_innerprod`, `qlabel`, and `dimension`. Exclude `size_byte`. Compare
field values, not object identity. In particular, compare the stored inverse
metric matrices and factors themselves, not only their represented product.
Reduce inverse-metric matrix entries to relatively prime integers, absorbing
the common factor into the rational scalar. Only inverse metrics may be
rational. Existing databases are removed before validation; old entries are
not canonical references. Hard-coded recipes define acyclic dependencies,
resolved recursively; all required fundamentals are saved before canonical comparisons.
Do not accept an intertwiner as a substitute for field equality. Different
multiplicity copies have different ambient embeddings but must produce the
same `Irep` fields.

For example, extracting the SO(8) adjoint `(0, 1, 0, 0)` from the vector square
first produces `(2, 0, 0, 0)`. Compare that incidental 35-dimensional irrep
with the canonical maximal-component construction of `(2, 0, 0, 0)` from two
vector fundamentals. Apply the same check to all nonfundamental components arising
in the C, D, E6, and E7 fundamental recipes.

Validation bookkeeping stays outside the production hot path. Report deferred
or skipped targets explicitly; fundamental-bootstrap validation is complete
only when every recorded candidate has passed. Avoid retaining full ambient
embeddings after extraction when the normalized irrep data suffices.

## Validation by independent construction paths

Use ordinary per-family functions returning concrete tuples. Each validation
rule contains the target q-label, two ordered source pairs, deterministic copy
indices, and rank restrictions. Keep construction recipes and validation rules
separate so a failed comparison identifies both paths.

Let `A` denote the adjoint and `F_k*` the dual fundamental. The adjoint occurs
in `R tensor R*` for every nontrivial irrep R of a simple Lie algebra: the Lie
algebra action embeds the adjoint in `End(R)`. Hard-code the following baseline
and alternative paths, then extend the alternative column to every fundamental
not already covered by the baseline.

| Family | Target A | Baseline sources | Alternative sources |
| --- | --- | --- | --- |
| A_r, r≥3 | ``\omega_1+\omega_r`` | `F_1, F_r` | `F_k, F_{r+1-k}`, k=2,…,r−1 |
| C_r, r≥2 | ``2\omega_1`` | `F_1, F_1` | `F_k, F_k`, k=2,…,r |
| D_r, r≥4 | ``\omega_2`` | `F_1, F_1` | `F_k, F_k*`, k=2,…,r |
| E6 | ``\omega_2`` (78) | `F_1, F_6` | `F_3, F_5`; `F_2, F_2`; `F_4, F_4` |
| E7 | ``\omega_1`` (133) | `F_7, F_7` | `F_k, F_k`, k=1,…,6 |

Encode duality explicitly: A reverses nodes; C fundamentals are self-dual;
D half-spin nodes are exchanged for odd r and self-dual for even r, while
other D fundamentals are self-dual; E6 exchanges 1↔6 and 3↔5, fixing 2 and 4;
all E7 fundamentals are self-dual. Avoid duplicate comparisons related only
by reversing the same source pair. Check dual seed weight sets separately.

For A1 and A2, compare the baseline adjoint from `F_1 tensor F_1*` with the
adjoint from `A tensor A`. Apply these rules through the C1/D3 adapters. For
D2 validate each A1 factor separately, with the other factor trivial; do not
apply the simple-algebra rule to a single combined adjoint.

These are construction-consistency checks. Mathematical correctness is also
covered by existing tests, including the Sp(6) irrep dimension test, and the
separate large-scale tests. Do not duplicate those suites here. Performance
policy for large exceptional validation runs is deferred.

For SO(8), construct `(0, 0, 2, 0)` independently from:

- `(0, 0, 1, 0)` tensor `(0, 0, 1, 0)`;
- `(0, 1, 0, 0)` tensor `(0, 1, 0, 0)`.

Both contain this 35-dimensional irrep. Repeat for the other chirality,
`(0, 0, 0, 2)`, using its half-spin square and the adjoint square.
Vector tensor adjoint is not an alternative for these targets:
`8 tensor 28 = 8 + 56 + 160`.

Build both results independently without allowing the target cache to return
the first result for the second path. Use the same Demazure and normalization
conventions and require equality of all `Irep` fields, as above.

## Completion

Use small cases in existing test files. Start with SO(8): verify all three
8-dimensional minuscule seeds and construction of the 28-dimensional
target from the vector square with no preexisting target cache entry. Check the
highest vector is killed by raising operators, the resulting dimensions and
weights, generator relations, and metric compatibility. Verify unit norm of
the maximal-weight state and uniform rescaling of all Gram blocks. Compare the target with
the existing exterior-square realization through an intertwining transformation
during migration, before removing that backend.

Check that seed vertices have one-element `Vector{Int}` words, seed symbol
ranges do not collide and are independent of request order, and tensor words
preserve factor order and their weights under raising/lowering.
Run canonical-identity checks for every incidental nonfundamental irrep and copy generated during
fundamental construction, then the independent-path checks above.
Cover one repeated-copy target and an absent-target request. Check deterministic
copy selection, correct recursive recipe dispatch, and cache reload. Exercise
E6/E7 minuscule seeds and representative nonmaximal fundamental components;
verify removal of the alternating-product irrep backend and its callers.
Retain existing CG/recoupling checks, type-stability inspection, and relevant
allocation measurements. No new rank-check gate is required.

## Implementation entry points

- `initialize_symmetry(S, RT)` validates S and saves only the trivial irrep.
- `fundamental_recipes(S)` returns recipe structs indexed by node;
  `getNsave_irep` executes them and recursively requests ordered sources.
- `get_minuscule_irep(S, RT, k)` constructs a seed without database access.
- `get_irep_frommw` performs exact preceding-multiplet extraction; optional
  `copy_index` selects a repeated copy and `candidates` records incidental
  nonfundamental results outside the canonical cache. Optional `copies` collects
  all occurrences of the target.
- `validate_fundamental_ireps(S, RT)` saves required fundamentals, then
  compares every recorded copy with its canonical result.
- `irrep_validation_rules(S)` returns `Vector{IrrepValidationRule}`;
  `validate_fundamental_products(S, RT; rules=...)` checks every target copy
  and returns the number of copies checked.
- `same_irep_fields(a, b)` compares every representation field except
  `size_byte`; it does not perform a change of basis.

The alternating-product backend and bulk `get_fundamental_ireps` entry point
are removed. Nontrivial B/G2 irrep requests are unavailable until their
virtual-crystal backend is implemented; there is no legacy seed fallback.
The `minuscule-v2` cache namespace prevents reuse of old basis data.

## Reference

Bump and Schilling, *Crystal Bases: Representations and Combinatorics* (2017),
user-supplied PDF at
`/home/lurlurlur/Downloads/Daniel Bump, Anne Schilling - Crystal Bases_ Representations and Combinatorics-World Scientific (2017).pdf`.

Relevant passages: §5.4 and §5.6, Tables 5.2–5.6.
