# Part 3a: Minuscule irreps and fundamental construction recipes

Status: design proposal. See the [part 3 overview](crystal_basis_integration_plan.md)
for dependencies and shared requirements.

## Goal and implementation order

Hard-code how every fundamental irrep is obtained, following the construction
catalogue in Bump and Schilling, Tables 5.2–5.6 and the classical constructions
of §5.4. Each entry selects an existing seed, a minuscule seed, or a specified
component of an ordered tensor product. Do not infer fundamentals by assuming
they are the largest component of exterior or tensor powers.

Implement two building blocks:

1. Construct a minuscule `Irep` directly from its Weyl orbit, crystal edges,
   numerical generators, and metric.
2. Generalize `get_irep_frommw` to construct a requested irrep from two known
   input irreps, even when the requested highest weight is not maximal in their
   tensor product. Follow the highest-weight selection and multiplet workflow
   used for three-leg Clebsch–Gordan coefficients, without requiring an existing
   output irrep.

Then wire the hard-coded fundamental recipes into `get_fundamental_ireps` and
`getNsave_irep`. Retain `Irep`, `CrystalVertex`, exact arithmetic, and part 2's
Demazure traversal with frozen source counts and stopping at known vertices.
E6/E7 minuscule seeds belong here; F4's seed comes from part 3b and E8's seed
from part 3c. Register their recipes now and enable them as the seeds become ready.

## First example: SO(7)

| Target q-label | Dimension | Hard-coded construction |
| --- | ---: | --- |
| `(1, 0, 0)` | 7 | Existing defining representation |
| `(0, 1, 0)` | 21 | Select this component of `(1, 0, 0)` tensor `(1, 0, 0)` |
| `(0, 0, 1)` | 8 | Minuscule spin construction |

The largest component of the vector square has label `(2, 0, 0)` and dimension
27. The target 21-dimensional component therefore needs its own highest-weight
vector inside the ambient product weight space. Passing a different `mw` while
keeping the scalar initial state `RT[1]` does not perform that selection.

## Conventions and label adapters

Use the book's node numbering for the catalogue below. For E6, E7, E8 the chain
is `1–3–4–5–6[–7[–8]]`, with node 2 attached to 4. For D, the fork has terminal
nodes `r-1,r`. For F4, nodes 1,2 are long and 3,4 short; for G2, node 1 is short.

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
Document that meaning and update existing parity restrictions consistently. Treat B1=A1, C1=A1, D3=A3, D2=A1⊕A1
separately; SO(2) is Abelian rather than a member of the simple B/D catalogue.

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

For B_r spin, the weights are all ``\frac12(\pm e_1\pm\cdots\pm e_r)``,
giving dimension ``2^r``. For D_r, the two sign-parity classes give the two
half-spin modules of dimension ``2^{r-1}``. In the existing doubled SO storage,
these coordinates are integers ±1. Use the weight/node adapter to associate
parity with the correct terminal q-label.

## Spin vertices and highest words

Extend the existing word-only vertex representation with these symbols in this
part. Keep the current vector alphabet and its integer encoding compatible.

| Family | Fundamental node | Highest word |
| --- | --- | --- |
| B_r | r | `[s(+,...,+)]` |
| D_r | r−1 | `[s(+,...,+,-)]` |
| D_r | r | `[s(+,...,+)]` |

Each `s(...)` is one atomic symbol, encoded by one integer. These are true spin
fundamentals, distinct from the existing vector-derived doubled-spin labels.

For spin symbols ``s(\sigma_1,\ldots,\sigma_r)``, with each sign ±1, the
orthogonal-coordinate weight is ``\frac12\sum_j\sigma_j e_j``. Encode them
in a disjoint range beyond the existing vector alphabet, using a deterministic
sign-bit encoding with explicit integer-overflow checks. Symbols in both D
chiralities have distinct sign patterns, so no chirality identifier is needed.

- In B and D, for i<r, ``f_i`` changes adjacent signs `(+,-)` at positions
  i,i+1 into `(-,+)` and is undefined otherwise.
- In B, ``f_r`` changes the final `+` to `-`. All ``2^r`` sign patterns occur.
- In D, ``f_r`` changes the final pair `(+,+)` to `(-,-)`. Even and odd numbers
  of minus signs give the two separate components of size ``2^{r-1}``.
- Raising reverses these changes. Local spin epsilon/phi are 0 or 1. Convert
  weights to the existing doubled, reversed SO z-coordinates: the stored
  tuple is ``(\sigma_r,\ldots,\sigma_1)``. Validate transitions against `getdz`.

## Complete fundamental construction catalogue

Write ``C_\lambda(X)`` for one connected component of crystal X whose highest
vertex has weight λ. It never means that all of X is irreducible. Specify tensor
factor order and the signature convention; the book and other software may use
opposite conventions. Check a small tensor example before copying highest tuples.

### Classical families

| Type | Every fundamental crystal | Numerical realization |
| --- | --- | --- |
| A_r (`SU{r+1}`) | ``B(\omega_k)``, 1≤k≤r, are all minuscule; equivalently ``C_{\omega_k}(V^{\otimes k})`` | ``\bigwedge^k V``, dimension ``\binom{r+1}{k}`` |
| B_r (`SO{2r+1}` at Lie-algebra level) | Vector tensor components for 1≤k<r; minuscule spin for k=r | ``\bigwedge^k V`` for k<r, dimension ``\binom{2r+1}{k}``; spin dimension ``2^r`` |
| C_r (`Sp{2r}`) | Minuscule vector seed; ``C_{\omega_k}(V^{\otimes k})`` for all 1≤k≤r | Primitive part of ``\bigwedge^k V`` (kernel of symplectic contraction), dimension ``\binom{2r}{k}-\binom{2r}{k-2}`` |
| D_r (`SO{2r}` at Lie-algebra level) | Minuscule vector and two half-spin seeds; vector tensor components for 1≤k≤r−2 | ``\bigwedge^k V`` for k≤r−2, dimension ``\binom{2r}{k}``; both spin dimensions ``2^{r-1}`` |

A binomial coefficient with negative lower index is zero. The B top exterior
power has highest weight ``2\omega_r``; the D top exterior power splits into
``2\omega_{r-1}`` and ``2\omega_r``. These are useful existing seeds, but not
the spin fundamental modules. Section 5.4, Proposition 5.11 explains this gap.

### Exceptional families

The following is a dependency catalogue derived from §5.6, Tables 5.2–5.6.
Each tensor entry means the selected component, not the full tensor space.
Dimensions are included as independent implementation acceptance targets.

| Type | Node → dimension | Construction |
| --- | --- | --- |
| E6 | 1→27; 6→27 | Two minuscule seeds M=``B(\omega_1)``, M*=``B(\omega_6)`` |
| E6 | 2→78 | ``C_{\omega_2}(M\otimes M^*)`` |
| E6 | 3→351; 5→351 | ``C_{\omega_3}(M^{\otimes2})``; ``C_{\omega_5}((M^*)^{\otimes2})`` |
| E6 | 4→2,925 | ``C_{\omega_4}(M^{\otimes3})`` |
| E7 | 7→56 | Minuscule seed M=``B(\omega_7)`` |
| E7 | 1→133; 6→1,539 | Respective components of ``M^{\otimes2}`` |
| E7 | 2→912; 5→27,664 | Respective components of ``M^{\otimes3}`` |
| E7 | 3→8,645; 4→365,750 | Respective components of ``M^{\otimes4}`` |
| F4 | 4→26 | Virtual E6 seed V from part 3b |
| F4 | 1→52; 3→273 | Respective components of ``V^{\otimes2}`` |
| F4 | 2→1,274 | ``C_{\omega_2}(V^{\otimes3})`` |
| G2 | 1→7 | Existing numerical seed; virtual D4 crystal as a uniform construction/check |
| G2 | 2→14 | ``C_{\omega_2}(B(\omega_1)^{\otimes2})``; use the existing alternating-square realization as a comparison |

For E8, register the incremental dependencies documented in
[part 3c](crystal_basis_adjoint_plan.md); their seed is unavailable until that part.

### Hard-coded recipes and dependency order

Implement per-family/per-node methods using ordinary functions and concrete
labels/tuples. Each recipe specifies:

- target fundamental q-label and expected dimension;
- seed kind, or two ordered source labels for component extraction;
- selected component/copy and its highest crystal vertex or a fixed rule to
  obtain it in the library's tensor convention.

Convert the table's tensor powers into explicit staged binary constructions.
Record intermediate labels and copy choices rather than searching tensor powers
for a construction at runtime. Where the book gives a highest tuple, translate
it to the local factor order and verify its weight and vanishing raising actions.
A recipe requiring a later seed stays explicitly unavailable until that seed
backend is implemented.

Resolve dependencies recursively and cache each completed irrep. Distinguish
actual unit-label fundamentals from convenient doubled-spin bootstrap modules.
A fundamental cache miss must execute its own recipe rather than regenerate
all fundamentals and hope the requested one was saved. Keep recipes queryable
without enumeration; the exceptional dimension tables are not an instruction
to allocate every fundamental eagerly.

## Generalized tensor-component construction

The intended boundary is a function accepting two known `Irep`s, a target
q-label, and a deterministic copy-selection convention, returning a new `Irep`.
Its final name can be chosen during implementation. Generalize the numerical
core of `get_irep_frommw`; retain its current call sites as the maximal-component
special case so existing exterior-product and recursive constructions still work.

### 1. Build the ambient product

Reuse `tensor_prod_reps` and its sparse lowering blocks, weight ranges, and
`WeightInfo`. Preserve ordered input factors. Compute the target weight using
`qlab2mwz`; do not choose it from the largest ambient weight.

### 2. Select a numerical highest-weight vector

At the target weight, find a nonzero vector annihilated by all ambient raising
operators. Obtain raising maps from the known factor representations and their
exact metrics/normalization. The common kernel isolates highest vectors without
requiring numerical representations of the larger product components.

Reuse the CG3 construction's exact metric operations, deterministic candidate
ordering, and normalization where applicable. Its `get_mwstate` currently
projects away previously constructed multiplets; it cannot be called unchanged
without those subspaces. Either adapt that machinery to the target raising
kernel, or provide the required preceding subspaces without assuming the unknown
output irrep. Prefer direct target-kernel extraction to unnecessary construction
of unrelated components.

If the target occurs more than once, select a copy by a fixed exact kernel-basis
ordering and normalization (including sign). Reject an absent target or an
out-of-range copy. A crystal highest tuple identifies a combinatorial component;
it is not generally the numerical highest-weight pure tensor. Use its crystal
component to schedule lowering and a selected numerical highest vector to embed
an isomorphic copy; do not claim a canonical crystal-to-kernel copy matching.

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
to express ambient lowering in the generated basis. Construct inverse metrics
and return `Irep`. Preserve primitive rescaling, sign conventions, integrality,
and overflow handling; let existing numerical failures propagate. Save the
result only after successful construction.

The generated basis columns are the embedding data corresponding to a CG
multiplet. Returning a full CG tensor is not required for fundamental bootstrap.
Once the new irrep exists, ordinary CG3 construction may use it normally.

## Integration sequence

1. Implement and integrate direct minuscule seeds, including SO(7) spin.
2. Generalize the selected-component builder and extract SO(7) `(0, 1, 0)`
   from two defining irreps without preloading the target.
3. Encode the classical and exceptional fundamental recipes in fixed dependency
   order, including the staged binary forms of the table's tensor powers.
4. Replace the blanket exterior-power bootstrap with per-target recipe dispatch.
   Retain existing realizations as comparison cases or explicit seed recipes.
5. Keep general dominant-weight recursion through the maximal Cartan component
   of known fundamentals. Update label adapters, duality, and cache conventions.

## Completion

Use small cases in existing test files. Start with SO(7): verify the 7-dimensional
defining seed, 8-dimensional spin seed, and construction of the 21-dimensional
target from the vector square with no preexisting target cache entry. Check the
highest vector is killed by raising operators, the resulting dimensions and
weights, generator relations, and metric compatibility. Compare the target with
its existing exterior-square realization through an intertwining transformation.

Cover one repeated-copy target and an absent-target request. Check deterministic
copy selection, correct recursive recipe dispatch, and cache reload. Exercise
E6/E7 minuscule seeds and representative nonmaximal fundamental components;
F4/E8 recipes become executable when their later seed backends are ready.
Retain existing CG/recoupling checks, type-stability inspection, and relevant
allocation measurements. No new rank-check gate is required.

## Reference

Bump and Schilling, *Crystal Bases: Representations and Combinatorics* (2017),
user-supplied PDF at
`/home/lurlurlur/Downloads/Daniel Bump, Anne Schilling - Crystal Bases_ Representations and Combinatorics-World Scientific (2017).pdf`.

Relevant passages: §5.4 and §5.6, Tables 5.2–5.6.
