# Part 3: Numerical representations, symmetry integration, and caches

Status: design proposal, split from the source review at LurCGT commit `d1f4176`.
No representation code or cache data is changed by these documents.

This is one of three implementation plans:

1. [Unified crystal representation](crystal_basis_irrep_plan.md).
2. [Introduce reduced words and string patterns](crystal_basis_string_plan.md).
3. [Numerical representations, symmetry integration, and caches](crystal_basis_integration_plan.md).

## Scope and dependencies

Complete the work beyond crystal construction and string coordinates: exact
numerical seed realizations, ambient highest-weight extraction, public symmetry
coverage, duality/conjugation, and consistent cache migration. Retain `Irep`
and its sparse weight-block storage.

Part 1 replaces tableau vertices with tensor-word vertices and fundamental
embeddings within the existing traversal. The extended seed catalogue below remains staged integration work.
Part 2 supplies ordered string labels and proposed numerical schedules. Numerical realization
and API design can proceed alongside those parts, but production activation
requires their relevant acceptance gates. Introduce the cache namespace before
writing any changed-basis representation or dependent tensor data.

## Current construction and its limitations

| Source | Current behavior | Proposed change |
| --- | --- | --- |
| [`src/Symmetry/crystal.jl`](../src/Symmetry/crystal.jl) | `Tableau`, `Crystal_ops`, letter transitions, and a tensor-signature lowering rule encode crystals through vector alphabets. `det_getvector` asks whether an entire lowering string contains unseen tableaux. | Replace the tableau vertex type with tensor-word vertices while retaining the traversal. Fundamental vertices have tensor embeddings; string schedules belong to part 2. |
| [`src/clebsch.jl`](../src/clebsch.jl), `get_irep_frommw` | Traverses tableaux and numerical vectors together; the first accepted route determines the vector. Columns within a weight are subsequently assigned by iterating `tabsNvecs`, a dictionary. | Part 1 retains this traversal with the replacement vertex type. Part 2 separately introduces string ordering and specified lowering schedules. |
| Same file, `get_fundamental_ireps` | Saves the trivial and defining modules, then extracts highest-weight pieces of exterior powers of the defining module. Even orthogonal types also extract the second middle-exterior-power piece. | Replace the assumption that this supplies all fundamentals with a per-node construction catalogue. Retain exterior powers as efficient numerical realizations where applicable. |
| Same file, `getNsave_irep`, `split_qlabel`, `is_fundamental` | Builds other modules recursively from smaller highest weights. Orthogonal recursion includes doubled spin labels and parity restrictions. | Distinguish actual fundamental weights from convenient vector-derived bootstrap modules; include genuine spin seeds before enabling all dominant labels. |
| [`src/irep.jl`](../src/irep.jl) | `Irep` stores integer sparse lowering blocks, weight ranges, exact Gram matrices, and inverse-metric factors. `def_altprod` supplies numerical exterior-power realizations. | Keep these structures and operations; audit integrality whenever the basis changes. |
| [`src/Symmetry/Symmetry.jl`](../src/Symmetry/Symmetry.jl) | Supported non-Abelian families are SU, SO, Sp, G2. `get_dualq` reverses SU labels and otherwise returns the same label. | Add exceptional root data in a later stage; derive duality from the opposition involution of the longest Weyl element. |
| [`src/sqlite_io.jl`](../src/sqlite_io.jl) | Irrep memory keys use `(S, RT, qlabel)`; the basis prescription is not represented there. | Version the basis convention consistently across representations and all dependent CG/recoupling caches. |

These are source-level findings, not newly reproduced numerical failures.
The existing method is already crystal-guided; the proposed improvement is a
specified mathematical path and broader seed coverage, not the introduction of
crystals for the first time.

Specific gaps relevant to this work:

- `SO.jl:preprocess` rejects odd terminal spin labels in type B and odd sums of
  terminal labels in type D. Vector tensor powers cannot supply spin modules:
  their weights lie in a different lattice coset. Doubling z-coordinates already
  accommodates half-integral weights, but does not construct their modules.
- `is_fundamental` recognizes unit labels even though the orthogonal exterior
  bootstrap does not create spin fundamentals. A cache miss can therefore reach
  a bootstrap/reload path that does not produce the requested representation.
- `fundamental_qlabels` exists for SU, Sp, G2, but not SO.
- `get_irep_frommw` starts with `RT[1]` and assumes a one-dimensional ambient
  highest-weight space. Extracting an arbitrary fundamental constituent of a
  tensor product requires an actual highest-weight vector in a potentially
  larger ambient weight space, not merely a requested `mw`.
- `is_valid_tableau` defaults to true. A connected component of a proven normal
  tensor crystal is a clearer correctness boundary than unspecified tableau
  admissibility checks. This observation alone does not show current components
  are incorrect.
- The current duality default needs extension for type D of odd rank (exchange
  the two spin nodes) and E6 (exchange 1↔6, 3↔5). Review `getNsave_1jsym`, whose
  source already flags spin support, and the associated conjugation conventions.

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
they do not generally descend to ``\mathrm{SO}(N)``. Decide and document whether
`SO{N}` denotes the Lie algebra/simply connected representation category or the
SO group category before widening its accepted labels. This proposal does not
silently change that public contract. Treat B1=A1, C1=A1, D3=A3, D2=A1⊕A1
separately; SO(2) is Abelian rather than a member of the simple B/D catalogue.

## Additional seed constructions

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

For the E8 seed use the adjoint crystal: one vertex per root and **eight distinct
zero-weight vertices**, with the arrows of §5.5. For numerical matrices use the
adjoint action in a Chevalley basis, including its structure constants and metric;
do not collapse zero-weight multiplicity or substitute adjacency matrices.

For virtualization record ambient type, node orbits ``O_i``, scale factors
``\gamma_i``, and

```math
\Psi(\omega_i)=\gamma_i\sum_{j\in O_i}\widehat\omega_j,
\qquad
\widetilde f_i^{\rm virt}=\prod_{j\in O_i}
\widehat{\widetilde f}_j^{\gamma_i}.
```

Use the analogous raising operators. Orbit nodes commute. At every retained
vertex check equality of ambient string lengths within each orbit, divisibility
by ``\gamma_i``, closure under both virtual operators, and the weight map.
These are the alignment and normality requirements in §§5.1–5.3, not optional
checks on a visually folded graph.

The book's explicit exceptional constructions (§5.6, printed p. 70) are:

| Target seed | Ambient crystal and highest vertex | Virtual lowering operators in target node order |
| --- | --- | --- |
| F4, ``B(\omega_4)`` (26) | E6 ``B(\widehat\omega_1)\otimes B(\widehat\omega_6)``, highest⊗highest | ``\widehat f_2^2,\ \widehat f_4^2,\ \widehat f_3\widehat f_5,\ \widehat f_1\widehat f_6`` |
| G2, ``B(\omega_1)`` (7) | D4 ``B(\widehat\omega_1)\otimes B(\widehat\omega_3)\otimes B(\widehat\omega_4)``, all highest | ``\widehat f_1\widehat f_3\widehat f_4,\ \widehat f_2^3`` |

Hats in this table mean **crystal** operators. In F4 the factors are
``(2,2,1,1)``; in G2 they are ``(1,3)``. These are the book's similarity
embeddings, not naive restriction/branching of the 27 or 8 dimensional module.
In particular, a product of classical ambient lowering matrices is generally
not a target Chevalley generator: its weight change and commutators differ.

The separate F4 numerical realization and its acceptance gate are specified in
[numerical seed realizations](#numerical-seed-realizations).
Keep the existing validated G2 numerical seed.

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
| F4 | 4→26 | Virtual E6 seed V described above |
| F4 | 1→52; 3→273 | Respective components of ``V^{\otimes2}`` |
| F4 | 2→1,274 | ``C_{\omega_2}(V^{\otimes3})`` |
| G2 | 1→7 | Existing numerical seed; virtual D4 crystal as a uniform construction/check |
| G2 | 2→14 | ``C_{\omega_2}(B(\omega_1)^{\otimes2})``; numerically retain extraction from the alternating square |

E8 has no minuscule seed. Set A=``B(\omega_8)`` (the 248-dimensional adjoint).
The book gives explicit highest tuples, equivalently the following incremental
recipes. Here ``v_\beta`` is the adjoint vertex of root weight β, and ``u_j``
is the selected highest vertex of ``B(\omega_j)``. Use the book's tensor order.

| Target node | Highest tuple generating its component | Dimension |
| --- | --- | ---: |
| 8 | Adjoint highest root vertex ``v_{\omega_8}`` | 248 |
| 1 | ``v_{\omega_1-\omega_8}\otimes u_8`` | 3,875 |
| 7 | ``v_{\omega_7-\omega_8}\otimes u_8`` | 30,380 |
| 2 | ``v_{\omega_2-\omega_1}\otimes u_1`` | 147,250 |
| 6 | ``v_{\omega_6-\omega_7}\otimes u_7`` | 2,450,240 |
| 3 | ``v_{\omega_3-\omega_2}\otimes u_2`` | 6,696,000 |
| 5 | ``v_{\omega_5-\omega_6}\otimes u_6`` | 146,325,270 |
| 4 | ``v_{\omega_4-\omega_5}\otimes u_5`` | 6,899,079,264 |

Each displayed difference is a root, with the string lengths listed in Table
5.4. Verify that all raising operators kill each tuple before enumeration.
This records a construction for every node without claiming eager storage of
all E8 fundamentals is practical. Node 4 alone has billions of vertices.

### Turning the catalogue into executable recipes

Use ordinary functions and concrete tuples/arrays initially; no new hierarchy,
wrapper, or registry struct is needed. Store target label, source factors,
construction kind, fold parameters if applicable, expected dimension, and a
highest-vertex certificate. For recipes without an explicit highest tuple,
search ordered factor vertices subject to the target weight and signature
conditions, and select the lexicographically first highest tuple. Persist the
resulting certificate, rather than repeating a whole tensor search per request.
Multiplicity of the target component requires an explicit copy-selection rule.

Numerical highest-weight extraction is specified in
[numerical highest-weight extraction](#numerical-highest-weight-extraction).

For general λ, tensor the fundamental crystals in a fixed node order and take
the component of the tensor of their highest vertices. Numerically the analogous
Cartan component has highest weight equal to the sum and multiplicity one. This
preserves the intent of existing recursive construction once all seeds exist.

Enumerate only the requested component. Use staged products and precomputed
weight ranges; avoid allocating the entire tensor Cartesian product or dense
Kronecker operators. A complete recipe catalogue and a completely materialized
catalogue are different deliverables. Add dimension/memory estimates before
large exceptional constructions and keep recipes queryable without enumeration.

## Numerical seed realizations

Follow the numerical column of the
[fundamental catalogue](#complete-fundamental-construction-catalogue).
Retain exterior-power realizations where applicable. The unit-coefficient
Chevalley shortcut is specific to minuscule modules; general crystal adjacency
does not determine lowering matrices or metrics. E8 requires adjoint action in
a Chevalley basis, including structure constants and the zero-weight metric.

For F4 numerical matrices, a separate explicit realization is therefore required.
A general exact highest-weight-module construction is a suitable fallback:
use Chevalley–Serre relations, a PBW ordering and the highest-weight relations
``E_i v_\lambda=0``, ``H_i v_\lambda=\lambda_i v_\lambda``,
``F_i^{\lambda_i+1}v_\lambda=0``; compute the contravariant form and quotient
its radical weight by weight. This is substantial work and must be validated
before claiming virtual crystals alone supply a new `Irep` backend. A specialized
26-dimensional F4 realization can replace this fallback if its normalization and
relations are proved. Keep the existing validated G2 numerical seed.

## Numerical highest-weight extraction

Numerically, a crystal highest tuple labels a component but is not generally a
classical highest-weight pure tensor. Solve the common kernel of the raising
maps in its ambient weight space and choose a deterministic normalized vector
(or basis of the kernel when there are multiple copies). Extend extraction to
accept that vector. Existing `get_mwstate`/multiplet machinery is a useful starting
point; `get_irep_frommw`'s scalar initial vector is insufficient here.

For general highest weights, preserve the recursive Cartan-component approach
once the required fundamental seeds exist. Use fixed factor/node order, staged
products, and precomputed weight ranges. Avoid full tensor Cartesian products
and dense Kronecker operators. For multiple copies, make the relationship
between the chosen crystal certificate and numerical kernel basis explicit.

## Exact arithmetic and storage

Apply [part 2's numerical basis gate](crystal_basis_string_plan.md#a-crystal-path-is-not-a-numerical-representation-matrix),
including its monomial-basis justification, exact rank/metric checks, deterministic
rescaling, and explicit integrality/overflow audit for `Irep.Sl`. A crystal graph
or reconstruction theorem alone does not justify storing new numerical matrices.

## Symmetry contracts and integration

1. Distinguish actual fundamental weights from vector-derived bootstrap labels;
   provide complete `fundamental_qlabels` and make splitting/bootstrap logic
   consistent with available seeds.
2. Resolve whether `SO{N}` denotes the Lie algebra/simply connected category or
   the SO group category before accepting spin labels. Preserve external q-label
   meaning through the node adapters and handle the low-rank exceptions in
   the conventions below. Do not silently broaden the existing public contract.
3. Derive duality from the opposition involution of the longest Weyl element.
   Include exchange of D odd-rank spin nodes and E6 nodes 1↔6, 3↔5 in the
   catalogue convention. Update `getNsave_1jsym` and conjugation together with
   spin support.
4. Add exceptional root-data/backend integration in stages. Keep the public
   symmetry API unchanged initially; new symmetry types or public exports need
   a separate design justification and approval.

## Cache migration and production activation

Introduce a consistent basis-version namespace for both memory and SQLite data.
Include the chosen reduced word, within-weight ordering, and normalization
prescription in the versioned convention. Basis changes affect CGT, 1j, F/R/X
symbols and derived decompositions; a change to irrep cache keys alone is
insufficient. Retain old caches separately or perform explicit basis
transformations. Never mix gauges.

Implement the namespace before activating changed basis construction. Enable
validated families/recipes incrementally; neither graph dimension checks nor
availability of a recipe certifies numerical backend support. Keep large E8
recipes queryable without eager materialization, and require dimension/memory
estimates before expensive requests.

## Validation and completion gates

- Check Weyl dimensions and weight multiplicities, exact rank and nondegenerate
  metric in every constructed block, Chevalley commutators, and Serre relations
  with the documented Cartan convention. Verify the metric-adjoint relation
  after mapping existing generator normalizations to Chevalley normalization.
- Check spin dimensions, distinct D chiralities, coordinate round trips, duality
  of odd-rank D spinors and E6's 27 pair, and the singlet in a module times its
  dual. Compare the G2 numerical representation with its existing seed.
- Verify that changed discovery orders produce identical ordered data and that
  cache reload preserves the basis version. Compare old/new numerical matrices
  by intertwining transformations, not assumed entrywise equality.
- Run existing CG orthogonality and recoupling tests; preserve stored w-matrix
  column isometry, CGT leg ordering, and downstream Telum invariants.

Use `@code_warntype` or JET at payload function barriers and measure allocation
and runtime impact with existing benchmarks where applicable. Keep per-weight
payload writes independent and avoid shared mutable accumulators in heavy fills.

Completion requires validated numerical representations for each enabled recipe,
a documented public label contract, coherent duality/conjugation support, and
versioned dependent caches. Catalogue completeness alone does not establish
these conditions.

## Shared implementation policy

Part 1 proposes only the requested replacement vertex type and documents its
justification; it preserves the existing traversal and operator-data structure. Use ordinary functions and concrete tuples/arrays around
them; additional types require their own justification and approval. Every new named Julia function
needs a docstring. Keep containers concrete, preallocate output slots, and make
heavy per-vertex/per-weight work independent where possible.

For behavior-changing implementation, add small cases to existing relevant
files such as `test/runtests.jl` and `test/g2_symmetry_smoketest.jl`; justify a new
file only if no existing file fits. Reproduce any bug before fixing it. For pure
internal refactors follow the prescribed DMRG-only runner; new behavior follows
the normal relevant-test policy.

## References and validation provenance

See [part 1](crystal_basis_irrep_plan.md#references-and-provenance) for the
chapter 2 foundation and part 2 for the string coordinate justification. Numerical realization gates, storage choices, API
integration, and cache migration are engineering recommendations from the
original source review.

The seed catalogue retains the original review's provenance. Part 1 separately
records the reference passages consulted for its rewrite. No numerical tests,
benchmarks, or documentation builds were run for these documentation edits.

### Seed catalogue provenance

Primary reference: Daniel Bump and Anne Schilling, *Crystal Bases:
Representations and Combinatorics*, World Scientific (2017), the user-supplied
PDF at `/home/lurlurlur/Downloads/Daniel Bump, Anne Schilling - Crystal Bases_ Representations and Combinatorics-World Scientific (2017).pdf`.

- §§5.1–5.3: weight embeddings, virtual operators and alignment.
- §5.4, Proposition 5.11 and Proposition 5.14: classical fundamentals, missing
  spin modules, minuscule construction and optional B-spin virtualization.
- §5.5: adjoint crystal, including separate zero-weight vertices.
- §5.6, Tables 5.2–5.6 and Theorems 5.17–5.19: exceptional recipes. Printed
  pages 69–70 (PDF pages 81–82) were visually checked in the original review
  because text extraction loses superscripts in the folding formulas.
- §11.1, printed pp. 144 onward: string definition and operator order;
  Corollary 13.8: good-word result for normal highest-weight crystals.

Supplementary implementation references: [Sage crystal component selection](https://doc.sagemath.org/html/en/thematic_tutorials/lie/crystals.html),
[Sage minuscule representations](https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/quantum_groups/representations.html),
and [Geck, *Minuscule weights and Chevalley groups*](https://arxiv.org/abs/1603.07179).
The sequencing, storage choices, numerical realization gates and cache migration
above are engineering recommendations from this review, not claims that the book
provides a ready-made classical matrix algorithm.

