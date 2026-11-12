# Part 3a cleanup: symmetry-owned recipe and validation structs

Status: implemented. This document supplements [the minuscule plan](crystal_basis_minuscule_plan.md)
and supersedes its function-only recipe representation. It does not change
weight conventions, integer crystal words, Demazure traversal, or numerical
component extraction.

## Objective

Represent each supported symmetry's fundamental construction rules and
validation rules as concrete data structures. `getNsave_irep` must consume a
recipe object on a fundamental cache miss. Symmetry files define the data;
shared code interprets it without family branches.

Construction, crystal words, and validation share one inspectable source of
recipe data. The old sequence/source tables are removed; `minuscule_nodes`
derives its result from the catalogue.

Use a nonparametric abstract `FundamentalRecipe` with concrete subtypes for
each construction kind, and return a
`Vector{FundamentalRecipe}` for each symmetry. Likewise, return a
`Vector{IrrepValidationRule}` for validation. Do not add catalogue wrapper
types or specialize recipe types on rank, symmetry, or parameter schemas.
Family-specific catalogue-producing functions remain in the symmetry files.

## Data structures

The following documented types define the storage.

```julia
abstract type FundamentalRecipe end

struct MinusculeRecipe <: FundamentalRecipe
    target::Vector{Int}
    dimension::Int
end

struct AdjointRecipe <: FundamentalRecipe
    target::Vector{Int}
    dimension::Int
end

struct ProductRecipe <: FundamentalRecipe
    target::Vector{Int}
    sources::NTuple{2,Vector{Int}}
    copy_index::Int
    dimension::Int
end

struct IrrepValidationRule
    target::Vector{Int}
    sources::NTuple{2,Vector{Int}}
end
```

The revised [Part 3b plan](crystal_basis_virtual_plan.md) retains
`VirtualRecipe <: FundamentalRecipe` with target, dimension, ambient symmetry,
and operator words. Those words construct the target crystal only. The numerical interpreter uses target-string coefficients and metrics,
including F4's prescribed zero-weight block, and sparse commutator validation.
The ambient numerical interpreter has been removed; the recipe type and
catalogue container are unchanged.

`fundamental_recipes(S)` returns a `Vector{FundamentalRecipe}`, indexed by
fundamental node. Dependencies are requested recursively; there is no order
field. There are no rank/symmetry/payload type parameters or catalogue wrappers.
The abstract vector element type is intentional; each concrete subtype's fields
remain concretely typed. No kind tag or unused fields are needed.

Example SO(8) entries:

```julia
FundamentalRecipe[
    MinusculeRecipe([1,0,0,0], 8),
    ProductRecipe([0,1,0,0], ([1,0,0,0], [1,0,0,0]), 1, 28),
    MinusculeRecipe([0,0,1,0], 8),
    MinusculeRecipe([0,0,0,1], 8),
]
```

A product's sources belong to the target symmetry and retain their factor
order. The concrete recipe type selects the construction, not the mathematical
name of the result: SO(8)'s adjoint above uses `ProductRecipe`.
The `AdjointRecipe` numerical backend remains deferred to Part 3c.

Keep recipe labels as vectors to avoid creating a different metadata type for
every rank. At the existing numerical API boundary, check label length and
convert to `NTuple{nzops(S),Int}`. Rank specialization belongs to the existing
numerical representation code, not to recipe storage. Do not convert labels
inside per-weight or per-matrix-element loops.

`dimension` is the expected dimension, calculated with checked integer
arithmetic for rank-dependent families. Overflow is an error. Recipes contain
no `Irep`, crystal vertices, numerical blocks, closures, database handles, or
numeric coefficient type `RT`.

A validation rule independently extracts its target from one ordered source
pair and compares it with the target's canonical construction. Multiple rules
with the same target express different construction paths. Every occurrence of the target in that product must be validated. There is
no copy-index field in a validation rule. These product checks supplement
validation of every incidental copy generated during fundamental construction.

### Storage and specialization

The subtype represents a construction invariant: minuscule and adjoint seed
recipes have no sources, while product recipes have exactly two ordered sources
and a selected copy. Subtypes remove unused fields and kind-tag checks. This
small hierarchy is explicitly requested; it is not a general graph framework.

Indexing `Vector{FundamentalRecipe}` produces an abstractly typed value. Accept
dynamic dispatch once at the recipe execution boundary, outside numerical
inner loops. Use methods such as `construct_recipe(S, RT, recipe::ProductRecipe)`
and corresponding seed methods. Inside each method, concrete fields and the
existing rank-specific numerical routines retain their normal specialization.

Do not add type parameters or heterogeneous tuples to eliminate this setup
cost. Inspect inference inside the concrete execution methods, rather than
requiring a catalogue lookup to return a statically known subtype.

The label/source vectors are mutable even though their containing struct is
immutable. Return fresh owned metadata initially, and treat it as read-only
throughout construction. Copy labels if they must be modified; do not mutate
recipe labels during canonical splitting. No global catalogue cache is needed for this cleanup.

## Catalogue constructors and invariants

Provide two entry points:

```julia
fundamental_recipes(S)       # Vector{FundamentalRecipe}
irrep_validation_rules(S)    # Vector{IrrepValidationRule}
```

These functions construct metadata only. They must not enumerate Weyl orbits,
build crystal words, read databases, or construct numerical representations.

Validate the catalogue before interpretation:

1. Exactly one entry exists per fundamental node, with rank-length labels.
2. Entry k has the unit target label at k and a positive dimension.
3. Every entry uses a supported recipe subtype; product copy indices are
   positive. Source count is enforced by the product field type.
4. Product sources have the target symmetry's rank and nonnegative labels.
   Repeated sources are allowed. The dependency graph must be acyclic.
5. Validation labels have the target symmetry's rank and nonnegative components.
   No multiplicity copy can be excluded by a validation-rule field.

The hard-coded dependencies define a well-founded construction sequence.
Establish that invariant when writing and reviewing the catalogues; do not add
runtime cycle detection, an active-request registry, or a separately stored
execution order. Product-target occurrence and multiplicity still require the
existing crystal/numerical checks.

B/G2 remain unavailable until virtual-crystal work. Do not restore a legacy
fallback. C1, D2, and D3 retain their explicit existing low-rank handling;
SO(2) is not a simple D-family catalogue.

## Concrete construction catalogues

Write `F_k` for the unit Dynkin label at node k. Each arrow below becomes one
product entry with copy index 1. Every seed becomes a minuscule entry.

| Family | Seed entries | Product entries |
| --- | --- | --- |
| SU(r+1) | Every `F_k` | None |
| Sp(2r) | `F_1` | `F_k ← F_{k-1} tensor F_1`, k=2,…,r |
| SO(2r), r≥4 | `F_1`, `F_{r-1}`, `F_r` | `F_k ← F_{k-1} tensor F_1`, k=2,…,r−2 |
| E6 | `F_1`, `F_6` | `F_2 ← F_1 tensor F_6`; `F_3 ← F_1 tensor F_1`; `F_5 ← F_6 tensor F_6`; `F_4 ← F_3 tensor F_1` |
| E7 | `F_7` | `F_1 ← F_7 tensor F_7`; `F_6 ← F_7 tensor F_7`; `F_2 ← F_1 tensor F_7`; `F_5 ← F_6 tensor F_7`; `F_3 ← F_1 tensor F_1`; `F_4 ← F_5 tensor F_7` |

Use the dimension formulas and exceptional dimension tables in the minuscule
plan to populate each entry. D2 has minuscule entries at nodes 1 and 2;
D3 has them at nodes 1, 2, and 3; C1 has only node 1.

For SO(8), the complete data is:

| Node | Target | Kind | Ordered sources | Copy | Dimension |
| --- | --- | --- | --- | ---: | ---: |
| 1 | `(1,0,0,0)` | Minuscule | — | — | 8 |
| 2 | `(0,1,0,0)` | Product | `(1,0,0,0)`, `(1,0,0,0)` | 1 | 28 |
| 3 | `(0,0,1,0)` | Minuscule | — | — | 8 |
| 4 | `(0,0,0,1)` | Minuscule | — | — | 8 |

Requesting node 2 recursively obtains node 1 and then extracts the adjoint.
Nodes 3 and 4 are not constructed for that request. Reading the catalogue
itself constructs no numerical representations.

## `getNsave_irep` execution

Keep a cache hit as a fast return without building a catalogue. For a cache
miss:

1. Validate the symmetry and label. Initialize and save the trivial irrep when
   necessary, using `initialize_symmetry`.
2. If the target is fundamental, obtain its entry from `fundamental_recipes(S)`.
3. Recursively request each required source through `getNsave_irep`, preserving
   its position in the ordered source pair. Load completed dependencies from
   the normal cache/database; do not build unrelated fundamentals.
4. Dispatch on the concrete recipe subtype. Minuscule entries use the
   orbit constructor; virtual and adjoint entries use their dedicated backends.
   Product entries obtain their ordered source irreps, construct the product,
   and call `get_irep_frommw` with the target and construction copy index.
5. Check the result's target label and dimension. Retain exact integrality,
   uniform highest-weight normalization, and primitive inverse-metric storage.
6. Save and return the successfully constructed target. Never save incidental
   components under canonical cache keys.

Nonfundamental requests continue through the deterministic canonical split
and maximal-component construction. Their fundamental dependencies use the
same catalogue interpreter. No separate seed/product switch belongs in
`getNsave_irep` beyond invoking that interpreter.

The interpreter must accept optional incidental-result collection for
validation without changing production cache behavior. Metadata can be passed
through shared helpers to avoid repeatedly allocating catalogues during one
construction operation; keep it out of the numerical inner loops.

## Crystal construction uses the same catalogue

`mw_crystal_qlabel`, `_build_fundamental_word`, and minuscule alphabet setup
must obtain seed membership and ordered fundamental dependencies from the
catalogue. Preserve the existing integer symbol assignments and fixed
Demazure order. Symbol encoding, root coordinates, and duality remain methods
in the symmetry files; they are not recipe fields.

Catalogue construction must not call these crystal functions. This direction
of dependency avoids a catalogue → crystal → catalogue initialization cycle.

## Validation-rule catalogues

Migrate the existing finite product rules into `Vector{IrrepValidationRule}`.
Preserve their order, targets, and source factor order. Expand validation
from the former copy-1 checks to all occurrences of each target. A shared
helper may fill the adjoint-in-fundamental-times-dual entries, but the resulting
catalogue must contain explicit labels; no numerical decomposition is used to
invent rules at runtime.

For SO(8), retain the adjoint checks and these entries in particular:

| Target | Ordered sources |
| --- | --- |
| `(0,0,2,0)` | `(0,0,1,0)`, `(0,0,1,0)` |
| `(0,0,2,0)` | `(0,1,0,0)`, `(0,1,0,0)` |
| `(0,0,0,2)` | `(0,0,0,1)`, `(0,0,0,1)` |
| `(0,0,0,2)` | `(0,1,0,0)`, `(0,1,0,0)` |

Each extraction runs independently of the target cache; the canonical
reference is requested separately through `getNsave_irep`. Extend the component
extraction interface to visit every target highest-weight state, rather than
returning after its first copy. Establish completion from exhausted target
multiplicity; do not treat an arbitrary construction exception as end of
copies. An absent target is a failed validation rule, not a vacuous success.

Keep two complementary validation procedures:

1. **Incidental components:** request the fundamentals selected for validation
   through the recursive interpreter, recording every generated nonfundamental
   irrep and every copy. Enumeration of the selected targets is not a separate
   dependency schedule. After all required fundamentals are saved, obtain
   canonical references and compare each
   candidate. Reference construction must not recursively launch validation.
2. **Explicit product rules:** consume the vector of validation-rule entries and execute
   every target copy for each selected rule against its canonical reference. Allow an explicit
   subset for large exceptional runs, reporting what was actually checked.

Both procedures require equal `Irep` types and equal `Sl`, `Sz`, `innerprod`,
`inv_innerprod`, `qlabel`, and `dimension`. Exclude only `size_byte`. Compare
stored inverse matrices and factors after making matrix entries relatively
prime. Only inverse metrics may be rational; reject nonintegral results for
integer representation storage. An intertwiner is not a replacement for this
field-equality requirement.

Existing dimension tests, including Sp(6), and external large-scale tests
remain responsible for their existing mathematical checks. Do not duplicate
those suites in the catalogue implementation.

## Cleanup and file placement

Define the abstract recipe type, currently required concrete subtypes, and
validation entry before including the symmetry files, next to the common
symmetry declarations in `Symmetry/Symmetry.jl`. Their fields contain only
recipe metadata; they must not depend on `Irep` or `CrystalVertex`, which are
defined later. Add the virtual subtype with the next step.

Keep catalogue-producing methods and all family-specific data in the existing
symmetry files, with shared metadata checks beside the recipe types in
`Symmetry/Symmetry.jl`. Keep recipe execution and component extraction in
`irep.jl`, crystal traversal and word caches in `Symmetry/crystal.jl`, and
cross-construction checks in `irep_validation.jl`. Reserve `minuscule.jl` for
minuscule orbits, seeds, and alphabets. Preserve existing metadata access and
numerical APIs.

Replace or retire the following duplicated sources of rule data:

- `fundamental_sequence`: remove it; recursive dependencies replace it.
- `minuscule_nodes`: derive seed nodes from catalogue entries.
- `fundamental_sources`: read the selected product entry's ordered sources.
- `fundamental_validation_rules`: replace tuple-only output with the validation
  catalogue entry point and update consumers.
- `fundamental_qlabels`: derive labels from the construction catalogue where it
  supplies the same existing ordering; do not silently change an established
  caller-visible order.

Small compatibility accessors may delegate to the catalogue if existing callers
need them. They must not maintain independent tables. An empty
`fundamental_recipes(S)` catalogue denotes an unavailable numerical backend;
do not maintain a separate symmetry-support predicate. Do not reintroduce `get_fundamental_ireps`, alternating products,
or B/G2 legacy seeds.

## Recursive construction and validation details

### Recursive dependency construction

The catalogue author supplies acyclic dependencies. Recursive requests follow
those dependencies and save completed results. For `V tensor V`, construct V
once and reuse the cached result for the second factor. Shared dependencies in
different branches work the same way. No runtime cycle machinery is required.

### Dependency execution versus tensor-factor order

Recursive availability does not change the mathematical source pair. If a
recipe says `X tensor Y`, pass X as the first factor and Y as the second even
if Y was already cached or constructed earlier. Sorting labels or swapping
sources would change the ambient basis and copy-selection convention.

### Construct one copy, validate every copy

A construction recipe selects one deterministic highest-weight copy to return
one `Irep`. A validation rule does not make that selection: enumerate every
occurrence of its target in the product, using the existing exact projection
against preceding components and copies. Compare each resulting `Irep` with
the canonical reference, excluding only `size_byte`.

Complete validation only when the target's remaining multiplicity is exhausted.
An absent target fails the rule. An arithmetic or construction error also fails;
it must not be caught and interpreted as exhaustion. Report the number of
copies checked. This expands the former copy-1 validation behavior and needs
focused coverage using the existing repeated-copy test.

### Keep candidates separate from references

Save the required fundamentals first, retaining incidental nonfundamental
results outside the canonical cache. Then obtain canonical references through
normal recursive requests. A candidate must never populate the cache entry
that is later read as its own reference. Canonical-reference construction must
not recursively initiate another validation run.

Virtual-crystal details are deferred to the next step. The metadata cleanup
must neither implement them nor claim that the planned virtual subtype enables them.

## Implementation and acceptance sequence

1. Add the documented abstract recipe type, concrete subtypes, validation
   entry type, and simple catalogue invariant checks.
2. Populate construction and validation catalogues in each supported symmetry
   file, including low-rank cases and explicit dimensions.
3. Route fundamental cache misses through one shared interpreter.
4. Make crystal construction and validation consume the same catalogue data.
5. Remove redundant dispatch tables and update documentation and callers.
6. Review include order, concrete subtype fields, execution-method inference,
   cache behavior, and unchanged symbol/basis ordering.

Use existing test files. Add focused coverage for catalogue invariants and
on-demand consumption: request SO(8)'s adjoint from an empty database, verify
its vector dependency and result, and verify that unrelated half-spin seeds
were not saved. Check repeated/shared dependencies reuse completed results.
Extend existing repeated-copy coverage to verify that a single validation rule
checks every copy. Cover invalid labels/copy indices and unavailable adjoint execution without
pretending the deferred numerical backends exist.
Reuse existing field-equality, normalization, and dimension coverage.

This is intended to preserve numerical behavior. Follow the repository's
internal-refactor policy for unchanged numerical paths: run the DMRG runner
and check the Nkeep=20 and 50 energies against
`-82.37493703992092` and `-82.39703156828021`. Tests of newly exposed catalogue
validation/error semantics are a separate behavior change and should be
limited to those contracts. Do not add a broad new integration suite or run a
Telum documentation build.

Completion means family rules exist as data, every fundamental request follows
that data, and both validation procedures consume the vectors of rule structs.
Moving methods between files alone is not completion. A metadata-only cleanup
must not change the numerical cache namespace; if validation detects a basis
change, investigate it rather than masking it with cache invalidation.
