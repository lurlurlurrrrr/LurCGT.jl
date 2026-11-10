# Part 3a cleanup: symmetry-owned recipe and validation structs

Status: implementation proposal. Part 3a is not complete until this cleanup is
implemented. This document supplements [the minuscule plan](crystal_basis_minuscule_plan.md)
and supersedes its function-only recipe representation. It does not change
weight conventions, integer crystal words, Demazure traversal, or numerical
component extraction.

## Objective

Represent each supported symmetry's fundamental construction rules and
validation rules as concrete data structures. `getNsave_irep` must consume a
recipe object on a fundamental cache miss. Symmetry files define the data;
shared code interprets it without family branches.

Currently, `fundamental_sequence`, `minuscule_nodes`, `fundamental_sources`,
and `fundamental_validation_rules` distribute the rules across executable
methods. Merely moving those methods to symmetry files does not finish the
cleanup: construction, crystal words, and validation need one inspectable
source of recipe data.

Use shared struct definitions with a separate concretely typed instance for
each symmetry, such as `FundamentalRecipes{SO{8},4}`. Separate named types such
as `SURecipes` and `SORecipes` are unnecessary because their data layout and
interpreter are the same. All family-specific catalogue construction stays in
`Symmetry/SU.jl`, `Sp.jl`, `SO.jl`, `E6.jl`, and `E7.jl`.

## Proposed data structures

The following signatures specify the intended storage. Each new type and
constructor must receive a Julia-compatible docstring during implementation.

```julia
struct FundamentalRecipe{R}
    target::NTuple{R,Int}
    minuscule::Bool
    sources::NTuple{2,NTuple{R,Int}}
    copy_index::Int
    dimension::Int
end

struct FundamentalRecipes{S<:NonabelianSymm,R}
    entries::Vector{FundamentalRecipe{R}}
    order::Vector{Int}
end

struct IrrepValidationRule{R}
    target::NTuple{R,Int}
    sources::NTuple{2,NTuple{R,Int}}
    copy_index::Int
end

struct IrrepValidationRules{S<:NonabelianSymm,R}
    products::Vector{IrrepValidationRule{R}}
end
```

`entries[k]` is the recipe for fundamental node k. `order` is an explicit
permutation of those node indices in construction order. Keep entry indexing
and execution order separate: E6 and E7 do not construct in ascending node
order.

For a minuscule entry, `sources` contains two zero labels and `copy_index=0`;
these are unused fields, not instructions to multiply trivial irreps. Its
unit target label identifies the seed node. For a product entry, sources are
two ordered fundamental labels and `copy_index>=1`. Never sort the source
pair. Current recipes need no nonfundamental intermediate source entries.

`dimension` is the expected dimension, calculated with checked integer
arithmetic for rank-dependent families. Overflow is an error. Recipes contain
no `Irep`, crystal vertices, numerical blocks, closures, database handles, or
numeric coefficient type `RT`.

A validation rule independently extracts its target from one ordered source
pair and compares it with the target's canonical construction. Multiple rules
with the same target express different construction paths. Their copy indices
are explicit; these selected product checks do not replace validation of every
incidental multiplicity copy generated during fundamental construction.

### Why structs are appropriate here

- Function-only dispatch cannot expose a complete recipe as data without
  reconstructing information from several separate methods. The catalogue
  structs make dependencies, ordering, dimensions, and copy selection explicit.
- Entry structs enforce rank-length labels and homogeneous container elements.
  Family parameters on catalogues prevent accidental reuse across symmetries.
- Concrete vectors provide predictable storage and iteration. Shared numerical
  algorithms do not specialize on a different closure or recipe type per node.
- Check inference for catalogue lookup and the interpreter with
  `@code_warntype` or `@inferred`, including a runtime node index.
- This is a small data-driven replacement for the existing dispatch tables,
  not a new abstract hierarchy or general dependency-graph framework.

Struct fields containing vectors remain mutable through those vectors.
Initially return fresh, owned metadata from the catalogue constructors and
pass it through a construction operation without mutation. Do not introduce
an `Any`-valued global catalogue cache. Caching can be considered separately
if metadata allocation is measured to matter.

## Catalogue constructors and invariants

Provide two entry points:

```julia
fundamental_recipes(S)       # FundamentalRecipes{S,R}
irrep_validation_rules(S)    # IrrepValidationRules{S,R}
```

These functions construct metadata only. They must not enumerate Weyl orbits,
build crystal words, read databases, or construct numerical representations.

Validate the catalogue before interpretation:

1. `R == nzops(S)` and exactly one entry exists per fundamental node.
2. Entry k has the unit target label at k and a positive dimension.
3. `order` contains every node exactly once.
4. Seed entries have the specified unused source/copy fields.
5. Product sources are valid fundamental labels and occur earlier than their
   target in `order`; copy indices are positive. Repeated sources are allowed.
6. Validation labels have the correct rank and nonnegative components, and
   validation copy indices are positive.

These checks reject dependency cycles before numerical work. Product-target
occurrence and multiplicity still require the existing crystal/numerical
checks; a structurally valid catalogue alone cannot establish them.

B/G2 remain unavailable until virtual-crystal work. Do not restore a legacy
fallback. C1, D2, and D3 retain their explicit existing low-rank handling;
SO(2) is not a simple D-family catalogue.

## Concrete construction catalogues

Write `F_k` for the unit Dynkin label at node k. Each arrow below becomes one
product entry with copy index 1. Every seed becomes a minuscule entry.

| Family | Seed entries | Product entries | Stored order |
| --- | --- | --- | --- |
| SU(r+1) | Every `F_k` | None | `1:r` |
| Sp(2r) | `F_1` | `F_k ← F_{k-1} tensor F_1`, k=2,…,r | `1:r` |
| SO(2r), r≥4 | `F_1`, `F_{r-1}`, `F_r` | `F_k ← F_{k-1} tensor F_1`, k=2,…,r−2 | `[1,r-1,r,2,…,r-2]` |
| E6 | `F_1`, `F_6` | `F_2 ← F_1 tensor F_6`; `F_3 ← F_1 tensor F_1`; `F_5 ← F_6 tensor F_6`; `F_4 ← F_3 tensor F_1` | `[1,6,2,3,5,4]` |
| E7 | `F_7` | `F_1 ← F_7 tensor F_7`; `F_6 ← F_7 tensor F_7`; `F_2 ← F_1 tensor F_7`; `F_5 ← F_6 tensor F_7`; `F_3 ← F_1 tensor F_1`; `F_4 ← F_5 tensor F_7` | `[7,1,6,2,5,3,4]` |

Use the dimension formulas and exceptional dimension tables in the minuscule
plan to populate each entry. D2 has two minuscule entries in order `[1,2]`;
D3 has three in order `[1,2,3]`; C1 has one in order `[1]`.

For SO(8), the complete data is:

| Node | Target | Kind | Ordered sources | Copy | Dimension |
| --- | --- | --- | --- | ---: | ---: |
| 1 | `(1,0,0,0)` | Minuscule | unused | 0 | 8 |
| 2 | `(0,1,0,0)` | Product | `(1,0,0,0)`, `(1,0,0,0)` | 1 | 28 |
| 3 | `(0,0,1,0)` | Minuscule | unused | 0 | 8 |
| 4 | `(0,0,0,1)` | Minuscule | unused | 0 | 8 |

The catalogue stores `order = [1,3,4,2]`, without constructing any of these
representations when the catalogue is requested.

## `getNsave_irep` execution

Keep a cache hit as a fast return without building a catalogue. For a cache
miss:

1. Validate the symmetry and label. Initialize and save the trivial irrep when
   necessary, using `initialize_symmetry`.
2. If the target is fundamental, obtain its entry from `fundamental_recipes(S)`.
3. Find only the requested entry's dependency closure, and execute those entries
   in the catalogue's stored order. Load already available dependencies from
   the normal cache/database. Do not build unrelated fundamentals.
4. For a seed entry, call the existing minuscule orbit constructor. For a
   product entry, obtain its ordered source irreps, construct the product, and
   call `get_irep_frommw` with the target maximal weight and stored copy index.
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

Migrate the existing finite product rules into `IrrepValidationRules{S,R}`.
Preserve their order, targets, source factor order, and copy index 1. A shared
helper may fill the adjoint-in-fundamental-times-dual entries, but the resulting
catalogue must contain explicit labels; no numerical decomposition is used to
invent rules at runtime.

For SO(8), retain the adjoint checks and these entries in particular:

| Target | Ordered sources | Copy |
| --- | --- | ---: |
| `(0,0,2,0)` | `(0,0,1,0)`, `(0,0,1,0)` | 1 |
| `(0,0,2,0)` | `(0,1,0,0)`, `(0,1,0,0)` | 1 |
| `(0,0,0,2)` | `(0,0,0,1)`, `(0,0,0,1)` | 1 |
| `(0,0,0,2)` | `(0,1,0,0)`, `(0,1,0,0)` | 1 |

Each extraction runs independently of the target cache; the canonical
reference is requested separately through `getNsave_irep`.

Keep two complementary validation procedures:

1. **Incidental components:** construct and save fundamentals in catalogue
   order, recording every generated nonfundamental irrep and every copy.
   After fundamentals are saved, obtain canonical references and compare each
   candidate. Reference construction must not recursively launch validation.
2. **Explicit product rules:** consume the validation-rule struct and execute
   each selected rule against its canonical reference. Allow an explicit
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

Define the four shared structs before including the symmetry files, next to
the common symmetry declarations in `Symmetry/Symmetry.jl`. Their fields need
only symmetry types and integer labels, so they must not depend on `Irep` or
`CrystalVertex`, which are defined later.

Keep catalogue-producing methods and all family-specific data in the existing
symmetry files. Keep the generic interpreter and validation executors in
`minuscule.jl`. Preserve existing metadata access and numerical APIs.

Replace or retire the following duplicated sources of rule data:

- `fundamental_sequence`: read `catalogue.order`.
- `minuscule_nodes`: derive seed nodes from catalogue entries.
- `fundamental_sources`: read the selected product entry's ordered sources.
- `fundamental_validation_rules`: replace tuple-only output with the validation
  catalogue entry point and update consumers.
- `fundamental_qlabels`: derive labels from the construction catalogue where it
  supplies the same existing ordering; do not silently change an established
  caller-visible order.

Small compatibility accessors may delegate to the catalogue if existing callers
need them. They must not maintain independent tables. The unsupported-backend
check may stay a lightweight symmetry method so cache hits need not allocate
catalogues. Do not reintroduce `get_fundamental_ireps`, alternating products,
or B/G2 legacy seeds.

## Implementation and acceptance sequence

1. Add documented concrete struct definitions and catalogue invariant checks.
2. Populate construction and validation catalogues in each supported symmetry
   file, including low-rank cases and explicit dimensions.
3. Route fundamental cache misses through one shared interpreter.
4. Make crystal construction and validation consume the same catalogue data.
5. Remove redundant dispatch tables and update documentation and callers.
6. Review include order, concrete container types, runtime-index inference,
   cache behavior, and unchanged symbol/basis ordering.

Use existing test files. Add focused coverage for catalogue invariants and
on-demand consumption: request SO(8)'s adjoint from an empty database, verify
its vector dependency and result, and verify that unrelated half-spin seeds
were not saved. Check invalid dependency order is rejected before numerical
construction and that validation consumes stored copy indices. Reuse existing
field-equality, repeated-copy, normalization, and dimension coverage.

This is intended to preserve numerical behavior. Follow the repository's
internal-refactor policy for unchanged numerical paths: run the DMRG runner
and check the Nkeep=20 and 50 energies against
`-82.37493703992092` and `-82.39703156828021`. Tests of newly exposed catalogue
validation/error semantics are a separate behavior change and should be
limited to those contracts. Do not add a broad new integration suite or run a
Telum documentation build.

Completion means family rules exist as data, every fundamental request follows
that data, and both validation procedures consume the catalogue structures.
Moving methods between files alone is not completion. A metadata-only cleanup
must not change the numerical cache namespace; if validation detects a basis
change, investigate it rather than masking it with cache invalidation.
