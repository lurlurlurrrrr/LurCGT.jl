# Part 3: Extend numerical irreps to all exceptional families

Status: Parts 3a and 3b are implemented; Part 3c remains planned. Parts 1 and 2
provide tensor-word vertices and Demazure lowering schedules. This document
is the overview and shared policy.

## Implementation parts

1. [Part 3a: Minuscule irreps and fundamental construction recipes](crystal_basis_minuscule_plan.md).
   Build minuscule seeds (including E6, E7, and classical spin), hard-code the
   recipes for all fundamental irreps, and generalize `get_irep_frommw` to
   select nonmaximal tensor components without requiring an existing output irrep.
2. [Part 3b: Irreps for new symmetries using virtual crystals](crystal_basis_virtual_plan.md).
   Obtain B-type, G2, and F4 target crystals through virtual recipes, then
   construct numerical irreps directly from target strings with binomial norms
   and integer lowering coefficients. Use the prescribed F4 zero-weight metric
   and validate sparse self-commutators and symmetry-wide coefficient agreement.
   This replaces the ambient numerical-word construction. Odd-SO spin seeds,
   including SO(3), continue to use the minuscule method.
3. [Part 3c: E8 adjoint crystal and numerical seed](crystal_basis_adjoint_plan.md).
   Build the adjoint crystal and exact adjoint matrices, then enable the
   incremental recipes for E8 fundamentals.

Parts 3b and 3c both consume machinery from 3a; neither depends on the other.
Together they provide a construction route for G2, F4, E6, E7, and E8.
Numerical backend support requires validated matrices as well as crystal
vertices. Large E8 recipes need not be materialized eagerly.

## Existing foundation

[Part 1](crystal_basis_irrep_plan.md) supplies `CrystalVertex` and tensor-word
operators. [Part 2](crystal_basis_string_plan.md) supplies fixed reduced words,
frozen source counts, first-discovery numerical vectors, and insertion-ordered
columns. Do not reintroduce string-coordinate extraction or a separate rank-check
gate. Retain existing exact numerical checks and let construction failures
propagate.

G2, F4, E6, and E7 now have seed and product numerical construction. E8
numerical construction remains deferred to Part 3c. Check Cartan data, node
numbering, longest-element words, and duality against Part 3a's conventions. In particular, E6 is not covered by a blanket self-duality
rule. Keep `Irep` and its sparse weight-block storage.

## Shared numerical and integration requirements

- Preserve exact Gram matrices, coefficient solves, normalization, and explicit
  integrality/overflow handling. Minuscule edge coefficients and virtual crystal
  operators have different numerical meanings; follow each part's construction.
- Make actual fundamental labels, bootstrap dependencies, and label splitting
  consistent. For product component extraction, obtain numerical highest vectors
  from product raising kernels;
  a crystal highest tuple is not generally a classical highest pure tensor.
- Maintain fixed factor/node order and deterministic component-copy selection.
  Enumerate requested components using staged products and precomputed weight
  ranges rather than full tensor Cartesian products or dense Kronecker matrices.
- Document the requested Spin/Lie-algebra interpretation of SO labels when
  enabling spin modules. Extend opposition duality and conjugation together, including E6 and
  odd-rank D spin-node exchange.
- Validate dimensions, weight multiplicities, generator relations, metric
  compatibility, and existing CG/recoupling identities for enabled backends.
  Preserve w-matrix column isometry and downstream Telum invariants.
- Version memory and SQLite caches consistently for every basis change. The
  unchanged minuscule namespace is `minuscule-v2`; B/G2/F4 use
  `crystal-seed-v1`. Invalidate affected irrep and dependent
  CG/recoupling databases when changing conventions. Never mix old/new bases.
- Inspect type stability at payload barriers and measure allocation/runtime
  impact where benchmarks exist. Keep heavy per-weight writes independent.

## Shared implementation policy

Use ordinary functions and concrete tuples/arrays with the existing types.
Additional types require justification and approval. Every new named Julia
function needs a docstring. Follow repository guidance for minimal regression
coverage in existing tests and the prescribed DMRG-only validation of pure
internal refactors. Do not run documentation builds unless requested.

## References

The three plans preserve the seed constructions from the original integration
proposal, based on Bump and Schilling, *Crystal Bases: Representations and
Combinatorics* (2017), chapter 5. See each plan for the relevant sections and
[part 1's provenance](crystal_basis_irrep_plan.md#references-and-provenance).
Numerical realization, storage, and cache requirements are implementation
requirements; the book's crystal constructions do not themselves provide every
classical matrix realization.
