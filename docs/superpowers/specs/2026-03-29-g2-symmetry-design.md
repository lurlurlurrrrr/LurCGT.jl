# G2 Symmetry Design

## Goal

Add core `G2` symmetry support to `LurCGT` so the package can construct the defining irrep, derive the second fundamental irrep from the alternating square of the defining irrep, and generate higher irreps through the existing crystal-based machinery.

## Scope

In scope:
- Implement the symmetry-specific hooks required by the existing nonabelian irrep/crystal pipeline.
- Expose `G2` from the public module API.
- Add focused tests for core representation generation and dimensions.

Out of scope for this pass:
- `Int128` frontier or sum-slice helper updates.
- New generic abstractions for exceptional Lie groups.
- Refactoring shared tableau machinery beyond what `G2` strictly needs.

## Existing Integration Points

`LurCGT` already models nonabelian symmetries through symmetry-specific methods consumed by shared code in:
- `src/irep.jl` for defining/trivial irreps and alternating products.
- `src/clebsch.jl` for building irreps from maximal weights and tensor products.
- `src/Symmetry/crystal.jl` for highest-weight tableaux and crystal traversal.

`G2` is already declared as a symmetry type in `src/Symmetry/Symmetry.jl` and included there, but `src/Symmetry/G2.jl` is currently empty.

## G2 Data Model

Treat `G2` as a fixed rank-2 nonabelian symmetry:
- `defirepdim(G2) = 7`
- `nlops(G2) = 2`
- `nzops(G2) = 2`
- `totxt(G2) = "G2"`
- `isvalidsymm(G2) = true`

### Defining Crystal Alphabet

Use the defining crystal states in the order:

`1, 2, 3, 0, -3, -2, -1`

The crystal lowering operators are:
- `f1`: `1 -> 2`, `3 -> 0`, `0 -> -3`, `-2 -> -1`
- `f2`: `2 -> 3`, `-3 -> -2`

The inverse maps define the raising operators through the existing shared helper.

### Defining Weights

Assign the defining weights so that:
- the maximal weight of the defining irrep `(1, 0)` is `(1, 1)`
- `f1` lowers weight by `(2, 0)`
- `f2` lowers weight by `(-3, 1)`

This fixes the seven defining weights as:

| state | weight |
| --- | --- |
| `1` | `(1, 1)` |
| `2` | `(-1, 1)` |
| `3` | `(2, 0)` |
| `0` | `(0, 0)` |
| `-3` | `(-2, 0)` |
| `-2` | `(1, -1)` |
| `-1` | `(-1, -1)` |

These weights will be returned by `getsz_def(G2)` and `getsz_def_vec(G2)` in the crystal order needed by the shared code.

### Maximal-Weight / qlabel Conversion

Use the Dynkin-label to maximal-weight map:

`qlab2mwz(G2, (a, b)) = (a, a + 2b)`

This satisfies:
- `(1, 0) -> (1, 1)`
- `(0, 1) -> (0, 2)`

The inverse conversion is:
- `a = z1`
- `b = (z2 - z1) / 2`

So `getqlabel(G2, (z1, z2))` must assert that `z2 - z1` is even and return `(z1, (z2 - z1) ÷ 2)`.

## Tableau / Crystal Construction

The current rank-2 tableau shape builder in `src/Symmetry/crystal.jl` already converts `qlabel = (a, b)` into shape `[a + b, b]`, which matches the intended highest-weight tableau shape.

Because of that, the first pass should reuse the existing tableau construction unchanged and only provide `G2` implementations for:
- `charlist`
- `crystal_chars_map`
- `get_fops_std`
- `mw_column`
- `less_weight`

For `G2`, the highest-weight tableau for `qlabel = (a, b)` must use shape `[a + b, b]` with columns `[1]` and `[1, 2]`. The entries must strictly increase within each column, so the implementation should treat `[1, 2]` as the two-box maximal column rather than `[2, 1]`.

## Fundamental Irreps

The shared `get_fundamental_ireps` path saves:
- the trivial irrep
- the defining irrep
- alternating products `Λ^n(def)` for `n = 2:nzops(S)`

For `G2`, `nzops(G2) = 2`, so this automatically computes only `Λ^2(7)`.

The design relies on the given fact that the second fundamental irrep `(0, 1)` is obtained as the irreducible component extracted from `Λ^2(7)`. The existing `def_altprod` plus `get_irep_frommw` pipeline should therefore be sufficient with no special-case code beyond the G2 symmetry hooks.

## Public API Change

Export `G2` from `src/LurCGT.jl` alongside the existing symmetry types so callers can request:
- `getNsave_irep(G2, BigInt, (1, 0))`
- `getNsave_irep(G2, BigInt, (0, 1))`

## Testing Strategy

Add a focused `G2` testset to `test/runtests.jl` that covers:

1. Symmetry metadata:
- `isvalidsymm(G2)`
- `defirepdim(G2) == 7`
- `nlops(G2) == 2`
- `nzops(G2) == 2`

2. Defining crystal data:
- exact defining weights
- exact `getdz(G2, 1)` and `getdz(G2, 2)`
- exact lowering maps in `get_fops_std(G2)`

3. Weight conversion:
- `qlab2mwz(G2, (1, 0)) == (1, 1)`
- `qlab2mwz(G2, (0, 1)) == (0, 2)`
- round-trip `getqlabel(G2, qlab2mwz(G2, q)) == q` for a small sample

4. Fundamental irreps:
- `dimension(getNsave_irep(G2, BigInt, (1, 0))) == 7`
- `dimension(getNsave_irep(G2, BigInt, (0, 1))) == 14`

5. Dimension regression:
- enumerate all nonnegative `qlabel = (a, b)` whose Weyl-dimension value is at most `13090`
- assert the sorted list of computed irrep dimensions matches exactly:

`1, 7, 14, 27, 64, 77, 77, 182, 189, 273, 286, 378, 448, 714, 729, 748, 896, 924, 1254, 1547, 1728, 1729, 2079, 2079, 2261, 2926, 3003, 3289, 3542, 4096, 4914, 4928, 4928, 5005, 5103, 6630, 7293, 7371, 7722, 8372, 9177, 9660, 10206, 10556, 11571, 11648, 12096, 13090`

This regression gives broad coverage without tying the design to the unrelated `Int128` helper utilities.

## Risks And Mitigations

Risk: the existing tableau machinery may rely on admissibility rules that happen to hold for `SU`/`Sp`/`SO` but not for `G2`.

Mitigation: the first implementation will validate the generated dimensions against the supplied list. If dimensions diverge, the likely next step is a `G2`-specific tableau admissibility or lowering-path refinement rather than a full rewrite.

Risk: `getqlabel` may be called on weights that are not dominant.

Mitigation: keep the same pattern as the current symmetry modules: enforce only the algebraic integrality conditions needed for inversion and let the existing irrep-generation flow determine which dominant outputs are valid.

## Implementation Summary

The implementation should be a narrow extension:
- fill in `src/Symmetry/G2.jl`
- export `G2` in `src/LurCGT.jl`
- add a dedicated `G2` testset in `test/runtests.jl`

No shared-architecture changes are planned unless the dimension regression proves they are necessary.
