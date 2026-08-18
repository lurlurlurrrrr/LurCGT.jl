# LurCGT.jl

LurCGT is a Julia package for computing Clebsch-Gordan coefficients and related
representation-theory data for non-Abelian symmetries. Clebsch-Gordan
coefficients describe how tensor products of irreducible representations split
into irreducible components; they are a fundamental ingredient in
symmetry-aware tensor-network calculations.

LurCGT constructs this data with integer arithmetic, rather than floating-point
arithmetic, so its computations are exact. This avoids round-off errors and
makes the generated symmetry data reproducible, which is particularly useful
when composing Clebsch-Gordan tensors and derived F/R/X symbols.

LurCGT is the Clebsch-Gordan backend of the non-Abelian tensor network library
Telum. It provides the symmetry types, irreducible representations,
Clebsch-Gordan tensors, F/R/X-symbol machinery, decomposition helpers, and
SQLite-backed storage that Telum uses to perform symmetry-aware tensor-network
operations.

Most users should use LurCGT through Telum: a normal user may not need to read
or call the functions in LurCGT directly. The public API is also available for
Telum developers and for users who need to work with symmetry data directly.

## Quick Start

After installing the package, import it and choose a symmetry family. For
example, this queries basic information about the fundamental representation of
`SU(2)`:

```julia
using LurCGT

S = SU{2}
fundamental = (1,)          # highest-weight label of the 2-dimensional irrep

dimension(S, fundamental)  # 2
get_dualq(S, fundamental)  # (1,), since all SU(2) irreps are self-dual.

# Couple two spin-1/2 irreps to the spin-0 irrep, then materialize it.
exact_cgt = LurCGT.getNsave_cg3(SU{2}, BigInt, ((1,), (1,)), [(0,)])[(0,)]
cgt, qlabels, directions = to_float(exact_cgt)  # Float64 by default
size(cgt)  # (2, 2, 1, 1): two inputs, one output, one OM channel
```

For a symmetry-aware tensor-network calculation, use Telum; it calls this
backend as needed. Direct LurCGT calls use symmetry family *types* (such as
`SU{2}`) and q-label tuples in the conventions of that symmetry family.

## Installation

Install the registered package from the General registry with:

```julia
using Pkg
Pkg.add("LurCGT")
```

For local development:

```julia
using Pkg
Pkg.develop(path=".")
```

## Testing

```julia
using Pkg
Pkg.test()
```

## Public API

The public API is intended for Telum integrations and for users working
directly with symmetry data. Available symmetry types are
Abelian: Z{N}, U1
Non-Abelian: SU{N}, SO{N} (without spin representation), Sp{2N}, G2
Irreducible representations are identified by q-label tuples; their exact
convention depends on the selected symmetry family.

A Clebsch-Gordan tensor (CGT) is the change-of-basis data that relates a tensor
product of irreducible representations to its decomposition into irreducible
components. A CGT can have multiple admissible decomposition paths; these are
distinguished by an outer-multiplicity index. LurCGT represents a CGT with
sorted tuples of upper and lower q-labels, which label its outgoing and
incoming legs, respectively.

For more information on CGTs and X-symbols, see the papers available at
[SciPost submission 202405_00027v2](https://scipost.org/submissions/scipost_202405_00027v2/)
and [Physical Review Research 2, 023385](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.023385).

The `getNsave_*` functions load a previously generated result when available,
or compute it and cache it for later use. The higher-level cached operations
are:

- `getNsave_CGTperm`: outer-multiplicity basis transform for a CGT leg permutation.
- `getNsave_Conjperm`: outer-multiplicity permutation needed for CGT conjugation.
- `getNsave_CGTSVD` and `getNsave_CGTQR`: CGT-side basis changes for SVD- and QR-style splits.
- `getNsave_Xsymbol`: basis transform for contracting two CGT canonical bases.

```julia
using LurCGT

?getNsave_Xsymbol
?getNsave_CGTSVD
```
