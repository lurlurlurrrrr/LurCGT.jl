# Getting started

Install LurCGT from Julia's General registry:

```julia
using Pkg
Pkg.add("LurCGT")
```

Symmetry families are represented by Julia types. For example, `SU{2}` uses a
one-element Dynkin q-label tuple: `(1,)` is the spin-1/2 irrep and `(0,)` is
the trivial irrep.

```julia
using LurCGT

S = SU{2}
fundamental = (1,)
dimension(S, fundamental)  # 2
get_dualq(S, fundamental)  # (1,)
```

## Floating sparse Clebsch-Gordan tensors

`to_float` converts an exact Clebsch-Gordan tensor to floating sparse storage.
This example obtains an exact three-leg CGT from LurCGT's cache, generating it
when necessary, then materializes it for numerical work.

```julia
exact_cgt = LurCGT.getNsave_cg3(SU{2}, BigInt, ((1,), (1,)), [(0,)])[(0,)]
cgt, qlabels, directions = to_float(exact_cgt)  # Float64 by default

size(cgt)       # (2, 2, 1, 1)
qlabels         # ((1,), (1,), (0,))
directions      # ('+', '+', '-')
```

The first three axes are the irrep basis axes in the requested input/input/output
order. The trailing axis enumerates outer-multiplicity channels. For the
spin-0 channel of two spin-1/2 irreps, there is one such channel.

Use `normalize=true` when the numeric tensor must include LurCGT's stored
outer-multiplicity normalization factors.
