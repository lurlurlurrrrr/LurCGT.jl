# Getting started

Install LurCGT from Julia's General registry:

```julia
using Pkg
Pkg.add("LurCGT")
```

Symmetry families are represented by Julia types. For example,
``\mathrm{SU}(2)`` uses a one-element Dynkin q-label tuple: `(1,)` is the
spin-1/2 irrep and `(0,)` is the trivial (1-dimensional) irrep.

```julia
using LurCGT

S = SU{2}
fundamental = (1,) # 2S
dimension(S, fundamental)  # 2
get_dualq(S, fundamental)  # (1,)
```

## Floating sparse Clebsch-Gordan tensors

`to_float` converts an exact Clebsch-Gordan tensor (CGT) to a floating-point
sparse array. The coefficients are generated exactly, then converted to the
requested floating-point type. This makes the sparse array convenient for
inspecting coefficients or using them in numerical code.

The following example couples two spin-``1 / 2`` ``\mathrm{SU}(2)`` irreps:
``S = 1 / 2 \otimes S = 1 / 2 \to S = 0 \oplus S = 1``. The two spin-``1 / 2``
irreps are the input legs, and the singlet and triplet irreps are the output
legs.

```@example spin_half_coefficients
using LurCGT

cg3s = LurCGT.getNsave_cg3(SU{2}, BigInt, ((1,), (1,)), [(0,), (2,)])
singlet, qlabels, directions = to_float(cg3s[(0,)])
triplet, _, _ = to_float(cg3s[(2,)])

size(singlet)  # (2, 2, 1, 1)
qlabels         # ((1,), (1,), (0,))
directions      # ('+', '+', '-')
```

The first three axes are the irrep basis axes in the requested input/input/output
order. The trailing axis enumerates outer-multiplicity channels. It has length
one in this two-input, one-output ``\mathrm{SU}(2)`` example. See
[CGT and related symbols](@ref) for the general case.

Use `normalize=true` when the numeric tensor must include LurCGT's stored
outer-multiplicity normalization factors; it is the default.

## Singlet and triplet coefficients

For the basis order ``\lvert \uparrow \rangle, \lvert \downarrow \rangle`` on
each input and ``\lvert 1, 1 \rangle, \lvert 1, 0 \rangle,
\lvert 1, -1 \rangle`` on the triplet output, the sparse-array values agree
with the usual coupled states:

```math
\lvert 0, 0 \rangle =
  \frac{\lvert \uparrow\downarrow \rangle -
        \lvert \downarrow\uparrow \rangle}{\sqrt{2}},
\qquad
\begin{aligned}
\lvert 1, 1 \rangle &= \lvert \uparrow\uparrow \rangle, \\
\lvert 1, 0 \rangle &=
  \frac{\lvert \uparrow\downarrow \rangle +
        \lvert \downarrow\uparrow \rangle}{\sqrt{2}}, \\
\lvert 1, -1 \rangle &= \lvert \downarrow\downarrow \rangle.
\end{aligned}
```

```@example spin_half_coefficients
singlet_coeffs = Array(singlet[:, :, 1, 1])
triplet_coeffs = Array(triplet[:, :, :, 1])

expected_singlet = [0 1; -1 0] / sqrt(2)
expected_triplet = cat(
    [1 0; 0 0],
    [0 1; 1 0] / sqrt(2),
    [0 0; 0 1];
    dims=3,
)

(
    singlet = singlet_coeffs,
    triplet = triplet_coeffs,
    singlet_matches = singlet_coeffs ≈ expected_singlet,
    triplet_matches = triplet_coeffs ≈ expected_triplet,
)
```

The final two values are `true`, confirming that the floating sparse tensors
use these singlet and triplet coefficients. The relative signs fix LurCGT's
chosen basis convention.
