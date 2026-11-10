"""
    LurCGT

Exact Clebsch-Gordan and recoupling support for the symmetry-adapted tensor
library. Its public API accepts symmetry *types* such as `SU{2}` or `U1`, not
instances; q-labels are integer tuples whose width is `nzops(S)`.
"""
module LurCGT

using LinearAlgebra
using SparseArrays
using SparseArrayKit
using Combinatorics
using TensorOperations
using Nemo

"""
    VERBOSE::Base.RefValue{Int}

Global default verbosity for LurCGT operations. Set `VERBOSE[]` to `0` for
quiet operation or a positive integer to enable progressively detailed output.
An individual function's `verbose` keyword still overrides this setting.
"""
const VERBOSE = Ref{Int}(0)

include("Base.jl")

# The exported API is the symmetry interface consumed directly by Telum.
# Cache/storage and CGT-construction internals remain accessible as `LurCGT.name`.
export Z, U1, SU, SO, Sp, G2, F4, E6, E7, E8
export Symmetry, AbelianSymm, NonabelianSymm

export add_qn, decompose_irop, decompose_space, dimension
export get_CGTom, get_IROP, get_conj_perm, get_dualq
export getNsave_CGTperm, getNsave_CGTSVD, getNsave_CGTQR
export getNsave_Xsymbol, getNsave_Conjperm
export getNsave_omlist, getNsave_validout, isabelian
export nzops, remove_zeros, totxt, transf_basis!
export to_float
export VERBOSE

end
