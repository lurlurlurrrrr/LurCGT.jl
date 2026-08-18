module LurCGT

using LinearAlgebra
using SparseArrays
using SparseArrayKit
using Combinatorics
using TensorOperations
using Nemo

include("Base.jl")

# The exported API is the symmetry interface consumed directly by Telum.
# Cache/storage and CGT-construction internals remain accessible as `LurCGT.name`.
export Z, U1, SU, SO, Sp, G2
export Symmetry, AbelianSymm, NonabelianSymm

export add_qn, decompose_irop, decompose_space, detect_1j, dimension
export get_CGTom, get_IROP, get_conj_perm, get_dualq
export getNsave_CGTperm, getNsave_CGTSVD, getNsave_CGTQR
export getNsave_Xsymbol, getNsave_Conjperm
export getNsave_omlist, getNsave_validout, isabelian
export nzops, remove_zeros, totxt, transf_basis!
export to_float

end
