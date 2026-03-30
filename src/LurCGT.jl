module LurCGT

using LinearAlgebra
using SparseArrays
using SparseArrayKit
using Combinatorics
using TensorOperations
using Nemo

include("Base.jl")

export Z, U1, SU, SO, Sp, G2
export Symmetry, AbelianSymm, NonabelianSymm

export Irep, CGT, Fsymbol, Rsymbol, Xsymbol, CGTperm, CGTSVD, CG3Flip
export OMList, ValidOuts, FTree

export getNsave_irep, getNsave_cg3, getNsave_Fsymbol, getNsave_Rsymbol
export getNsave_omlist, getNsave_validout, getNsave_CGTperm, getNsave_CGTSVD, getNsave_Xsymbol, getNsave_1jsym

export load_cg3_float
export add_qn, isabelian, nzops, get_dualq, totxt, dimension

export contract_newcg3, contract_om, contract_CG3s
export FTree2arr, contract_arrs, get_canonical_basis
export decompose_space, get_IROP, decompose_irop
export remove_zeros, get_CGTom
export detect_1j, get_conj_perm, transf_basis!
export update_fixedint_irrep_catalog, fixedint_dimension_chunks
export fixedint_canonical_pairs, run_fixedint_cgt_chunk
export fixedint_catalog_path, fixedint_chunk_result_path
export collect_fixedint_plot_cells

end
