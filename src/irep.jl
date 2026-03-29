# Sl : lowering operators, Sz : z-operator
# norm2s : square norm values of basis element.
# For defining irep, values of innerprod are all [[1]], 1*1 matrix with entry 1
# NL(NZ) means number of lowering(z-) operators

# For larger irep, they are positive integer
# and matrix element of Sl, Sz are integer w.r.t such basis
# This may not be irreducible representation, for example, product of two irreps
struct Irep{S<:Symmetry, NL, NZ, RT}
	Sl::NTuple{NL, Dict{NTuple{NZ, Int}, SparseMatrixCSC{Int}}}
	Sz::Dict{NTuple{NZ, Int}, Tuple{Int, Int}}
    innerprod::Dict{NTuple{NZ, Int}, Matrix{RT}}
    inv_innerprod::Dict{NTuple{NZ, Int}, Tuple{Matrix{RT}, Rational{BigInt}}}

	qlabel::NTuple{NZ, Int}
    dimension::Int
    # Cached memory footprint in bytes, computed once at construction for LRU cache eviction by size
    size_byte::Int
    
    function Irep{S, NL, NZ, RT}(Sl, Sz, innerprod, inv_innerprod, qlabel, dimension, size_byte::Int=0) where {S, NL, NZ, RT}
        if size_byte == 0
            obj = new{S, NL, NZ, RT}(Sl, Sz, innerprod, inv_innerprod, qlabel, dimension, 0)
            size_byte = Base.summarysize(obj)
        end
        new{S, NL, NZ, RT}(Sl, Sz, innerprod, inv_innerprod, qlabel, dimension, size_byte)
    end
end

Base.isless(r1::Irep{S, NL, NZ, RT}, r2::Irep{S, NL, NZ, RT}) where {S, NL, NZ, RT} =
r1.qlabel < r2.qlabel

# Custom pretty printing for Irep type
totxt(rep::Irep{S, NL, NZ, RT}) where {S<:NonabelianSymm, NL, NZ, RT} = 
    "irep($(totxt(S)), $(rep.qlabel))"

# get max weight from other fields and complete Irep struct
function Irep(::Type{S}, 
    ::Type{RT}, 
    Sl, 
    Sz, 
    innerprod, 
    inv_innerprod,
    dim) where {S<:Symmetry, RT<:Number}
    maxz = sorted_zvals(Sz)[1]
    @assert nlops(S) == length(Sl)
    @assert nzops(S) == length(maxz)
    Irep{S, nlops(S), nzops(S), RT}(Sl, Sz, innerprod, 
        inv_innerprod, getqlabel(S, maxz), dim)
end

dimension(::Type{S}, q::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ} =
dimension(getNsave_irep(S, BigInt, q))

dimension(::Type{S}, q::NTuple{NZ, Int}) where {S<:AbelianSymm, NZ} = 1

dimension(rep::Irep) = rep.dimension

function dim_from_innerprod(innerprod::Dict{NTuple{NZ, Int}, Matrix{RT}}) where {NZ, RT}
    dim = 0
    # Sum the dimensions of each inner product matrix)
    for (_, mat) in innerprod dim += size(mat, 1) end 
    return dim
end

function get_identity_innerprod(::Type{S}, 
    ::Type{RT}, 
    Sz) where {S<:NonabelianSymm, RT<:Number}
    inner_prod = Dict{NTuple{nzops(S), Int}, Matrix{RT}}()
    for (z, (a, b)) in Sz
        sz = b - a + 1
        inner_prod[z] = sparse(RT, I, sz, sz)  # Identity matrix for inner product
    end
    return inner_prod
end

function get_identity_invinprod(::Type{S},
    ::Type{RT},
    Sz) where {S<:NonabelianSymm, RT<:Number}
    inv_innerprod = Dict{NTuple{nzops(S), Int}, Tuple{Matrix{RT}, Rational{RT}}}()
    for (z, (a, b)) in Sz
        sz = b - a + 1
        inv_innerprod[z] = (sparse(RT, I, sz, sz), Rational{BigInt}(1//1))  
    end
    return inv_innerprod
end

""" 
get defining representation
Lowering operators are stored in block sparse matrix
"""
function getdefirep(::Type{S}, ::Type{RT}) where {S<:NonabelianSymm, RT<:Number}
	# getsl : return tuple of lowering operators (each are SparseMatrix)
	# getsz : return tuple of diagonals of z-operators 
	Sl = getsl_def(S, RT); Sz = getsz_def(S)
    inner_prod = get_identity_innerprod(S, RT, Sz)
    inv_innerprod = get_identity_invinprod(S, RT, Sz)
    # For defining representation, Tmats is equal to innerprod mat
    dim = dim_from_innerprod(inner_prod)
	return Irep(S, RT, Sl, Sz, inner_prod, inv_innerprod, dim)
end

# Get trivial (1-dimensional) representation
function gettrivirep(::Type{S}, ::Type{RT}) where {S<:NonabelianSymm, RT<:Number}
    Sl = getsl_triv(S, RT); Sz = getsz_triv(S)
    inner_prod = get_identity_innerprod(S, RT, Sz)
    inv_innerprod = get_identity_invinprod(S, RT, Sz)
    # For trivial representation, Tmats is equal to innerprod mat
    dim = dim_from_innerprod(inner_prod)
    return Irep(S, RT, Sl, Sz, inner_prod, inv_innerprod, dim)
end

# Get the alternating product representation of the defining representation
# The first argument is always the defining representation
# TODO: Read the code and check if this is working correctly
# For SU(3), it works correctly
function def_altprod(irep::Irep{S, NL, NZ, RT}, n::Int) where {S, NL, NZ, RT}
    d_orig = dimension(irep)
    @assert 1 <= n <= d_orig
    @assert irep.qlabel == ntuple(i -> i == 1 ? 1 : 0, NZ)

    basis_weight, basis_local_idx, basis_norm2, local_to_global =
        get_altprod_basis_metadata(irep)
    weight_dict, basis_norms = build_altprod_basis(basis_weight, basis_norm2, n)
    basis_lookup = get_altprod_basis_lookup(weight_dict)

    Sz_new = Dict{NTuple{NZ, Int}, Tuple{Int, Int}}()
    zlst = sorted_zvals(weight_dict)
    mw = zlst[1]

    nd = 1
    for k in zlst
        nstates = length(weight_dict[k])
        Sz_new[k] = (nd, nd + nstates - 1)
        nd += nstates
    end
    Sl_new = build_altprod_lowering(
        irep,
        weight_dict,
        basis_lookup,
        basis_weight,
        basis_local_idx,
        local_to_global,
    )
    inner_prod = build_altprod_innerprod(basis_norms)
    return Sl_new, Sz_new, mw, inner_prod
end

function get_altprod_basis_metadata(irep::Irep{S, NL, NZ, RT}) where {S<:NonabelianSymm, NL, NZ, RT}
    d_orig = dimension(irep)
    basis_weight = Vector{NTuple{NZ, Int}}(undef, d_orig)
    basis_local_idx = Vector{Int}(undef, d_orig)
    basis_norm2 = Vector{RT}(undef, d_orig)
    local_to_global = Dict{NTuple{NZ, Int}, Vector{Int}}()

    for z in sorted_zvals(irep.Sz)
        a, b = irep.Sz[z]
        diag_entries = get_diagonal_entries(irep.innerprod[z])
        @assert length(diag_entries) == b - a + 1
        local_to_global[z] = collect(a:b)
        for (local_idx, global_idx) in enumerate(a:b)
            basis_weight[global_idx] = z
            basis_local_idx[global_idx] = local_idx
            basis_norm2[global_idx] = diag_entries[local_idx]
        end
    end

    return basis_weight, basis_local_idx, basis_norm2, local_to_global
end

function get_diagonal_entries(mat::Matrix{RT}) where RT
    @assert size(mat, 1) == size(mat, 2)
    diag_entries = Vector{RT}(undef, size(mat, 1))
    for j in axes(mat, 2), i in axes(mat, 1)
        if i == j
            diag_entries[i] = mat[i, j]
        else
            @assert iszero(mat[i, j])
        end
    end
    return diag_entries
end

function build_altprod_basis(
    basis_weight::Vector{NTuple{NZ, Int}},
    basis_norm2::Vector{RT},
    n::Int,
) where {NZ, RT}
    weight_dict = Dict{NTuple{NZ, Int}, Vector{Vector{Int}}}()
    basis_norms = Dict{NTuple{NZ, Int}, Vector{RT}}()

    for nb in combinations(1:length(basis_weight), n)
        new_zval = ntuple(_ -> 0, NZ)
        new_norm = one(RT)
        for bidx in nb
            new_zval = new_zval .+ basis_weight[bidx]
            new_norm *= basis_norm2[bidx]
        end
        if !haskey(weight_dict, new_zval)
            weight_dict[new_zval] = Vector{Vector{Int}}()
            basis_norms[new_zval] = RT[]
        end
        push!(weight_dict[new_zval], collect(nb))
        push!(basis_norms[new_zval], new_norm)
    end

    return weight_dict, basis_norms
end

function get_altprod_basis_lookup(weight_dict::Dict{NTuple{NZ, Int}, Vector{Vector{Int}}}) where NZ
    basis_lookup = Dict{NTuple{NZ, Int}, Dict{Tuple{Vararg{Int}}, Int}}()
    for (z, basis_vecs) in weight_dict
        basis_lookup[z] = Dict(Tuple(basis) => i for (i, basis) in enumerate(basis_vecs))
    end
    return basis_lookup
end

function build_altprod_lowering(
    irep::Irep{S, NL, NZ, RT},
    weight_dict::Dict{NTuple{NZ, Int}, Vector{Vector{Int}}},
    basis_lookup::Dict{NTuple{NZ, Int}, Dict{Tuple{Vararg{Int}}, Int}},
    basis_weight::Vector{NTuple{NZ, Int}},
    basis_local_idx::Vector{Int},
    local_to_global::Dict{NTuple{NZ, Int}, Vector{Int}},
) where {S<:NonabelianSymm, NL, NZ, RT}
    Sl_new = Tuple(Dict{NTuple{NZ, Int}, SparseMatrixCSC{RT}}() for _=1:NL)

    for lop in 1:NL
        dz = getdz(S, lop)
        for (ow, old_basis_vecs) in weight_dict
            nw = ow .- dz
            haskey(weight_dict, nw) || continue
            mat = spzeros(RT, length(weight_dict[nw]), length(old_basis_vecs))
            row_lookup = basis_lookup[nw]
            for (col, obasis) in enumerate(old_basis_vecs)
                coeffs = get_altprod_column(
                    irep,
                    lop,
                    dz,
                    obasis,
                    row_lookup,
                    basis_weight,
                    basis_local_idx,
                    local_to_global,
                )
                for (row, coeff) in coeffs
                    mat[row, col] = coeff
                end
            end
            Sl_new[lop][ow] = mat
        end
    end

    return Sl_new
end

function get_altprod_column(
    irep::Irep{S, NL, NZ, RT},
    lop::Int,
    dz::NTuple{NZ, Int},
    obasis::Vector{Int},
    row_lookup::Dict{Tuple{Vararg{Int}}, Int},
    basis_weight::Vector{NTuple{NZ, Int}},
    basis_local_idx::Vector{Int},
    local_to_global::Dict{NTuple{NZ, Int}, Vector{Int}},
) where {S<:NonabelianSymm, NL, NZ, RT}
    coeffs = Dict{Int, RT}()

    for pos in eachindex(obasis)
        src_global = obasis[pos]
        src_weight = basis_weight[src_global]
        haskey(irep.Sl[lop], src_weight) || continue

        op_block = irep.Sl[lop][src_weight]
        target_weight = src_weight .- dz
        haskey(local_to_global, target_weight) || continue

        src_local = basis_local_idx[src_global]
        target_globals = local_to_global[target_weight]
        @assert size(op_block, 1) == length(target_globals)
        @assert size(op_block, 2) >= src_local

        for target_local in eachindex(target_globals)
            coeff = op_block[target_local, src_local]
            iszero(coeff) && continue

            new_basis = copy(obasis)
            new_basis[pos] = target_globals[target_local]
            length(unique(new_basis)) == length(new_basis) || continue

            perm = sortperm(new_basis)
            sorted_basis = new_basis[perm]
            row = row_lookup[Tuple(sorted_basis)]
            coeffs[row] = get(coeffs, row, zero(RT)) + coeff * RT(permutation_sign(perm))
        end
    end

    return coeffs
end

function build_altprod_innerprod(basis_norms::Dict{NTuple{NZ, Int}, Vector{RT}}) where {NZ, RT}
    inner_prod = Dict{NTuple{NZ, Int}, Matrix{RT}}()
    for (z, norms) in basis_norms
        mat = zeros(RT, length(norms), length(norms))
        for i in eachindex(norms)
            mat[i, i] = norms[i]
        end
        inner_prod[z] = mat
    end
    return inner_prod
end

function Base.show(io::IO, rep::Irep{S}) where {S<:Symmetry}
    print(io, "$(totxt(S)) Irep: qlabel=$(rep.qlabel), dimension=$(rep.dimension)")
end

symm(rep::Irep{S}) where {S<:Symmetry} = S
