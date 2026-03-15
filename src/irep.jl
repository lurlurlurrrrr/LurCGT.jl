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
function def_altprod(irep::Irep{S, NL, NZ, RT}, n::Int) where {S, NL, NZ, RT}
    # dimension of the original representation
    d_orig = dimension(irep)
    Sz_orig = sorted_zvals(irep.Sz)
    # Is it defining representation?
    @assert irep.qlabel == ntuple(i -> i == 1 ? 1 : 0, NZ)

    # Get a basis of alternating product space and weight of them
    # TODO: For other kind of symmetry such as Sp(2n), this should be generalized
    new_basis = collect(combinations(1:d_orig, n))
    weight_dict = Dict{NTuple{NZ, Int}, Vector{Vector{Int}}}()
    norm2s_new = RT[]

    Sl_numform = get_Sl_numform(irep)

    for nb in new_basis
        newnorm = RT(1)
        new_zval = NTuple{NZ, Int}(0 for _ in 1:NZ)
        for bidx in nb
            new_zval = new_zval .+ Sz_orig[bidx]  # accumulate z-operator values
        end
        if !haskey(weight_dict, new_zval)
            weight_dict[new_zval] = Vector{Vector{Int}}()
        end
        push!(weight_dict[new_zval], nb)
        push!(norm2s_new, newnorm)
    end

    Sz_new = Dict{NTuple{NZ, Int}, Tuple{Int, Int}}()
    zlst = sorted_zvals(weight_dict)
    mw = zlst[1]

    nd = 1
    for k in zlst
        nstates = length(weight_dict[k])
        Sz_new[k] = (nd, nd + nstates - 1)
        nd += nstates
    end

    display(Sz_new)
    
    # Get new lowering operator of the new representation
    Sl_new = Tuple(Dict{NTuple{NZ, Int}, SparseMatrixCSC{RT}}() for _=1:NL)
    for lop=1:NZ
        op_numform = Sl_numform[lop]
        dz = getdz(S, lop) # Get the change of weight for the i-th lowering op
        # For each weight of resulting representation
        for ow in keys(weight_dict) # old weight
            # Get the basis corresponding to the weight
            nstates = length(weight_dict[ow])
            nw = ow .- dz  # new weight
            if !haskey(weight_dict, nw) continue end

            m, n = length(weight_dict[nw]), length(weight_dict[ow])
            spmat = spzeros(RT, m, n)
            for (i, obasis) in enumerate(weight_dict[ow])
                for (j, nbasis) in enumerate(weight_dict[nw])
                    spmat[j, i] = get_matelem(op_numform, obasis, nbasis)
                    # If nbasis can be obtained from obasis by applying lowering operator,
                end
            end
            Sl_new[lop][ow] = spmat
        end
    end
    inner_prod = get_identity_innerprod(S, RT, Sz_new)
    return Sl_new, Sz_new, mw, inner_prod
end

# Assume that rep is defining representation
function get_Sl_numform(rep::Irep{S, NL, NZ, RT}) where {S<:NonabelianSymm, NL, NZ, RT}
    Sl_numform = Vector{Dict{Int, Int}}()
    for lop=1:NL
        dz = getdz(S, lop)
        dict_op = Dict{Int, Int}()
        for (ow, v) in rep.Sl[lop]
            # Assume all weight space is one-dimensional
            i1, i2 = rep.Sz[ow]; @assert i1 == i2
            # Assume that the matrix is 1*1, and its element is 1
            @assert size(v, 1) == 1 && size(v, 2) == 1 && v[1, 1] == 1
            nw = ow .- dz
            j1, j2 = rep.Sz[nw]; @assert j1 == j2
            dict_op[i1] = j1
        end
        push!(Sl_numform, dict_op)
    end
    return Sl_numform
end

function get_matelem(op_numform::Dict{Int, Int}, 
        obasis::Vector{Int}, 
        nbasis::Vector{Int})
    # Check if nbasis can be obtained from obasis by applying lowering operator
    # The elements of obasis and nbasis are all distinct
    @assert length(obasis) == length(nbasis)
    s1 = Set(obasis); s2 = Set(nbasis)
    s12 = setdiff(s1, s2); s21 = setdiff(s2, s1)
    if length(s12) != 1 || length(s21) != 1 return 0 end

    a, b = first(s12), first(s21)
    if get(op_numform, a, -1) != b return 0 end

    i = findfirst(x -> x == a, obasis)
    obasis_copy = obasis[:]
    obasis_copy[i] = b  # replace a with b
    perm = sortperm(obasis_copy)
    return permutation_sign(perm)
end

# What element is different between two vectors
# The output is three numbers i, a, b
# i is the index of the element in v1 that is not in v2
# a is the element of v1 that is not in v2
# b is the element of v2 that is not in v1
function find_diff(v1::Vector{Int}, v2::Vector{Int})
    s1 = Set(v1); s2 = Set(v2)
    s12 = setdiff(s1, s2)
    s21 = setdiff(s2, s1)
    @assert length(s12) == 1 && length(s21) == 1
    a, b = first(s12), first(s21)
    return findfirst(x -> x == a, v1), a, b
end

function Base.show(io::IO, rep::Irep{S}) where {S<:Symmetry}
    print(io, "$(totxt(S)) Irep: qlabel=$(rep.qlabel), dimension=$(rep.dimension)")
end

symm(rep::Irep{S}) where {S<:Symmetry} = S