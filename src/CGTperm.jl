struct CGTperm{S<:NonabelianSymm, U, D, N, NZ}
    perm_arr::Array{Float64, 2}
    upsp::NTuple{U, NTuple{NZ, Int}}  # outgoing spaces
    dnsp::NTuple{D, NTuple{NZ, Int}}  # incoming spaces
    perm::NTuple{N, Int}  # permutation of N legs
    size_byte::Int

    function CGTperm{S, U, D, N, NZ}(perm_arr, upsp, dnsp, perm, size_byte::Int=0) where {S<:NonabelianSymm, U, D, N, NZ}
        if size_byte == 0
            obj = new{S, U, D, N, NZ}(perm_arr, upsp, dnsp, perm, 0)
            size_byte = Base.summarysize(obj)
        end
        new{S, U, D, N, NZ}(perm_arr, upsp, dnsp, perm, size_byte)
    end
end

function perm_FTrees!(FTrees_dict::Dict{NTuple{NZ, Int}, Vector{FTree{S, N, NZ}}}, 
    perm::NTuple{M, Int}) where {S<:NonabelianSymm, N, M, NZ}

    for (_, ftree_vec) in FTrees_dict
        for i in 1:length(ftree_vec)
            #println(ftree_vec[i])
            #println(typeof(ftree_vec[i]))
            ftree_vec[i] = permute!(ftree_vec[i], perm)
        end
    end
end

function fill_CGTperm_matrix!(::Type{S},
    CGTperm_arr::Array{Float64, 2},
    CGT_oms::CGTom{S, NZ},
    CGT_FTrees_up::Dict{NTuple{NZ, Int}, Vector{FTree{S, U, NZ}}},
    CGT_FTrees_dn::Dict{NTuple{NZ, Int}, Vector{FTree{S, D, NZ}}}) where {S<:NonabelianSymm, U, D, NZ}

    for i in 1:CGT_oms.totalOM
        csp, upidx, dnidx = getominfo(CGT_oms, i)
        up_ftree = CGT_FTrees_up[csp][upidx]
        dn_ftree = CGT_FTrees_dn[csp][dnidx]
        omlist_id = getNsave_omlist(S, (csp,) ,csp)
        id_ftree1 = create_unit_FTree(omlist_id, 1)
        result = contN2canonical(up_ftree, id_ftree1, copy(id_ftree1), dn_ftree, (U+1,), (1,))
        CGTperm_arr[:, i] = to_vector(result, CGT_oms)
    end
end

function remove_zeros(::Type{S}, 
    spaces::Tuple{}, 
    perm::Tuple{}) where {S<:NonabelianSymm}

    NZ = nzops(S)
    zq = Tuple(0 for _=1:NZ)
    return (zq,), (1,)
end

function remove_zeros(::Type{S}, spaces::NTuple{M, NTuple{NZ, Int}}, 
    perm::NTuple{M, Int}) where {S<:NonabelianSymm, M, NZ}

    @assert nzops(S) == NZ
    zq = Tuple(0 for _=1:NZ)
    nz = findlast(i->spaces[i]==zq, 1:M)
    if isnothing(nz) return spaces, perm end
    spcs = spaces[nz+1:end]; perm_ = perm[nz+1:end] .- nz
    if isempty(spcs) spcs, perm_ = (zq,), (1,) end
    return spcs, perm_
end

getNsave_CGTperm(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    perm::NTuple{N, Int};
    save=true) where {S<:AbelianSymm, U, D, NZ, N} = nothing


function getNsave_CGTperm(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    perm::NTuple{N, Int};
    save=true) where {S<:NonabelianSymm, U, D, NZ, N}

    @assert issorted(upsp) && issorted(dnsp)
    perm_up = perm[1:U]; perm_dn = Tuple(i-U for i in perm[U+1:end])
    upsp_, perm_up_ = remove_zeros(S, upsp, perm_up)
    dnsp_, perm_dn_ = remove_zeros(S, dnsp, perm_dn)
    perm_ = (perm_up_..., [i+length(upsp_) for i in perm_dn_]...)
    # If permutation becomes identity after removing zeros, return nothing
    if issorted(perm_) return nothing end
    getNsave_CGTperm_std(S, upsp_, dnsp_, perm_; save=save)
end

function getNsave_CGTperm_std(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    perm::NTuple{N, Int};
    save=true) where {S<:NonabelianSymm, U, D, NZ, N}

    spaces = (upsp..., dnsp...)
    # permute spaces must be the same as original spaces
    @assert spaces == Tuple(spaces[i] for i in perm)
    for i in 1:U @assert perm[i] <= U end
    for i in 1:D @assert perm[U+i] > U end
    perm_up = perm[1:U]; perm_dn = Tuple(i-U for i in perm[U+1:end])
    @assert NZ == nzops(S)
    
    # Try to load from HDF5, use it if exists
    loaded = load_CGTperm_sqlite(S, upsp, dnsp, perm)
    if !isnothing(loaded) return loaded end

    CGT_oms = get_CGTom(S, upsp, dnsp)
    CGTperm_arr = zeros(Float64, CGT_oms.totalOM, CGT_oms.totalOM)

    CGT_FTrees_up = Dict{NTuple{NZ, Int}, Vector{FTree{S, U, NZ}}}()
    CGT_FTrees_dn = Dict{NTuple{NZ, Int}, Vector{FTree{S, D, NZ}}}()
    for csp in CGT_oms.spaces
        fill_FTrees!(CGT_FTrees_up, upsp, csp)
        fill_FTrees!(CGT_FTrees_dn, dnsp, csp)
    end

    perm_FTrees!(CGT_FTrees_up, perm_up)
    perm_FTrees!(CGT_FTrees_dn, perm_dn)

    fill_CGTperm_matrix!(S, CGTperm_arr, CGT_oms, CGT_FTrees_up, CGT_FTrees_dn)
    # If outer multiplicity == 1, no need to normalize
    # since division and multiplication are canceled out
    if CGT_oms.totalOM > 1
        CGT_norms = get_CGT_norms(S, upsp, dnsp, CGT_oms, false; use1j=false)
        div_along_dim!(CGTperm_arr, CGT_norms, 1)
        mul_along_dim!(CGTperm_arr, CGT_norms, 2)
    end
    CGTperm_obj = CGTperm{S, U, D, N, NZ}(CGTperm_arr, upsp, dnsp, perm)

    if save save_CGTperm_sqlite(S, CGTperm_obj) end
    return CGTperm_obj
end