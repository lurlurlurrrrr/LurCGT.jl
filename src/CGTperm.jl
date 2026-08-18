"""
    CGTperm{S,U,D,N,NZ}

Cached outer-multiplicity basis transform for a non-Abelian CGT leg
permutation.

`S` is the symmetry type. `U` and `D` are the numbers of upper and lower CGT
legs. `N == U + D` is the number of permuted physical legs. `NZ` is the q-label
tuple width for `S`.

The transform is only needed when a permutation changes the canonical
fusion-tree basis inside a non-Abelian CGT. Abelian symmetries and permutations
that become identity after removing trivial zero legs use `nothing` instead.

# Fields

- `perm_arr`: dense matrix whose columns are source flattened OM basis states
  and whose rows are permuted canonical OM basis states.
- `upsp`: sorted upper/outgoing q-label tuple used to define the source CGT.
- `dnsp`: sorted lower/incoming q-label tuple used to define the source CGT.
- `perm`: requested physical-leg permutation in concatenated
  `(upsp..., dnsp...)` indexing. Standardized cached objects only permute upper
  legs among upper positions and lower legs among lower positions.
- `size_byte`: cached memory footprint for LRU eviction.
"""
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

"""
    perm_FTrees!(FTrees_dict, perm) -> nothing

Permute every fusion tree stored in a CGT canonical-basis dictionary.

`FTrees_dict[csp]` is the vector of `FTree` basis elements for intermediate
space `csp`. `perm` is an incoming-leg permutation for those trees, expressed
in the output-position convention consumed by `Base.permute!(::FTree, ...)`.
The function mutates each vector in place by replacing every basis tree with
the result of the exact R/F-symbol recoupling permutation.
"""
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

"""
    fill_CGTperm_matrix!(S, CGTperm_arr, CGT_oms, CGT_FTrees_up, CGT_FTrees_dn)
        -> nothing

Fill the OM-basis transform matrix for a CGT leg permutation.

`S` is the symmetry type. `CGTperm_arr` is the preallocated dense matrix to
fill; its shape must be `(CGT_oms.totalOM, CGT_oms.totalOM)`. `CGT_oms`
describes the flattened source outer-multiplicity ordering. `CGT_FTrees_up` and
`CGT_FTrees_dn` contain the already-permuted upper and lower fusion-tree basis
vectors, keyed by central q-label.

For source column `i`, the method looks up `(central space, upper OM index,
lower OM index)`, contracts the corresponding upper/down trees through a unit
identity tree on the central space, converts the result back to the canonical
OM vector, and writes that vector into column `i`.
"""
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

"""
    remove_zeros(::Type{S}, spaces, perm) -> (spaces_, perm_)

Remove leading trivial qlabel legs from a CGT permutation problem.

`S` supplies the zero q-label width. `spaces` is a tuple of q-labels from either
the upper or lower side of a CGT. `perm` is the side-local permutation over
those positions. Leading zero-q-label legs do not contribute nontrivial
non-Abelian permutation data, so the method strips the prefix ending at the last
leading trivial label and shifts the remaining permutation down by the removed
count.

If every space is trivial, the standardized representation is a single trivial
space with identity permutation. This keeps downstream code from having to
handle empty upper/lower fusion-tree problems.
"""
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

"""
    getNsave_CGTperm(::Type{S}, upsp, dnsp, perm; save=true) -> Union{Nothing,CGTperm}

Load or compute the OM-basis transform for a CGT leg permutation.

`S` is the symmetry type. `upsp` and `dnsp` are sorted q-label tuples for the
upper/outgoing and lower/incoming CGT leg groups. `perm` is the requested
permutation in concatenated `(upsp..., dnsp...)` indexing. `save` controls
SQLite persistence for generated non-Abelian transforms.

For Abelian symmetries the method returns `nothing`, because no nontrivial
outer-multiplicity basis transform exists. For non-Abelian symmetries, leading
zero q-labels are removed independently on the upper and lower sides before the
standard cache key is formed. If the standardized permutation is identity,
`nothing` is returned as an effective no-op.
"""
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

"""
    getNsave_CGTperm_std(::Type{S}, upsp, dnsp, perm; save=true) -> CGTperm

Load or compute a standardized nontrivial non-Abelian CGT permutation.

`S` is the symmetry type. `upsp` and `dnsp` are standardized sorted q-label
tuples with removable trivial leading labels already stripped. `perm` is the
standardized concatenated permutation and must preserve the upper/lower split:
upper output positions map to upper source legs and lower output positions map
to lower source legs. `save` controls persistence on a cache miss.

The function validates that applying `perm` leaves the q-label multiset in the
same order required by canonical CGT storage, then checks SQLite cache before
delegating to `computeNsave_CGTperm_std`.
"""
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
    return computeNsave_CGTperm_std(S, upsp, dnsp, perm; save)
end

"""
    computeNsave_CGTperm_std(::Type{S}, upsp, dnsp, perm; save=true) -> CGTperm

Compute a non-Abelian CGT permutation matrix from fusion-tree recouplings.

`S` is the symmetry type. `upsp` and `dnsp` define the standardized canonical
CGT fusion problem. `perm` permutes legs within the upper and lower groups.
`save` determines whether the completed object is written to SQLite.

The method builds canonical upper/down fusion-tree bases for every common
central space, applies side-local tree permutations, contracts each permuted
basis pair back through a central identity tree, and fills the dense transform
matrix column by column. When the total outer multiplicity is greater than one,
CGT norms are used to convert between normalized and unnormalized OM vector
conventions before optional persistence.
"""
function computeNsave_CGTperm_std(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    perm::NTuple{N, Int};
    save=true) where {S<:NonabelianSymm, U, D, NZ, N}

    perm_up = perm[1:U]; perm_dn = Tuple(i-U for i in perm[U+1:end])
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
