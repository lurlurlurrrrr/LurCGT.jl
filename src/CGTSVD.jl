"""
    CGTSVD{S,U,D,L,NZ}

Cached description of the CGT-only basis change used by the QSpace-side SVD
construction.

`CGTSVD` stores the unitary matrix that maps coefficients in the canonical CGT
basis of `(upsp, dnsp)` to the orthonormal split basis defined by `leftlegs`.
The split basis is grouped by intermediate bond sector `q`, with multiplicities
recorded in `bond_sps`.

# Fields
- `svd_arr`: unitary matrix taking canonical CGT coefficients to split-basis
  coefficients.
- `upsp`: sorted upper-leg irrep tuple of the original CGT.
- `dnsp`: sorted lower-leg irrep tuple of the original CGT.
- `leftlegs`: normalized leg selection in the combined ordering
  `(upsp..., dnsp...)`.
- `bond_sps`: ordered `(q, omL, omR)` blocks of the split basis.
- `size_byte`: cached object size used by the package cache machinery.
"""
struct CGTSVD{S<:NonabelianSymm, U, D, L, NZ}
    svd_arr::Array{Float64, 2}
    upsp::NTuple{U, NTuple{NZ, Int}}
    dnsp::NTuple{D, NTuple{NZ, Int}}
    leftlegs::NTuple{L, Int}
    bond_sps::Vector{Tuple{NTuple{NZ, Int}, Int, Int}}
    size_byte::Int

    function CGTSVD{S, U, D, L, NZ}(svd_arr, upsp, dnsp, leftlegs, bond_sps,
        size_byte::Int=0) where {S<:NonabelianSymm, U, D, L, NZ}
        if size_byte == 0
            obj = new{S, U, D, L, NZ}(svd_arr, upsp, dnsp, leftlegs, bond_sps, 0)
            size_byte = Base.summarysize(obj)
        end
        new{S, U, D, L, NZ}(svd_arr, upsp, dnsp, leftlegs, bond_sps, size_byte)
    end
end

"""
    getNsave_CGTSVD(S, upsp, dnsp, leftlegs; save=true, verbose=0)

Construct the CGT-side basis change associated with splitting a CGT into a
left part and a right part.

The input CGT is specified by its sorted upper-leg spaces `upsp` and sorted
lower-leg spaces `dnsp`. `leftlegs` is a list of leg indices in the combined
ordering `(upsp..., dnsp...)`; the selected legs are assigned to the left side
and the remaining legs to the right side.

For every admissible intermediate sector `q`, this routine builds the split
basis

`(left CGT with outgoing q) ⊗ (all-incoming 1j bridge with q and dual(q)) ⊗ (right CGT with outgoing dual(q))`

and re-expands that three-CGT contraction in the canonical CGT basis of the
original object. The returned matrix is the unitary transformation from
canonical CGT coefficients to this new split basis:

`coeff_split = obj.svd_arr * coeff_canonical`

The rows of `obj.svd_arr` are ordered first by increasing intermediate sector
`q`, then by the left/right outer-multiplicity ordering used internally by
`LurCGT`. The field `obj.bond_sps` stores this block structure as
`(q, omL, omR)` triples.

Returns `nothing` for Abelian symmetries.

# Arguments
- `S`: non-abelian symmetry type such as `SU{2}` or `SU{3}`.
- `upsp`: sorted tuple of upper-leg irreps of the CGT.
- `dnsp`: sorted tuple of lower-leg irreps of the CGT.
- `leftlegs`: non-empty proper subset of leg indices in `(upsp..., dnsp...)`.

# Keyword arguments
- `save=true`: cache the resulting `CGTSVD` object in the SQLite store.
- `verbose=0`: forwarded to internal helper routines when building the basis.

# Returns
- `CGTSVD`: object containing the unitary matrix `svd_arr`, the split
  specification, and the resulting bond-sector list `bond_sps`.
"""
getNsave_CGTSVD(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    leftlegs;
    save=true,
    verbose=0) where {S<:AbelianSymm, U, D, NZ} = nothing

# Validate and normalize the user-provided left-leg selection for `CGTSVD`.
# The result is a sorted tuple of unique leg indices in the combined ordering
# `(upsp..., dnsp...)`, and must be a non-empty proper subset of `1:total`.
function normalize_cgtsvd_leftlegs(leftlegs, total::Int)
    legs = sort!(collect(Int, leftlegs))
    length(legs) > 0 || throw(ArgumentError("CGTSVD requires at least one left leg"))
    length(legs) < total || throw(ArgumentError(
        "CGTSVD requires at most $(total) left legs for a rank-$total CGT"))
    all(1 .<= legs .<= total) || throw(ArgumentError(
        "CGTSVD left legs must lie between 1 and $total"))
    length(unique(legs)) == length(legs) || throw(ArgumentError(
        "CGTSVD left legs must not contain duplicates"))
    return Tuple(legs)
end

# Stable sorting is needed so an added bond leg keeps a deterministic slot even
# when its q-label already appears among the original legs.
stable_sort_tuple(spaces) = Tuple(sort!(collect(spaces); alg=MergeSort))

get_zeroq(::Type{S}) where {S<:Symmetry} = Tuple(0 for _ in 1:nzops(S))

function get_leftfirst_blockperm(spaces, leftmask)
    perm = Int[]
    i = 1
    while i <= length(spaces)
        j = i
        while j <= length(spaces) && spaces[j] == spaces[i]
            j += 1
        end
        for k in i:j-1
            if leftmask[k]
                push!(perm, k)
            end
        end
        for k in i:j-1
            if !leftmask[k]
                push!(perm, k)
            end
        end
        i = j
    end
    return Tuple(perm)
end

function get_cgtsvd_final_perm(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    leftlegs::NTuple{L, Int}) where {S<:NonabelianSymm, U, D, NZ, L}

    leftset = Set(leftlegs)
    leftmask_up = ntuple(i -> i in leftset, U)
    leftmask_dn = ntuple(i -> (U + i) in leftset, D)
    perm_up = get_leftfirst_blockperm(upsp, leftmask_up)
    perm_dn = get_leftfirst_blockperm(dnsp, leftmask_dn)
    return (perm_up..., (U .+ perm_dn)...)
end

# Split the original CGT legs into left and right parts while preserving the
# original incoming/outgoing directions of every physical leg.
function get_split_side_spaces(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    leftlegs::NTuple{L, Int}) where {S<:NonabelianSymm, U, D, NZ, L}

    leftset = Set(leftlegs)
    left_up = Tuple(upsp[i] for i in 1:U if i in leftset)
    left_dn = Tuple(dnsp[i] for i in 1:D if (U + i) in leftset)
    right_up = Tuple(upsp[i] for i in 1:U if !(i in leftset))
    right_dn = Tuple(dnsp[i] for i in 1:D if !((U + i) in leftset))
    return left_up, left_dn, right_up, right_dn
end

# Find the admissible intermediate sectors. Here the physical outgoing legs are
# dualized only for the sector-selection step; the actual basis construction
# below keeps the original leg directions intact.
function get_split_qs(::Type{S},
    left_up::NTuple{UL, NTuple{NZ, Int}},
    left_dn::NTuple{DL, NTuple{NZ, Int}},
    right_up::NTuple{UR, NTuple{NZ, Int}},
    right_dn::NTuple{DR, NTuple{NZ, Int}}) where {S<:NonabelianSymm, UL, DL, UR, DR, NZ}

    left_merge = stable_sort_tuple((left_up..., map(x -> get_dualq(S, x), left_dn)...))
    right_merge = stable_sort_tuple((right_up..., map(x -> get_dualq(S, x), right_dn)...))

    right_outs = Set(getNsave_validout(S, right_merge).out_spaces)
    return Tuple(q for q in getNsave_validout(S, left_merge).out_spaces if get_dualq(S, q) in right_outs)
end

# Construct the canonical-basis coefficient vectors for one fixed intermediate
# bond sector `q` using two X-symbol contractions:
# 1. contract the left CGT with the all-incoming 1j bridge,
# 2. contract the resulting intermediate CGT with the right CGT.
function get_split_sector_vectors(::Type{S},
    left_up::NTuple{UL, NTuple{NZ, Int}},
    left_dn::NTuple{DL, NTuple{NZ, Int}},
    right_up::NTuple{UR, NTuple{NZ, Int}},
    right_dn::NTuple{DR, NTuple{NZ, Int}},
    q::NTuple{NZ, Int},
    zeroq::NTuple{NZ, Int};
    verbose=0) where {S<:NonabelianSymm, UL, DL, UR, DR, NZ}

    dualq = get_dualq(S, q)
    center_up = stable_sort_tuple((q, dualq))
    center_dn = (zeroq,)

    left_dn_q = stable_sort_tuple((left_dn..., q))
    right_dn_dualq = stable_sort_tuple((dualq, right_dn...))
    interm_up = stable_sort_tuple((dualq, left_up...))
    interm_dn = left_dn

    qpos_center = findfirst(==(q), center_up)
    qpos_left = findlast(==(q), left_dn_q)
    dualqpos_interm = findfirst(==(dualq), interm_up)
    dualqpos_right = findfirst(==(dualq), right_dn_dualq)
    @assert !isnothing(qpos_center)
    @assert !isnothing(qpos_left)
    @assert !isnothing(dualqpos_interm)
    @assert !isnothing(dualqpos_right)

    X1 = getNsave_Xsymbol(S,
        center_up, center_dn,
        left_up, left_dn_q,
        (qpos_center,), (length(left_up) + qpos_left,);
        verbose, save=true)
    if isnothing(X1) || iszero(X1.xsym_arr) return nothing end

    X2 = getNsave_Xsymbol(S,
        interm_up, interm_dn,
        right_up, right_dn_dualq,
        (dualqpos_interm,), (length(right_up) + dualqpos_right,);
        verbose, save=true)
    if isnothing(X2) || iszero(X2.xsym_arr) return nothing end

    @assert size(X1.xsym_arr, 1) == 1
    omL = size(X1.xsym_arr, 2)
    omC = size(X1.xsym_arr, 3)
    omR = size(X2.xsym_arr, 2)
    total = size(X2.xsym_arr, 3)
    @assert size(X2.xsym_arr, 1) == omC

    return reshape(
        reshape(X1.xsym_arr, omL, omC) * reshape(X2.xsym_arr, omC, :),
        omL, omR, total)
end

# Assemble the orthonormal split basis directly in canonical CGT coordinates.
# The companion `bond_sps` stores each q-sector as `(q, omL, omR)`.
function get_split_basis_matrix(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    leftlegs::NTuple{L, Int};
    verbose=0) where {S<:NonabelianSymm, U, D, NZ, L}

    @assert issorted(upsp)
    @assert issorted(dnsp)
    @assert 1 <= L < U + D

    zeroq = get_zeroq(S)
    total = get_CGTom(S, upsp, dnsp).totalOM

    left_up, left_dn, right_up, right_dn = get_split_side_spaces(S, upsp, dnsp, leftlegs)
    qs = get_split_qs(S, left_up, left_dn, right_up, right_dn)

    # Rows of the final change-of-basis matrix in canonical CGT coordinates.
    basis_rows = zeros(Float64, total, total)
    # Sector bookkeeping: `(q, omL, omR)`.
    bond_sps = Tuple{NTuple{NZ, Int}, Int, Int}[]

    row = 1
    for q in qs
        sector_rows = get_split_sector_vectors(S, left_up, left_dn, right_up, right_dn, q, zeroq; verbose)
        isnothing(sector_rows) && continue

        omL, omR = size(sector_rows, 1), size(sector_rows, 2)
        push!(bond_sps, (q, omL, omR))
        dimq = Float64(dimension(getNsave_irep(S, BigInt, q)))
        sector_mat = reshape(sector_rows, omL * omR, total)
        basis_rows[row:row+omL*omR-1, :] .= dimq .* sector_mat
        row += omL * omR
    end

    @assert row == total + 1
    @assert sum(omL * omR for (_, omL, omR) in bond_sps) == total
    return basis_rows, bond_sps
end

function getNsave_CGTSVD(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    leftlegs;
    save=true,
    verbose=0) where {S<:NonabelianSymm, U, D, NZ}

    @assert issorted(upsp) && issorted(dnsp)
    upsp, dnsp, leftlegs = standardize_spaces_and_legs(S, upsp, dnsp, leftlegs, true)
    getNsave_CGTSVD_stan(S, upsp, dnsp, leftlegs; save, verbose)
end

function getNsave_CGTSVD_stan(::Type{S}, 
    upsp::NTuple{U, NTuple{NZ, Int}}, 
    dnsp::NTuple{D, NTuple{NZ, Int}}, 
    leftlegs; 
    save, 
    verbose) where {S<:NonabelianSymm, U, D, NZ}

    if length(leftlegs) == 0 return true 
    elseif length(leftlegs) == U + D return false end
    
    leftlegs_ = normalize_cgtsvd_leftlegs(leftlegs, U + D)
    loaded = load_CGTSVD_sqlite(S, upsp, dnsp, leftlegs_)
    if !isnothing(loaded) return loaded end

    svd_arr, bond_sps = get_split_basis_matrix(S, upsp, dnsp, leftlegs_; verbose)
    final_perm = get_cgtsvd_final_perm(S, upsp, dnsp, leftlegs_)
    if final_perm != Tuple(1:(U + D))
        cgtperm = getNsave_CGTperm(S, upsp, dnsp, final_perm; save=true)
        @assert !isnothing(cgtperm)
        svd_arr = svd_arr * cgtperm.perm_arr
    end

    obj = CGTSVD{S, U, D, length(leftlegs_), NZ}(svd_arr, upsp, dnsp, leftlegs_, bond_sps)
    if save save_CGTSVD_sqlite(obj) end
    return obj
end
