"""
    CGTQR{S,U,D,L,NZ}

Cached description of the CGT-only basis change used by a QR-style two-factor
split.

`CGTQR` stores a matrix whose rows are the concatenated `(left_om, right_om)`
split basis and whose columns are the original canonical CGT basis rows. Within
each `(q, omL, omR)` block, `left_om` varies fastest and `right_om` varies next,
matching `CGTSVD`.

# Fields
- `qr_arr`: basis-change matrix from canonical CGT coefficients to QR split coordinates.
- `upsp`, `dnsp`: sorted original upper/lower irrep tuples.
- `leftlegs`: normalized selected legs in combined `(upsp..., dnsp...)` indexing.
- `bond_sps`: ordered `(q, omL, omR)` split-basis block descriptors.
- `size_byte`: cached memory footprint for LRU eviction.
"""
struct CGTQR{S<:NonabelianSymm, U, D, L, NZ}
    qr_arr::Array{Float64, 2}
    upsp::NTuple{U, NTuple{NZ, Int}}
    dnsp::NTuple{D, NTuple{NZ, Int}}
    leftlegs::NTuple{L, Int}
    bond_sps::Vector{Tuple{NTuple{NZ, Int}, Int, Int}}
    size_byte::Int

    function CGTQR{S, U, D, L, NZ}(qr_arr, upsp, dnsp, leftlegs, bond_sps,
        size_byte::Int=0) where {S<:NonabelianSymm, U, D, L, NZ}
        if size_byte == 0
            obj = new{S, U, D, L, NZ}(qr_arr, upsp, dnsp, leftlegs, bond_sps, 0)
            size_byte = Base.summarysize(obj)
        end
        new{S, U, D, L, NZ}(qr_arr, upsp, dnsp, leftlegs, bond_sps, size_byte)
    end
end

"""
    getNsave_CGTQR(::Type{S<:AbelianSymm}, upsp, dnsp, leftlegs; save=true, verbose=0)

Return `nothing` for Abelian symmetries, where no nontrivial CGT basis-change
matrix is needed for the QR split.

`upsp` and `dnsp` are the upper and lower q-label tuples of the original CGT.
`leftlegs` identifies the visible legs assigned to the left QR factor.
`save` and `verbose` are accepted for API compatibility with the non-Abelian
method but have no effect for Abelian symmetries.
"""
getNsave_CGTQR(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    leftlegs;
    save=true,
    verbose=0) where {S<:AbelianSymm, U, D, NZ} = nothing

"""
    get_qr_split_qs(::Type{S}, left_up, left_dn, right_up, right_dn) -> Tuple

Find the bond q-labels that can appear in a non-Abelian CGT QR split.

`S` is the symmetry type. `left_up` and `left_dn` are the original upper/lower
spaces assigned to the left factor. `right_up` and `right_dn` are the spaces
assigned to the right factor. The left side treats lower spaces as dual charges,
while the right side treats the split bond `q` as an incoming leg, so the right
candidate set is computed from `right_dn` against dualized `right_up`.

The returned tuple preserves the order from the left-side valid-output list and
contains only q-labels valid on both sides of the split.
"""
function get_qr_split_qs(::Type{S},
    left_up::NTuple{UL, NTuple{NZ, Int}},
    left_dn::NTuple{DL, NTuple{NZ, Int}},
    right_up::NTuple{UR, NTuple{NZ, Int}},
    right_dn::NTuple{DR, NTuple{NZ, Int}}) where {S<:NonabelianSymm, UL, DL, UR, DR, NZ}

    left_merge = stable_sort_tuple((left_up..., map(x -> get_dualq(S, x), left_dn)...))
    # The right factor has q as an incoming leg: (q, right_up) -> right_dn.
    # Therefore q is selected from the residual charge of right_dn against right_up.
    right_merge = stable_sort_tuple((right_dn..., map(x -> get_dualq(S, x), right_up)...))

    right_qs = Set(getNsave_validout(S, right_merge).out_spaces)
    return Tuple(q for q in getNsave_validout(S, left_merge).out_spaces if q in right_qs)
end

"""
    get_qr_split_sector_vectors(::Type{S}, left_up, left_dn, right_up, right_dn, q; verbose=0)

Build the X-symbol block that maps one split-bond sector into the original CGT
outer-multiplicity basis.

`left_up`, `left_dn`, `right_up`, and `right_dn` describe the two QR factors.
`q` is the candidate bond irrep being inserted as an extra lower leg on the
left factor and an extra upper leg on the right factor. `verbose` is forwarded
to `getNsave_Xsymbol`.

The method constructs the augmented left/right spaces, records the positions of
the inserted bond legs, and requests the corresponding X-symbol. It returns the
X-symbol array for this sector, or `nothing` if the sector is invalid or exactly
zero.
"""
function get_qr_split_sector_vectors(::Type{S},
    left_up::NTuple{UL, NTuple{NZ, Int}},
    left_dn::NTuple{DL, NTuple{NZ, Int}},
    right_up::NTuple{UR, NTuple{NZ, Int}},
    right_dn::NTuple{DR, NTuple{NZ, Int}},
    q::NTuple{NZ, Int};
    verbose=0) where {S<:NonabelianSymm, UL, DL, UR, DR, NZ}

    left_dn_q = stable_sort_tuple((left_dn..., q))
    right_up_q = stable_sort_tuple((q, right_up...))

    qpos_left = findlast(==(q), left_dn_q)
    qpos_right = findfirst(==(q), right_up_q)
    @assert !isnothing(qpos_left)
    @assert !isnothing(qpos_right)

    X = getNsave_Xsymbol(S,
        left_up, left_dn_q,
        right_up_q, right_dn,
        (length(left_up) + qpos_left,), (qpos_right,);
        verbose, save=true)
    if isnothing(X) || iszero(X.xsym_arr) return nothing end
    return X.xsym_arr
end

"""
    get_qr_split_basis_matrix(::Type{S}, upsp, dnsp, leftlegs; verbose=0)
        -> (basis_rows, bond_sps)

Construct the full non-Abelian CGT basis-change matrix for a QR split.

`S` is the symmetry. `upsp` and `dnsp` are sorted upper and lower q-label tuples
for the original CGT basis. `leftlegs` is the normalized tuple of visible leg
positions placed on the left QR factor, using combined `(upsp..., dnsp...)`
numbering. `verbose` is forwarded to sector X-symbol construction.

`basis_rows` has one row for each split basis vector and one column for each
original CGT outer-multiplicity basis vector. Rows are grouped by bond sector;
inside a sector the left outer multiplicity varies fastest and the right outer
multiplicity varies next. Each sector block is multiplied by `sqrt(dim(q))`.
`bond_sps` records `(q, omL, omR)` for every nonzero split sector and therefore
defines how rows of `basis_rows` are interpreted.
"""
function get_qr_split_basis_matrix(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    leftlegs::NTuple{L, Int};
    verbose=0) where {S<:NonabelianSymm, U, D, NZ, L}

    @assert issorted(upsp)
    @assert issorted(dnsp)
    @assert 1 <= L < U + D

    total = get_CGTom(S, upsp, dnsp).totalOM

    left_up, left_dn, right_up, right_dn = get_split_side_spaces(S, upsp, dnsp, leftlegs)
    qs = get_qr_split_qs(S, left_up, left_dn, right_up, right_dn)

    basis_rows = zeros(Float64, total, total)
    bond_sps = Tuple{NTuple{NZ, Int}, Int, Int}[]

    row = 1
    for q in qs
        sector_rows = get_qr_split_sector_vectors(S, left_up, left_dn, right_up, right_dn, q; verbose)
        isnothing(sector_rows) && continue

        omL, omR = size(sector_rows, 1), size(sector_rows, 2)
        push!(bond_sps, (q, omL, omR))
        dimq = Float64(dimension(getNsave_irep(S, BigInt, q)))
        sector_mat = reshape(sector_rows, omL * omR, total)
        basis_rows[row:row+omL*omR-1, :] .= sqrt(dimq) .* sector_mat
        row += omL * omR
    end

    @assert row == total + 1
    @assert sum(omL * omR for (_, omL, omR) in bond_sps) == total
    return basis_rows, bond_sps
end

"""
    getNsave_CGTQR(::Type{S<:NonabelianSymm}, upsp, dnsp, leftlegs; save=true, verbose=0)
        -> CGTQR

Load or compute the CGT-only basis data for a symmetry-respecting QR split.

`S` is the non-Abelian symmetry type. `upsp` and `dnsp` are the sorted upper and
lower q-label tuples of the original CGT. `leftlegs` specifies which visible
legs are assigned to the left QR factor before standardization; accepted values
follow the same combined-leg convention as CGTSVD. `save` requests persistence
for cacheable lower-level CGT objects. `verbose` controls diagnostic output in
the X-symbol/basis construction path.

The method standardizes spaces and selected legs, then delegates to
`getNsave_CGTQR_stan`. The returned `CGTQR` maps original CGT outer-multiplicity
coordinates into the split `(q, left_om, right_om)` row basis.
"""
function getNsave_CGTQR(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    leftlegs;
    save=true,
    verbose=0) where {S<:NonabelianSymm, U, D, NZ}

    @assert issorted(upsp) && issorted(dnsp)
    upsp, dnsp, leftlegs = standardize_spaces_and_legs(S, upsp, dnsp, leftlegs, true)
    getNsave_CGTQR_stan(S, upsp, dnsp, leftlegs; save, verbose)
end

"""
    getNsave_CGTQR_stan(::Type{S}, upsp, dnsp, leftlegs; save, verbose) -> Union{Bool,CGTQR}

Compute CGTQR data after the spaces and selected legs have been standardized.

`S` is the non-Abelian symmetry. `upsp` and `dnsp` are the sorted canonical
upper/lower spaces. `leftlegs` is the standardized selected-leg tuple.
`save` is forwarded when permutation CGT data must be loaded or computed.
`verbose` controls diagnostics in split-basis construction.

Empty and full-left splits are represented by `true` and `false`, matching the
sentinel convention used by the CGTSVD split path. Otherwise the method builds
the split-basis matrix, applies the final CGT permutation needed to restore the
canonical original leg order, and returns a `CGTQR` object with cached bond
sector descriptors and memory footprint.
"""
function getNsave_CGTQR_stan(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    leftlegs;
    save,
    verbose) where {S<:NonabelianSymm, U, D, NZ}

    if length(leftlegs) == 0
        return true
    elseif length(leftlegs) == U + D
        return false
    end

    leftlegs_ = normalize_cgtsvd_leftlegs(leftlegs, U + D)
    qr_arr, bond_sps = get_qr_split_basis_matrix(S, upsp, dnsp, leftlegs_; verbose)
    final_perm = get_cgtsvd_final_perm(S, upsp, dnsp, leftlegs_)
    if final_perm != Tuple(1:(U + D))
        cgtperm = getNsave_CGTperm(S, upsp, dnsp, final_perm; save=true)
        @assert !isnothing(cgtperm)
        qr_arr = qr_arr * cgtperm.perm_arr
    end

    return CGTQR{S, U, D, length(leftlegs_), NZ}(qr_arr, upsp, dnsp, leftlegs_, bond_sps)
end
