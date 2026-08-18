# Contraction of CG3 and 1j-symbol

# Always the second incoming space of original CG3 is contracted to 1j-symbol
"""
Basis transform obtained by contracting one CG3 leg with a 1j tensor.

# Fields

- `incom_spaces`: two original CG3 incoming q-labels.
- `out_space`: original CG3 output q-label.
- `flip_mat`: dense map between original and flipped outer-multiplicity bases.
- `nfactor`: exact normalization factors for the map.
- `size_byte`: cached memory footprint for LRU eviction.
"""
struct CG3Flip{S<:NonabelianSymm, CT, NZ}
    # Incoming spaces of original CG3
    incom_spaces::NTuple{2, NTuple{NZ, Int}} 
    # Outgoing space of original CG3
    out_space::NTuple{NZ, Int}
    # The matrix encodes the contraction of CG3 and 1j-symbol
    flip_mat::Array{CT, 2}
    orig_nfactor::Vector{Rational{CT}}
    new_nfactor::Vector{Rational{CT}}

    size_byte::Int

    function CG3Flip{S, CT, NZ}(incom_spaces, out_space, flip_mat, 
        orig_nfactor, new_nfactor, size_byte::Int=0) where {S<:NonabelianSymm, CT<:Number, NZ}
        if size_byte == 0
            obj = new{S, CT, NZ}(incom_spaces, out_space, flip_mat, 
            orig_nfactor, new_nfactor, 0)
            size_byte = Base.summarysize(obj)
        end
        new{S, CT, NZ}(incom_spaces, out_space, flip_mat, 
        orig_nfactor, new_nfactor, size_byte)
    end
end

"""
    getNsave_cg3flip(::Type{S}, ::Type{RT}, ::Type{CT}, incom_spaces, out_space; verbose=0) -> CG3Flip

Load or compute the basis transform obtained by flipping one CG3 leg through a
1j tensor.

`S` selects the symmetry, `RT` is the exact irrep/1j coefficient type, and `CT`
is the stored CGT scalar type. `incom_spaces == (in1, in2)` are the original
CG3 input qlabels and `out_space` is the original output qlabel. On a cache
miss, the method loads the original CG3, loads and metric-conjugates the 1j
block for `in2`, loads the target CG3 `(out_space, dual(in2)) -> in1`,
contracts maximal-weight blocks, stores original and new normalization factors,
saves the transform, and returns it.
"""
# Obtain the CG3Flip object by contracting CG3 with 1j-symbol
function getNsave_cg3flip(::Type{S},
    ::Type{RT},
    ::Type{CT},
    incom_spaces::NTuple{2, NTuple{NZ, Int}},
    out_space::NTuple{NZ, Int};
    verbose=0) where {S<:NonabelianSymm, NZ, RT<:Number, CT<:Number}

    # Try to load from file, use it if exists
    @assert NZ == nzops(S)
    loaded = load_cg3flip_sqlite(S, CT, incom_spaces, out_space)
    if !isnothing(loaded) return loaded end

    # If not, generate it 
    in1, in2 = incom_spaces
    # Load original CG3
    oblk, ofac = load_cg3blk(S, CT, (in1, in2), [out_space])[out_space]
    dualin2 = get_dualq(S, in2)
    blk1j, fac1j = load_1jblk(S, RT, CT, in2)
    conjugate_1j!(S, blk1j, fac1j, in2, dualin2)
    nblk, nfac = load_cg3blk(S, CT, (out_space, dualin2), [in1])[in1]

    oblk_mw = mwpartof(S, oblk, (in1, in2, out_space), 1)
    nblk_mw = mwpartof(S, nblk, (out_space, dualin2, in1), 3)

    @assert length(ofac) == length(nfac) "Outer multiplicities do not match"
    om = length(ofac)
    flip_arr = contract_blks(S, oblk_mw, blk1j, nblk_mw, om; verbose)

    # Construct and save the generated struct
    cg3flip_struct = CG3Flip{S, CT, NZ}(incom_spaces, out_space, 
        flip_arr, ofac .* fac1j, nfac)
    save_cg3flip_sqlite(S, cg3flip_struct)
    return cg3flip_struct
end

"""
    cg3flip_rightnormalized(cg3flip::CG3Flip) -> (Matrix{BigInt}, Integer)

Convert a CG3 flip transform to an integer right-normalized matrix.

`cg3flip.flip_mat` is multiplied by `Diagonal(cg3flip.new_nfactor)`. The least
common multiple of denominators is returned as `nfac`, and the matrix is scaled
to integer entries for exact recoupling routines.
"""
function cg3flip_rightnormalized(cg3flip::CG3Flip{S, CT, NZ}) where {S<:NonabelianSymm, CT<:Number, NZ}
    @assert nzops(S) == NZ
    mat_rational = cg3flip.flip_mat * Diagonal(cg3flip.new_nfactor)
    nfac::Integer = lcm(denominator.(mat_rational))
    return Matrix{BigInt}(mat_rational * nfac), nfac
end

"""
    contract_blks(::Type{S}, oblk, blk1j, nblk, om; verbose=0)

Contract original CG3, 1j, and flipped CG3 maximal-weight blocks.

`oblk` is the original CG3 block dictionary, `blk1j` is the conjugated 1j block,
`nblk` is the new CG3 block dictionary, and `om` is the shared outer
multiplicity. The result is the dense `om x om` flip matrix.
"""
function contract_blks(::Type{S},
    oblk::Dict{NTuple{2, NTuple{NZ, Int}}, Array{CT, 3}},
    blk1j::Dict{NTuple{2, NTuple{NZ, Int}}, Array{CT}},
    nblk::Dict{NTuple{2, NTuple{NZ, Int}}, Array{CT, 3}},
    om::Int; verbose=0) where {S<:NonabelianSymm, CT<:Number, NZ}

    flip_arr = zeros(CT, om, om)

    for ((ok1, ok2), oarr) in oblk
        arr1j = blk1j[(ok1, .-ok1)]
        narr = nblk[(ok2, .-ok1)]

        @tensor cont_res[oom, om1j, nom] := 
        oarr[i1, i2, oom] * arr1j[i1, i3, om1j] * narr[i2, i3, nom]
        flip_arr += reshape(cont_res, om, om)
    end

    return flip_arr
end

"""
    conjugate_1j!(::Type{S}, blk, fac, in1, in2)

Apply product inner-product metrics to a 1j block dictionary in place.

`blk[(z1, z2)]` is contracted on both input axes with the inner-product matrices
of irreps `in1` and `in2`. `fac` must contain the single 1j normalization
factor; it is checked for consistency but not modified.
"""
function conjugate_1j!(::Type{S},
    blk::Dict{NTuple{2, NTuple{NZ, Int}}, Array{CT}},
    fac::Vector{Rational{CT}},
    in1::NTuple{NZ, Int},
    in2::NTuple{NZ, Int}) where {S<:NonabelianSymm, CT<:Number, NZ}

    @assert length(fac) == 1; facden = fac[1].den
    irep1 = getNsave_irep(S, CT, in1)
    irep2 = getNsave_irep(S, CT, in2)
    for k in keys(blk)
        z1, z2 = k
        m1 = irep1.innerprod[z1]
        m2 = irep2.innerprod[z2]
        @tensor newblk[j1, j2, o] := blk[k][i1, i2, o] * m1[i1, j1] * m2[i2, j2]
        blk[k] = newblk
    end
end
