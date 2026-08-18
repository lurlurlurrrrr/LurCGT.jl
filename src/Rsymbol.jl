# Overall structure of this file is similar to Fsymbol.jl

"""
R-move (exchange) coefficients for a pair of CGT fusion spaces.

# Fields

- `in`, `out`: q-labels defining the exchanged fusion channel.
- `rsym_mat`: dense R-move matrix in the flattened OM basis.
- `nfactor`: exact normalization factors for matrix entries/basis states.
- `size_byte`: cached memory footprint for LRU eviction.
"""
struct Rsymbol{S<:NonabelianSymm, CT, NZ}
    # Only need degenerate case, so only one q-label is needed
    in::NTuple{NZ, Int}
    out::NTuple{NZ, Int}
    rsym_mat::Array{CT, 2}
    nfactor::Vector{Rational{CT}}
    # Cached memory footprint in bytes, computed once at construction for LRU cache eviction by size
    size_byte::Int
    
    function Rsymbol{S, CT, NZ}(in, out, rsym_mat, nfactor, size_byte::Int=0) where {S, CT, NZ}
        if size_byte == 0
            obj = new{S, CT, NZ}(in, out, rsym_mat, nfactor, 0)
            size_byte = Base.summarysize(obj)
        end
        new{S, CT, NZ}(in, out, rsym_mat, nfactor, size_byte)
    end
end

"""
    getNsave_Rsymbol(::Type{S}, ::Type{CT}, in, out; verbose=0) -> Rsymbol

Load or compute the exchange matrix for the degenerate fusion `in x in -> out`.

`S` selects the symmetry, `CT` selects the stored scalar type, `in` is the
exchanged qlabel, and `out` is the fusion output qlabel. The method first checks
SQLite/cache storage. On a miss, it loads the relevant CG3 block, keeps the
maximal-weight output part, exchanges the first two inputs, applies the product
metric to the original block, contracts the overlap matrix, saves the result,
and returns an `Rsymbol`.
"""
function getNsave_Rsymbol(::Type{S},
    ::Type{CT},
    in::NTuple{NZ, Int},
    out::NTuple{NZ, Int};
    verbose=0) where {S<:NonabelianSymm, CT<:Number, NZ}

    @assert NZ == nzops(S)
    # If output is trivial representation, get it from 1jsym function
    if out == Tuple(0 for _ in 1:NZ)
        if verbose > 0 println("R-symbol from 1j-symbol") end
        getNsave_1jsym(S, BigInt, CT, in)
    end

    loaded = load_Rsymbol_sqlite(S, CT, in, out)
    if !isnothing(loaded) return loaded end

    blk, fac = load_cg3blk(S, CT, (in, in), [out])[out]
    blk_mw = mwpartof(S, blk, (in, in, out), 3)
    om = length(fac)
    rsym_mat = zeros(Float64, om, om)

    blk_permuted = permute12(blk_mw)
    conjugate_rsym!(S, blk_mw, in)
    for key in keys(blk_mw)
        if haskey(blk_permuted, key)
            @tensor mat2add[ν, μ] := 
            blk_permuted[key][in1, in2, μ] * blk_mw[key][in1, in2, ν]
            rsym_mat += mat2add
        end
    end
    rsym_struct = Rsymbol{S, CT, NZ}(in, out, rsym_mat, fac)
    save_Rsymbol_sqlite(rsym_struct)
    return rsym_struct
end

"""
    permute12(blk) -> Dict

Swap the first two input axes of a CG3 block dictionary.

`blk[(z1, z2)]` stores a rank-3 array with axes `(input1, input2, om)`. The
returned dictionary uses keys `(z2, z1)` and arrays permuted as `(2, 1, 3)`.
"""
function permute12(blk::Dict{NTuple{2, NTuple{NZ, Int}}, Array{CT, 3}}) where {CT<:Number, NZ}
    permuted_blk = Dict{NTuple{2, NTuple{NZ, Int}}, Array{CT, 3}}()
    for (key, value) in blk
        new_key = (key[2], key[1])
        permuted_blk[new_key] = permutedims(value, (2, 1, 3))
    end
    return permuted_blk
end

"""
    conjugate_rsym!(::Type{S}, blk, in)

Apply product inner-product metrics to an R-symbol CG3 block in place.

`blk` maps `(z1, z2)` weight pairs to rank-3 arrays. `in` selects the irrep used
for both inputs. For each block, the two input axes are contracted with
`irep.innerprod[z1]` and `irep.innerprod[z2]`.
"""
function conjugate_rsym!(::Type{S},
    blk::Dict{NTuple{2, NTuple{NZ, Int}}, Array{CT, 3}},
    in::NTuple{NZ, Int}) where {S<:NonabelianSymm, CT<:Number, NZ}

    irep = getNsave_irep(S, CT, in)
    for key in keys(blk)
        b = blk[key]; @assert ndims(b) == 3
        z1, z2 = key
        m1 = irep.innerprod[z1]
        m2 = irep.innerprod[z2]
        @tensor newblk[j1, j2, o] :=
            b[i1, i2, o] * m1[i1, j1] * m2[i2, j2]
        blk[key] = newblk
    end
end

"""
    rsym_rightnormalized(rsym::Rsymbol) -> (Matrix{BigInt}, Integer)

Convert an R-symbol to an integer right-normalized matrix.

`rsym.rsym_mat` is multiplied by `Diagonal(rsym.nfactor)`. The least common
multiple of denominators is pulled out as `nfac`, and the returned matrix is
integer-valued. This form is used by exact fusion-tree recoupling.
"""
function rsym_rightnormalized(rsym::Rsymbol{S, CT, NZ}) where {S<:NonabelianSymm, CT<:Number, NZ}
    @assert nzops(S) == NZ
    mat_rational = rsym.rsym_mat * Diagonal(rsym.nfactor)
    nfac = lcm(denominator.(mat_rational))
    return Matrix{BigInt}(mat_rational * nfac), nfac
end
