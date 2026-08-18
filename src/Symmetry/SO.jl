# For SO(N), the symmetries are greatly different for even and odd N
# Start from odd N first
"""
    isvalidsymm(::Type{SO{N}}) -> Bool

Return whether `SO{N}` is supported by this implementation.

`N` is the defining dimension. The code supports `SO{N}` for `N >= 4`; both
odd and even cases are handled, but several qlabel and tableau conventions
depend on the parity of `N`.
"""
isvalidsymm(::Type{SO{N}}) where N = 4 <= N

isSON(::Type{<:SO}) = true
isSON(::Any) = false

# dimension of defining irep
defirepdim(::Type{SO{N}}) where N = N

# Text expression of the symmetry.
totxt(::Type{SO{N}}) where N = "SO$(N)"

# number of raising operators (for SO(N), it is floor(N / 2))
nlops(::Type{SO{N}}) where N = div(N, 2)

# number of z-operators (for SO(N), it is floor(N / 2))
nzops(::Type{SO{N}}) where N = div(N, 2)

maxalt(::Type{SO{N}}) where N = div(N, 2)

getsl_triv(::Type{SO{N}}, ::Type{RT}) where {N, RT<:Number} =
Tuple(Dict{NTuple{div(N, 2), Int}, SparseMatrixCSC{RT}}() for _=1:div(N, 2))

"""
    getsz_triv(::Type{SO{N}}) -> Dict

Return the trivial SO(N) weight-sector dictionary.

The only sector is the all-zero tuple of length `div(N, 2)`, and its basis
index range is `(1, 1)`.
"""
getsz_triv(::Type{SO{N}}) where N = Dict(Tuple(0 for _=1:div(N, 2))=>(1, 1))

"""
    getsl_def(::Type{SO{N}}, ::Type{RT}) -> NTuple{div(N,2), Dict}

Build sparse lowering-operator blocks for the defining SO(N) representation.

`RT` is the scalar type of the sparse matrices. The returned tuple has one
dictionary per simple lowering operator. Odd and even `N` differ on the last
simple root, so the source z-weight sectors are chosen with the parity-aware
SO(N) convention.
"""
function getsl_def(::Type{SO{N}}, ::Type{RT}) where {N, RT<:Number} 
    NZ = div(N, 2); Nodd = (N % 2 == 1)
    def_sl = Tuple(Dict{NTuple{div(N, 2), Int}, SparseMatrixCSC{RT}}() for _=1:div(N, 2))
    def_sz = getsz_def(SO{N})
    def_zvals = sorted_zvals(def_sz)
    for i=1:NZ
        dz = getdz(SO{N}, i)
        val = spzeros(RT, 1, 1); val[1, 1] = RT(1)
        sw = def_zvals[!Nodd&&i==NZ ? i-1 : i]
        def_sl[i][sw] = val
        def_sl[i][.-(sw.-dz)] = val
    end
    return def_sl
end

"""
    getsz_def_vec_(::Type{SO{N}}) -> Vector{NTuple{div(N,2), Int}}

Return defining-representation weights in the internal doubled convention.

SO(N) spin representations require half-integral weights. This code stores all
SO(N) z-weights multiplied by two, so defining-vector weights are `+/-2` on one
coordinate, with an additional zero weight for odd `N`.
"""
# We need to consider spin representation, so weights are multiplied by 2
function getsz_def_vec_(::Type{SO{N}}) where N 
    szs = Vector{NTuple{div(N, 2), Int}}()
    for i in 1:N
        sz = zeros(Int, div(N, 2))
        if N % 2 == 1 && i == N push!(szs, Tuple(sz)); continue; end
        dv, rem = divrem(i+1, 2)
        sz[dv] = rem == 0 ? 2 : -2
        push!(szs, Tuple(sz))
    end
    sort!(szs; lt=rev_less, rev=true)
    return szs
end

getsz_def_vec(::Type{SO{N}}) where N = [collect(t) for t in getsz_def_vec_(SO{N})]

getsz_def(::Type{SO{N}}) where N = Dict(w => (i, i) for (i, w) in enumerate(getsz_def_vec_(SO{N})))


"""
    getdz(::Type{SO{N}}, lop::Int) -> NTuple{div(N,2), Int}

Return the z-weight change produced by SO(N) lowering operator `lop`.

`lop` is a simple-root index. For even SO(N), the last lowering operator uses a
different pair of defining weights than the previous simple roots; this helper
encodes that parity-specific convention.
"""
function getdz(::Type{SO{N}}, lop::Int) where N
    NZ = div(N, 2); Nodd = (N % 2 == 1)
    sz_lst = getsz_def_vec_(SO{N})
    if Nodd || lop < NZ return sz_lst[lop] .- sz_lst[lop+1] end
    return sz_lst[lop-1] .- sz_lst[lop+1]
end

"""
    crystal_chars_map(::Type{SO{N}}) -> Dict{Int, Int}

Map SO(N) tableau characters to defining-representation basis positions.

Positive characters represent the first half of the vector weights, negative
characters represent the dual half, and odd SO(N) includes the extra zero
character.
"""
function crystal_chars_map(::Type{SO{N}}) where N
    charmap = Dict{Int, Int}(); 
    NZ = div(N, 2); Nodd = N % 2
    for i in 1:NZ charmap[i] = i end
    for i in NZ+1:NZ*2 charmap[i-2*NZ-1] = i + Nodd end
    if Nodd == 1 charmap[0] = NZ+1 end
    return charmap
end

"""
    mw_column(::Type{<:SO{N}}, l::Int, aux::Bool) -> Vector{Int}

Return a maximal-weight tableau column for SO(N).

`l` is the column height. `aux` selects the alternate final column needed for
even SO(N) spinor-related conventions; when it applies, the top character is
negated to distinguish the auxiliary column.
"""
# This function should be modified
function mw_column(::Type{<:SO{N}}, l::Int, aux::Bool) where N
    col = collect(l:-1:1)
    if aux && l == div(N, 2) col[1] = -col[1] end
    return col
end


getdzs(::Type{SO{N}}) where N = [collect(getdz(SO{N}, i)) for i=1:nlops(SO{N})]

"""
    charlist(::Type{SO{N}}) -> Vector{Int}

Return the ordered SO(N) tableau character alphabet.

The list contains positive characters followed by negative characters. Odd
SO(N) appends `0` for the middle defining-vector weight.
"""
function charlist(::Type{SO{N}}) where N
    lst = vcat(collect(1:div(N, 2)), collect(-div(N, 2):-1))
    if N % 2 == 1 push!(lst, 0) end
    return lst
end

"""
    get_fops_std(::Type{SO{N}}) -> Vector{Dict{Int, Int}}

Return standard crystal lowering maps for SO(N) tableau characters.

Each dictionary describes one simple lowering operator. The last dictionary is
parity dependent: odd SO(N) lowers through the zero character, while even SO(N)
uses the two terminal vector characters.
"""
# f-operations defined for crystal of Tableau
function get_fops_std(::Type{SO{N}}) where N
    fops = Vector{Dict{Int, Int}}()
    NZ = div(N, 2); Nodd = (N % 2 == 1)
    for i in 1:NZ-1
        fop = Dict{Int, Int}()
        fop[i] = i + 1; fop[-i-1] = -i
        push!(fops, fop)
    end
    fop_last = Dict{Int, Int}()
    if Nodd fop_last[NZ] = 0; fop_last[0] = -NZ
    else fop_last[NZ-1] = -NZ; fop_last[NZ] = -NZ+1 end
    push!(fops, fop_last)
    return fops
end

"""
    qlab2mwz(::Type{SO{N}}, qlabel::NTuple{div(N,2), Int}) -> NTuple{div(N,2), Int}

Convert an SO(N) qlabel to maximal-weight z-coordinates.

`qlabel` uses the package's Dynkin-label convention. Odd and even SO(N) use
different formulas for the final coordinates, reflecting the different last
simple root and spinor-label structure.
"""
function qlab2mwz(::Type{SO{N}}, qlabel::NTuple{NZ, Int}) where {N, NZ}
    @assert div(N, 2) == NZ
    if N % 2 == 1
        z = fill(qlabel[NZ], NZ)
        for i in 1:NZ-1 z[NZ+1-i] += 2*sum(qlabel[i:NZ-1]) end
    else
        z = fill(sum(qlabel[NZ-1:NZ]), NZ)
        z[1] -= 2 * qlabel[NZ-1]
        for i in 1:NZ-2 z[NZ+1-i] += 2*sum(qlabel[i:NZ-2]) end
    end
    return Tuple(z)
end

"""
    getqlabel(::Type{SO{N}}, z::NTuple{div(N,2), Int}) -> NTuple{div(N,2), Int}

Convert maximal-weight z-coordinates back to an SO(N) qlabel.

`z` is expected to use the doubled-weight convention. Assertions check tuple
length and divisibility before `determine_last` reconstructs the parity-specific
final Dynkin labels.
"""
function getqlabel(::Type{SO{N}}, z::NTuple{NZ, Int}) where {N, NZ}
    Nodd = (N % 2 == 1)
    @assert NZ == nzops(SO{N})
    qlabel = zeros(Int, NZ)
    for i in 1:NZ-2
        wdiff = z[NZ+1-i] - z[NZ-i]
        @assert wdiff % 2 == 0
        qlabel[i] = div(wdiff, 2)
    end
    qlabel[NZ-1:NZ] = determine_last(SO{N}, z[1], z[2])
    return Tuple(qlabel)
end

"""
    determine_last(::Type{SO{N}}, z1::Int, z2::Int) -> Vector{Int}

Recover the final SO(N) qlabel coordinates from the first two z-coordinates.

`z1` and `z2` are doubled z-weight coordinates. Odd SO(N) returns the terminal
vector/spin coordinate pair; even SO(N) returns the two spinor-end coordinates.
"""
function determine_last(::Type{SO{N}}, z1::Int, z2::Int) where N
    Nodd = (N % 2 == 1)
    if Nodd
        a2 = z2 - z1; b = z1
        @assert a2 % 2 == 0; a = div(a2, 2)
    else
        a2, b2 = z1 + z2, z2 - z1
        @assert a2 % 2 == 0 && b2 % 2 == 0
        a = div(a2, 2); b = div(b2, 2)
    end
    return [a, b]
end

"""
    preprocess(::Type{SO{N}}, qlabel::Vector{Int}) -> Vector{Int}

Normalize an SO(N) qlabel before tableau construction.

`qlabel` is copied and transformed into the column-count convention used by the
crystal code. Odd SO(N) halves the final doubled spin coordinate. Even SO(N)
rewrites the last two coordinates into ordered spinor-column counts.
"""
function preprocess(::Type{SO{N}}, qlabel::Vector{Int}) where N
    Nodd = (N % 2 == 1); NZ = div(N, 2)
    @assert length(qlabel) == NZ
    res = copy(qlabel)
    if Nodd @assert qlabel[NZ] % 2 == 0; res[NZ] = div(qlabel[NZ], 2) 
    else 
        @assert (qlabel[NZ-1] + qlabel[NZ]) % 2 == 0
        mi, ma = minmax(qlabel[NZ-1], qlabel[NZ])
        res[NZ-1] = mi; res[NZ] = div(ma - mi, 2)
    end
    return res
end

get_auxarg(::Type{SO{N}}, qlabel::NTuple{NZ, Int}) where {N, NZ} = 
    N % 2 == 0 && qlabel[NZ-1] < qlabel[NZ]

"""
    less_weight(::Type{SO{N}}, w1::NTuple{div(N,2), Int}, w2::NTuple{div(N,2), Int}) -> Bool

Order SO(N) weight tuples for deterministic sector traversal.

`w1` and `w2` are compared lexicographically after reversing coordinate order,
matching the package's tableau-shape convention.
"""
# Weight comparision, inputs are a form of shape of the Young tableau
function less_weight(::Type{SO{N}}, w1::NTuple{NZ, Int}, w2::NTuple{NZ, Int}) where {N, NZ}
	@assert NZ == div(N, 2)
	return reverse(w1) < reverse(w2)
end
