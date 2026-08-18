"""
    isvalidsymm(::Type{G2}) -> Bool

Return whether the exceptional `G2` symmetry is supported.

`G2` has fixed rank two and a seven-dimensional defining representation, so no
size parameter is required.
"""
isvalidsymm(::Type{G2}) = true

"""
    totxt(::Type{G2}) -> String

Return the stable text key used for G2 file, folder, and database names.
"""
totxt(::Type{G2}) = "G2"

"""
    defirepdim(::Type{G2}) -> Int

Return the dimension of the defining G2 representation.

The package uses the standard seven-dimensional defining representation.
"""
defirepdim(::Type{G2}) = 7

"""
    nlops(::Type{G2}) -> Int
    nzops(::Type{G2}) -> Int

Return the number of simple lowering operators and z-weight coordinates.

Both are fixed at two for rank-two `G2`.
"""
nlops(::Type{G2}) = 2
nzops(::Type{G2}) = 2

maxalt(::Type{G2}) = 2

const G2_DEF_CHARS = (1, 2, 3, 0, -3, -2, -1)
const G2_DEF_WEIGHTS = (
    (1, 1),
    (-1, 1),
    (2, 0),
    (0, 0),
    (-2, 0),
    (1, -1),
    (-1, -1),
)
const G2_FOPS = (
    Dict(1 => 2, 3 => 0, 0 => -3, -2 => -1),
    Dict(2 => 3, -3 => -2),
)

getsl_triv(::Type{G2}, ::Type{RT}) where {RT<:Number} =
    Tuple(Dict{NTuple{2, Int}, SparseMatrixCSC{RT}}() for _=1:2)

"""
    getsz_triv(::Type{G2}) -> Dict

Return the trivial G2 weight-sector dictionary.

The only sector is z-weight `(0, 0)` with a one-dimensional basis range.
"""
getsz_triv(::Type{G2}) = Dict((0, 0) => (1, 1))

"""
    getsz_def(::Type{G2}, i::Int) -> NTuple{2, Int}

Return the z-weight of defining-representation basis state `i`.

`i` indexes the fixed seven-entry `G2_DEF_WEIGHTS` table.
"""
getsz_def(::Type{G2}, i::Int) = G2_DEF_WEIGHTS[i]

getsz_def(::Type{G2}) = Dict(getsz_def(G2, i) => (i, i) for i=1:defirepdim(G2))

getsz_def_vec(::Type{G2}) = [collect(getsz_def(G2, i)) for i=1:defirepdim(G2)]

"""
    getsl_def(::Type{G2}, ::Type{RT}) -> NTuple{2, Dict}

Build sparse lowering-operator blocks for the defining G2 representation.

`RT` is the scalar type of the sparse matrices. The returned tuple has one
dictionary for each simple lowering operator. Source sectors and transitions
come from `G2_FOPS`; the zero-weight block has coefficient `2` in the first
operator under the package's normalization.
"""
function getsl_def(::Type{G2}, ::Type{RT}) where {RT<:Number}
    def_sl = Tuple(Dict{NTuple{2, Int}, SparseMatrixCSC{RT}}() for _=1:2)
    char_weight = Dict(char => weight for (char, weight) in zip(G2_DEF_CHARS, G2_DEF_WEIGHTS))
    for (lop, fop) in enumerate(G2_FOPS)
        for (src, _) in fop
            val = spzeros(RT, 1, 1)
            val[1, 1] = RT(1)
            def_sl[lop][char_weight[src]] = val
        end
    end
    def_sl[1][(0, 0)] = RT[2;;]
    return def_sl
end

"""
    getdz(::Type{G2}, lop::Int) -> NTuple{2, Int}

Return the z-weight change for G2 lowering operator `lop`.

`lop` must be `1` or `2`. The returned values are the fixed simple-root
coordinate changes in this representation convention.
"""
function getdz(::Type{G2}, lop::Int)
    @assert 1 <= lop <= 2
    return lop == 1 ? (2, 0) : (-3, 1)
end

"""
    crystal_chars_map(::Type{G2}) -> Dict{Int, Int}

Map G2 tableau characters to defining-representation basis positions.

The mapping follows the fixed character ordering in `G2_DEF_CHARS`.
"""
crystal_chars_map(::Type{G2}) = Dict(char => i for (i, char) in enumerate(G2_DEF_CHARS))

"""
    mw_column(::Type{G2}, l::Int) -> Vector{Int}

Return a canonical maximal-weight tableau column for G2.

`l` must be `1` or `2`, matching the two fundamental representations. The
returned character list seeds tableau/crystal generation.
"""
function mw_column(::Type{G2}, l::Int)
    @assert 1 <= l <= 2
    return l == 1 ? [1] : [2, 1]
end

getdzs(::Type{G2}) = [collect(getdz(G2, i)) for i=1:nlops(G2)]

charlist(::Type{G2}) = collect(G2_DEF_CHARS)

"""
    get_fops_std(::Type{G2}) -> Vector{Dict{Int, Int}}

Return copies of the fixed G2 crystal lowering maps.

Copies are returned so callers may manipulate the dictionaries without mutating
the global `G2_FOPS` constants.
"""
get_fops_std(::Type{G2}) = [copy(G2_FOPS[1]), copy(G2_FOPS[2])]

"""
    qlab2mwz(::Type{G2}, qlabel::NTuple{2, Int}) -> NTuple{2, Int}

Convert a G2 highest-weight qlabel to maximal-weight z-coordinates.

`qlabel == (a, b)` uses fundamental-weight coordinates and maps to
`(a, a + 2b)` in the package's z-weight convention.
"""
function qlab2mwz(::Type{G2}, qlabel::NTuple{NZ, Int}) where {NZ}
    @assert NZ == 2
    a, b = qlabel
    return (a, a + 2 * b)
end

"""
    getqlabel(::Type{G2}, z::NTuple{2, Int}) -> NTuple{2, Int}

Convert maximal-weight z-coordinates back to a G2 qlabel.

`z` must satisfy `z[2] - z[1]` even, which is asserted before returning the
fundamental-weight coordinates.
"""
function getqlabel(::Type{G2}, z::NTuple{NZ, Int}) where {NZ}
    @assert NZ == 2
    z1, z2 = z
    @assert iseven(z2 - z1)
    return (z1, div(z2 - z1, 2))
end

"""
    less_weight(::Type{G2}, w1::NTuple{2, Int}, w2::NTuple{2, Int}) -> Bool

Order two G2 weight tuples for deterministic sector traversal.

Weights are compared lexicographically after reversing coordinate order,
matching the convention used by the other non-Abelian symmetry files.
"""
function less_weight(::Type{G2}, w1::NTuple{NZ, Int}, w2::NTuple{NZ, Int}) where {NZ}
    @assert NZ == 2
    return reverse(w1) < reverse(w2)
end

"""
    fundamental_qlabels(::Type{G2}) -> Vector{Tuple{Int, Int}}

Return the two fundamental G2 qlabels.
"""
fundamental_qlabels(::Type{G2}) = [(1, 0), (0, 1)]

"""
    getdefirep(::Type{G2}, ::Type{RT}) -> Irep{G2}

Construct the defining G2 irrep in the package's sparse block format.

`RT` chooses the scalar type. The returned `Irep` contains lowering blocks,
defining z-sector ranges, inner products, inverse inner products, qlabel
`(1, 0)`, and dimension `7`. The zero-weight inner product uses the special
G2 normalization stored explicitly in this method.
"""
function getdefirep(::Type{G2}, ::Type{RT}) where {RT<:Number}
    def_sl = getsl_def(G2, RT)
    def_sz = getsz_def(G2)
    innerprod = get_identity_innerprod(G2, RT, def_sz)
    innerprod[(0, 0)] = RT[2;;]
    inv_innerprod = get_identity_invinprod(G2, RT, def_sz)
    inv_innerprod[(0, 0)] = (RT[1;;], 1//2)
    dim = defirepdim(G2)
    return Irep{G2, 2, 2, RT}(def_sl, def_sz, innerprod, inv_innerprod, (1, 0), dim)
end
