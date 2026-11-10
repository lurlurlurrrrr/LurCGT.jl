# Sl : lowering operators, Sz : z-operator
# norm2s : square norm values of basis element.
# For defining irep, values of innerprod are all [[1]], 1*1 matrix with entry 1
# NL(NZ) means number of lowering(z-) operators

# For larger irep, they are positive integer
# and matrix element of Sl, Sz are integer w.r.t such basis
# This may not be irreducible representation, for example, product of two irreps
"""
    Irep

An irreducible representation of symmetry `S`, stored in a weight basis. `Sl`
contains lowering-operator blocks, `Sz` gives basis ranges per weight,
`innerprod` and `inv_innerprod` store the corresponding metric data, `qlabel`
is the highest-weight label, and `dimension` is total Hilbert-space dimension.
`size_byte` is an internal cache-eviction value. Constructed values are
immutable; use cache-building functions instead of changing fields.

# Fields

- `Sl`: lowering-operator sparse blocks grouped by weight.
- `Sz`: inclusive basis-index ranges for each weight.
- `innerprod`, `inv_innerprod`: metric blocks and exact inverse-normalization data.
- `qlabel`: highest-weight/Dynkin label.
- `dimension`: total representation dimension.
- `size_byte`: cached memory footprint for LRU eviction.
"""
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

"""
    isless(r1::Irep, r2::Irep) -> Bool

Order irreps of the same symmetry and scalar type by qlabel.

This supports deterministic sorting of cached/generated representation lists.
"""
Base.isless(r1::Irep{S, NL, NZ, RT}, r2::Irep{S, NL, NZ, RT}) where {S, NL, NZ, RT} =
r1.qlabel < r2.qlabel

# Custom pretty printing for Irep type
"""
    totxt(rep::Irep) -> String

Return a compact text label for a non-Abelian irrep instance.

The label includes the symmetry text key and highest-weight qlabel.
"""
totxt(rep::Irep{S, NL, NZ, RT}) where {S<:NonabelianSymm, NL, NZ, RT} = 
    "irep($(totxt(S)), $(rep.qlabel))"

# get max weight from other fields and complete Irep struct
"""
    Irep(::Type{S}, ::Type{RT}, Sl, Sz, innerprod, inv_innerprod, dim) -> Irep

Construct an irrep from weight-basis storage and infer its qlabel.

`Sl` is the tuple of lowering-operator block dictionaries, `Sz` maps z-weights
to inclusive global basis ranges, `innerprod` stores metric blocks,
`inv_innerprod` stores inverse metric blocks plus exact rational factors, and
`dim` is total representation dimension. The highest z-weight in `Sz` is
converted to the qlabel by `getqlabel(S, maxz)`.
"""
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

"""
    dimension(S, q)
    dimension(rep)

Return the representation dimension. In `dimension(S, q)`, `S` is a symmetry
family type and `q` is its q-label tuple; the representation is loaded or
computed for non-Abelian families and Abelian dimensions are one. In
`dimension(rep)`, the cached dimension field is returned without allocation.
"""
dimension(::Type{S}, q::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ} =
dimension(getNsave_irep(S, BigInt, q))

dimension(::Type{S}, q::NTuple{NZ, Int}) where {S<:AbelianSymm, NZ} = 1

dimension(rep::Irep) = rep.dimension

"""
    dim_from_innerprod(innerprod) -> Int

Compute representation dimension from metric block sizes.

`innerprod[z]` is a square matrix for one z-weight sector. The result is the sum
of row counts across all sectors.
"""
function dim_from_innerprod(innerprod::Dict{NTuple{NZ, Int}, Matrix{RT}}) where {NZ, RT}
    dim = 0
    # Sum the dimensions of each inner product matrix)
    for (_, mat) in innerprod dim += size(mat, 1) end 
    return dim
end

"""
    get_identity_innerprod(::Type{S}, ::Type{RT}, Sz)

Create identity metric blocks matching an `Sz` sector dictionary.

`S` selects qlabel width, `RT` chooses matrix element type, and `Sz[z] == (a,b)`
defines each sector dimension as `b - a + 1`.
"""
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

"""
    get_identity_invinprod(::Type{S}, ::Type{RT}, Sz)

Create inverse-identity metric blocks matching an `Sz` sector dictionary.

Each returned value is `(I, 1//1)`, using the same sector sizes as
`get_identity_innerprod`. The rational factor field keeps the shape compatible
with generated irreps whose inverse metrics need denominator factors.
"""
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
    getdefirep(::Type{S}, ::Type{RT}) -> Irep

Construct the defining representation of a non-Abelian symmetry.

`S` supplies symmetry-specific defining lowering blocks and z-sector ranges via
`getsl_def` and `getsz_def`; `RT` chooses metric scalar type. Inner products are
identity blocks for the defining representation.
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
"""
    gettrivirep(::Type{S}, ::Type{RT}) -> Irep

Construct the one-dimensional trivial representation of `S`.

The lowering-operator tuple and zero z-sector dictionary come from
`getsl_triv` and `getsz_triv`, and all metric blocks are identity blocks.
"""
function gettrivirep(::Type{S}, ::Type{RT}) where {S<:NonabelianSymm, RT<:Number}
    Sl = getsl_triv(S, RT); Sz = getsz_triv(S)
    inner_prod = get_identity_innerprod(S, RT, Sz)
    inv_innerprod = get_identity_invinprod(S, RT, Sz)
    # For trivial representation, Tmats is equal to innerprod mat
    dim = dim_from_innerprod(inner_prod)
    return Irep(S, RT, Sl, Sz, inner_prod, inv_innerprod, dim)
end

"""
    show(io::IO, rep::Irep)

Print a compact irrep summary containing symmetry name, qlabel, and dimension.
"""
function Base.show(io::IO, rep::Irep{S}) where {S<:Symmetry}
    print(io, "$(totxt(S)) Irep: qlabel=$(rep.qlabel), dimension=$(rep.dimension)")
end

"""
    symm(rep::Irep) -> Type{<:Symmetry}

Return the symmetry type parameter of an irrep.
"""
symm(rep::Irep{S}) where {S<:Symmetry} = S
