const WeightInfo{NZ} = Dict{Tuple{NTuple{NZ, Int}, NTuple{NZ, Int}}, Tuple{Int, Int}}

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

"""
    construct_recipe(S, RT, recipe::MinusculeRecipe; candidates=nothing, verbose)

Execute a minuscule seed recipe without saving it. Dispatch occurs outside the
numerical loops; the target label selects the atomic crystal alphabet.
"""
function construct_recipe(::Type{S}, ::Type{RT}, recipe::MinusculeRecipe;
        candidates=nothing, verbose=VERBOSE[]) where {S,RT<:Integer}
    _validate_recipe(S,recipe)
    return get_minuscule_irep(S,RT,findfirst(==(1),recipe.target))
end

"""
    construct_recipe(S, RT, recipe::ProductRecipe; candidates=nothing, verbose)

Recursively obtain ordered source irreps and extract the selected target copy.
Propagate incidental collection through dependencies without saving candidates.
"""
function construct_recipe(::Type{S}, ::Type{RT}, recipe::ProductRecipe;
        candidates=nothing, verbose=VERBOSE[]) where {S,RT<:Integer}
    _validate_recipe(S,recipe)
    left,right = map(q -> _recipe_label(S,q),recipe.sources)
    # Cache lookup erases coefficient-type inference; restore the requested storage type.
    rep_type = Irep{S,nlops(S),nzops(S),RT}
    x = getNsave_irep(S,RT,left; candidates,verbose)::rep_type
    y = getNsave_irep(S,RT,right; candidates,verbose)::rep_type
    Sl,Sz,_,wi = tensor_prod_reps(x,y;verbose)
    return get_irep_frommw(S,Sl,Sz,x.innerprod,y.innerprod,
        qlab2mwz(S,_recipe_label(S,recipe.target)),wi;
        copy_index=recipe.copy_index,candidates,verbose)
end

"""`construct_recipe(S, RT, recipe::AdjointRecipe)` rejects the deferred adjoint backend."""
function construct_recipe(::Type{S}, ::Type{RT}, recipe::AdjointRecipe;
        candidates=nothing, verbose=VERBOSE[]) where {S,RT<:Integer}
    _validate_recipe(S,recipe)
    throw(ArgumentError("The direct adjoint crystal backend is not implemented for $S"))
end

"""
    construct_fundamental(S, RT, k; candidates=nothing) -> Irep

Execute the requested hard-coded fundamental recipe, recursively loading only
its required sources. Optional candidates receives incidental nonfundamental
components; the caller saves the returned fundamental.
"""
function construct_fundamental(::Type{S}, ::Type{RT}, k::Int; candidates=nothing, verbose=VERBOSE[]) where {S,RT<:Integer}
    recipes = fundamental_recipes(S)
    validate_recipe_catalogue(S,recipes)
    1 <= k <= length(recipes) || throw(ArgumentError("Unknown fundamental node"))
    recipe = recipes[k]
    rep = construct_recipe(S,RT,recipe; candidates,verbose)
    rep.qlabel == _recipe_label(S,recipe.target) || error("Constructed recipe target differs")
    dimension(rep) == recipe.dimension || error("Constructed recipe dimension differs")
    return rep
end

"""
    _orthogonal_remainder(v, basis, metric) -> Vector

Project an integer candidate off mutually metric-orthogonal stored columns.
Clear denominators by integer cross multiplication and primitive rescaling;
no rational vector or nonexact division is introduced.
"""
function _orthogonal_remainder(v::Vector{RT}, basis::Vector{Vector{RT}}, metric::Matrix{RT}) where RT
    v = copy(v)
    for u in basis
        overlap = dot(u,metric*v)
        iszero(overlap) && continue
        norm = dot(u,metric*u)
        v = norm*v - overlap*u
        iszero(v) && return v
        v = div.(v,gcd(v))
    end
    return v
end

"""
    get_irep_frommw(S, Sl, Sz, metric1, metric2, mw, weight_info; copy_index=1, candidates=nothing, copies=nothing)

Extract the component with explicit highest weight `mw` from an ordered
product, constructing all preceding multiplets first. The caller converts its
target q-label with `qlab2mwz(S, qlabel)`. Select copies by exact projection of
standard basis vectors in index order. Optional `candidates` collects nonfundamental normalized irreps
without saving them. Optional `copies` collects every target copy before returning
the selected copy; errors propagate and absent targets are rejected.
"""
function get_irep_frommw(::Type{S},
        Sl::NTuple{NL,Dict{NTuple{NZ,Int},SparseMatrixCSC{RT}}},
        Sz::Dict{NTuple{NZ,Int},Tuple{Int,Int}},
        metric1::Dict{NTuple{NZ,Int},Matrix{RT}},
        metric2::Dict{NTuple{NZ,Int},Matrix{RT}},
        mw::NTuple{NZ,Int},
        winfo::Dict{NTuple{NZ,Int},WeightInfo{NZ}};
        copy_index::Int=1, candidates=nothing, copies=nothing, verbose=VERBOSE[]) where {S,NL,NZ,RT<:Integer}
    copy_index > 0 || throw(ArgumentError("Copy index must be positive"))
    weights = sorted_zvals(Sz)
    haskey(Sz,mw) || throw(ArgumentError("Target weight absent from product"))
    # The maximal Cartan component has a one-dimensional highest sector.
    if mw == first(weights) && copy_index == 1
        rep, _ = _irep_component(S,Sl,Sz,metric1,metric2,mw,RT[1],winfo; verbose)
        if candidates !== nothing && sum(rep.qlabel) > 1; push!(candidates,rep); end
        copies === nothing || push!(copies,rep)
        return rep
    end
    metrics = Dict{NTuple{NZ,Int},Matrix{RT}}()
    found = Dict(w => Vector{RT}[] for w in weights)
    target_copies = 0
    selected = nothing
    for w in weights
        rev_less(w,mw) && break
        d = Sz[w][2]-Sz[w][1]+1
        if !haskey(metrics,w)
            metrics[w] = get_innerprodmat(S,RT,w,Matrix{RT}(I,d,d),d,metric1,metric2,winfo; verbose)
        end
        for index in 1:d
            length(found[w]) == d && break
            v = zeros(RT,d); v[index] = 1
            v = _orthogonal_remainder(v,found[w],metrics[w])
            iszero(v) && continue
            all(>=(0),getqlabel(S,w)) || error("Residual product state is not highest weight")
            rep, bases = _irep_component(S,Sl,Sz,metric1,metric2,w,v,winfo; verbose)
            if candidates !== nothing && sum(rep.qlabel)>1; push!(candidates,rep); end
            if w == mw
                target_copies += 1
                if target_copies == copy_index
                    copies === nothing && return rep
                    selected = rep
                end
                copies === nothing || push!(copies,rep)
            end
            for z in sorted_zvals(bases)
                B = bases[z]; zd = size(B,1)
                if !haskey(metrics,z)
                    metrics[z] = get_innerprodmat(S,RT,z,Matrix{RT}(I,zd,zd),zd,metric1,metric2,winfo; verbose)
                end
                for j in axes(B,2)
                    u = _orthogonal_remainder(B[:,j],found[z],metrics[z])
                    iszero(u) && error("Dependent multiplet column")
                    push!(found[z],u)
                end
            end
        end
    end
    selected === nothing || return selected
    throw(ArgumentError("Target irrep or requested copy absent from product"))
end

"""
    same_irep_fields(a, b) -> Bool

Compare every mathematical Irep field exactly, including primitive inverse
metric factors. Ignore only the allocation-dependent `size_byte` field.
"""
function same_irep_fields(a::Irep, b::Irep)
    typeof(a) === typeof(b) || return false
    return a.Sl == b.Sl && a.Sz == b.Sz && a.innerprod == b.innerprod &&
        a.inv_innerprod == b.inv_innerprod && a.qlabel == b.qlabel && a.dimension == b.dimension
end
