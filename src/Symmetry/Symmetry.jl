"""Abstract supertype for symmetry-family types. Use the concrete family type, not an instance, as API arguments."""
abstract type Symmetry end
"""Abstract subtype for Abelian symmetry families. Their sectors have trivial Clebsch-Gordan coefficients and one q-label component."""
abstract type AbelianSymm <: Symmetry end
"""Abstract subtype for non-Abelian symmetry families, whose Clebsch-Gordan coefficients may be nontrivial and whose q-labels use Dynkin-style conventions."""
abstract type NonabelianSymm <: Symmetry end

"""`Z{N}` is the cyclic Abelian group of order `N`; `N` is a positive type parameter and q-label arithmetic is modulo `N`."""
abstract type Z{N} <: AbelianSymm end
"""`U1` is the Abelian charge-symmetry family; q-labels are one-element integer tuples representing unbounded additive charge."""
abstract type U1 <: AbelianSymm end
"""``\\mathrm{SU}(N)`` is the special-unitary non-Abelian family of rank `N - 1`; its Julia type parameter `N` selects the defining representation dimension."""
abstract type SU{N} <: NonabelianSymm end
"""``\\mathrm{Sp}(N)`` is the compact symplectic non-Abelian family; its Julia type parameter `N` and q-labels use this library's Dynkin convention."""
abstract type Sp{N} <: NonabelianSymm end
"""``\\mathrm{SO}(N)`` denotes the orthogonal Lie-algebra family; even `N` includes Spin(N) half-spin irreps, which need not descend to SO(N). Its Julia type parameter `N` is the vector dimension and q-labels are Dynkin labels."""
abstract type SO{N} <: NonabelianSymm end
"""`G2` is the exceptional rank-two non-Abelian symmetry family with two-component Dynkin q-labels."""
abstract type G2 <: NonabelianSymm end
"""`F4` is the exceptional rank-four non-Abelian symmetry family."""
abstract type F4 <: NonabelianSymm end
"""`E6` is the exceptional rank-six non-Abelian symmetry family."""
abstract type E6 <: NonabelianSymm end
"""`E7` is the exceptional rank-seven non-Abelian symmetry family."""
abstract type E7 <: NonabelianSymm end
"""`E8` is the exceptional rank-eight non-Abelian symmetry family."""
abstract type E8 <: NonabelianSymm end

"""
    totxt(::Type{S}) -> String where {S<:Symmetry}

Return the stable text key for symmetry family `S`, used in cache and database
names. Supported families include `Z{N}`, `U1`, `SU{N}`, `SO{N}`, `Sp{N}`, and
`G2`; for example, `totxt(SU{2})` returns `"SU2"`.

Specialized methods provide compact keys for the supported symmetry families.
For any other subtype of `Symmetry`, the general method returns `string(S)`.
"""
totxt(::Type{S}) where {S<:Symmetry} = string(S)

# Stable family tags keep symmetry type hashes distinct from Julia's Type hash.
Base.hash(::Type{Z{N}}, h::UInt) where N = hash((0, N), h)
Base.hash(::Type{U1}, h::UInt) = hash((1,), h)
Base.hash(::Type{SU{N}}, h::UInt) where N = hash((2, N), h)
Base.hash(::Type{Sp{N}}, h::UInt) where N = hash((3, N), h)
Base.hash(::Type{SO{N}}, h::UInt) where N = hash((4, N), h)
Base.hash(::Type{G2}, h::UInt) = hash((5,), h)
Base.hash(::Type{F4}, h::UInt) = hash((6,), h)
Base.hash(::Type{E6}, h::UInt) = hash((7,), h)
Base.hash(::Type{E7}, h::UInt) = hash((8,), h)
Base.hash(::Type{E8}, h::UInt) = hash((9,), h)

"""
    FundamentalRecipe

Abstract construction metadata for one fundamental irrep. Catalogues are
`Vector{FundamentalRecipe}`; dispatch occurs once at the construction boundary.
"""
abstract type FundamentalRecipe end

"""
    MinusculeRecipe(target::Vector{Int}, dimension::Int)

Describe a minuscule Weyl-orbit seed with its expected dimension. The owned
label is read-only after catalogue construction; no numerical data is stored.
"""
struct MinusculeRecipe <: FundamentalRecipe
    target::Vector{Int}
    dimension::Int
end

"""
    AdjointRecipe(target::Vector{Int}, dimension::Int)

Describe a direct adjoint seed. Numerical execution remains unavailable until
the adjoint backend is implemented. The owned target label is read-only.
"""
struct AdjointRecipe <: FundamentalRecipe
    target::Vector{Int}
    dimension::Int
end

"""
    ProductRecipe(target, sources::NTuple{2,Vector{Int}}, copy_index, dimension)

Describe the selected copy of a target in an ordered source product. Labels
are owned read-only integer vectors; source order is never sorted. Copy indices
are positive and refer to deterministic numerical highest-vector selection.
"""
struct ProductRecipe <: FundamentalRecipe
    target::Vector{Int}
    sources::NTuple{2,Vector{Int}}
    copy_index::Int
    dimension::Int
end

"""
    IrrepValidationRule(target::Vector{Int}, sources::NTuple{2,Vector{Int}})

Compare every target copy extracted from the ordered sources with the canonical
irrep. Labels are owned read-only vectors. An absent target fails validation.
"""
struct IrrepValidationRule
    target::Vector{Int}
    sources::NTuple{2,Vector{Int}}
end

"""
    fundamental_recipes(S) -> Vector{FundamentalRecipe}

Return owned fundamental construction metadata indexed by node. An empty
catalogue means no numerical construction backend is available for `S`.
Symmetry-specific methods supply supported catalogues without constructing irreps.
"""
fundamental_recipes(::Type) = FundamentalRecipe[]

"""`_unit_label(S, k)` returns the rank-length unit Dynkin label at k."""
_unit_label(::Type{S}, k::Int) where S = ntuple(i -> Int(i == k), nzops(S))

"""
    _recipe_label(S, label::Vector{Int}) -> NTuple{nzops(S),Int}

Check rank and nonnegative Dynkin entries before converting recipe metadata to
the existing numerical API. Does not mutate the label.
"""
function _recipe_label(::Type{S}, label::Vector{Int}) where S
    length(label) == nzops(S) || throw(ArgumentError("Wrong recipe label rank for $S"))
    all(>=(0),label) || throw(ArgumentError("Negative recipe label for $S"))
    return Tuple(label)::NTuple{nzops(S),Int}
end

"""`_validate_recipe(S, recipe)` checks the common seed metadata without construction."""
function _validate_recipe(::Type{S}, recipe::FundamentalRecipe) where S
    _recipe_label(S,recipe.target)
    sum(recipe.target) == 1 || throw(ArgumentError("Recipe target must be fundamental"))
    recipe.dimension > 0 || throw(ArgumentError("Recipe dimension must be positive"))
    return nothing
end

"""`_validate_recipe(S, recipe::ProductRecipe)` additionally checks ordered sources and copy index."""
function _validate_recipe(::Type{S}, recipe::ProductRecipe) where S
    invoke(_validate_recipe,Tuple{Type{S},FundamentalRecipe},S,recipe)
    foreach(q -> _recipe_label(S,q),recipe.sources)
    recipe.copy_index > 0 || throw(ArgumentError("Recipe copy index must be positive"))
    return nothing
end

"""
    validate_recipe_catalogue(S, recipes) -> nothing

Check one correctly indexed fundamental entry per node. Dependencies are
hard-coded and well founded; no runtime graph or cycle machinery is added.
"""
function validate_recipe_catalogue(::Type{S}, recipes::Vector{FundamentalRecipe}) where S
    length(recipes) == nzops(S) || throw(ArgumentError("Incomplete fundamental catalogue"))
    for (k,recipe) in enumerate(recipes)
        _validate_recipe(S,recipe)
        _recipe_label(S,recipe.target) == _unit_label(S,k) || throw(ArgumentError("Misindexed fundamental recipe"))
    end
    return nothing
end

"""
    _irrep_validation_rules(S, adjoint) -> Vector{IrrepValidationRule}

Build shared adjoint-in-fundamental-times-dual checks in historical label order.
Family methods supply the adjoint label and append their special cases.
"""
function _irrep_validation_rules(::Type{S}, adjoint::NTuple{R,Int}) where {S,R}
    rules = IrrepValidationRule[]
    for left in fundamental_qlabels(S)
        right = get_dualq(S,left)
        left > right && continue
        push!(rules,IrrepValidationRule(collect(adjoint),(collect(left),collect(right))))
    end
    return rules
end

include("abelian.jl")
include("SU.jl")
include("Sp.jl")
include("SO.jl")
include("G2.jl")
include("F4.jl")
include("E6.jl")
include("E7.jl")
include("E8.jl")

"""
    isabelian(S) -> Bool

Return whether symmetry family type `S` is Abelian.

`S` must be a subtype of `Symmetry`, for example `U1`, `Z{N}`, or `SU{N}`.
The result is `true` for `AbelianSymm` families and `false` for
`NonabelianSymm` families.
"""
isabelian(::Type{<:AbelianSymm}) = true
isabelian(::Type{<:NonabelianSymm}) = false

isvalidsymm(::Any) = false

nlops(::Any) = 0
"""
    nzops(S) -> Int

Return the q-label width for symmetry family type `S`.

`S` is a symmetry family type. The returned count is the number of commuting
weight (`z`) operators and therefore the required length of every q-label tuple
for that family.
"""
nzops(::Any) = 0

getsr(::Any) = error("Not implemented")
getszdiag(::Any) = error("Not implemented")

"""
    get_dualq(S, q) -> qdual

Return the q-label of the dual irrep of non-Abelian symmetry family `S`.
For `SU`, duality reverses Dynkin labels. E6 and odd-rank even-orthogonal
families specialize the node permutation; the remaining supported families
are self-dual.

`S` selects the non-Abelian symmetry family and `q` is one q-label tuple for
that family.
"""
get_dualq(::Type{S}, q::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ} =
    S<:SU ? reverse(q) : q

include("crystal.jl")
