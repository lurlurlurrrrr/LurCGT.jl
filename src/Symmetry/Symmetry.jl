"""Abstract supertype for LurCGT symmetry-family types. Use the concrete family type, not an instance, as API arguments."""
abstract type Symmetry end
"""Abstract subtype for Abelian symmetry families. Their sectors have trivial Clebsch-Gordan w-matrices and one q-label component."""
abstract type AbelianSymm <: Symmetry end
"""Abstract subtype for non-Abelian symmetry families, whose CGT sectors carry nontrivial w-matrices and Dynkin-style q-labels."""
abstract type NonabelianSymm <: Symmetry end

"""`Z{N}` is the cyclic Abelian group of order `N`; `N` is a positive type parameter and q-label arithmetic is modulo `N`."""
abstract type Z{N} <: AbelianSymm end
"""`U1` is the one-dimensional continuous Abelian symmetry family; q-labels are one-element integer tuples representing additive charge."""
abstract type U1 <: AbelianSymm end
"""`SU{N}` is the special-unitary non-Abelian family of rank `N - 1`; `N` selects the defining representation dimension."""
abstract type SU{N} <: NonabelianSymm end
"""`Sp{N}` is the compact symplectic non-Abelian family; `N` is its type parameter and q-labels use LurCGT's Dynkin convention."""
abstract type Sp{N} <: NonabelianSymm end
"""`SO{N}` is the special-orthogonal non-Abelian family; `N` is its defining-group dimension and q-labels use LurCGT's Dynkin convention."""
abstract type SO{N} <: NonabelianSymm end
"""`G2` is the exceptional rank-two non-Abelian symmetry family with two-component Dynkin q-labels."""
abstract type G2 <: NonabelianSymm end

# Stable family tags keep symmetry type hashes distinct from Julia's Type hash.
Base.hash(::Type{Z{N}}, h::UInt) where N = hash((0, N), h)
Base.hash(::Type{U1}, h::UInt) = hash((1,), h)
Base.hash(::Type{SU{N}}, h::UInt) where N = hash((2, N), h)
Base.hash(::Type{Sp{N}}, h::UInt) where N = hash((3, N), h)
Base.hash(::Type{SO{N}}, h::UInt) where N = hash((4, N), h)
Base.hash(::Type{G2}, h::UInt) = hash((5,), h)

include("abelian.jl")
include("SU.jl")
include("Sp.jl")
include("SO.jl")
include("G2.jl")

"""Return `true` when symmetry family type `S` is Abelian and `false` for non-Abelian families; `S` must be a subtype of `Symmetry`."""
isabelian(::Type{<:AbelianSymm}) = true
isabelian(::Type{<:NonabelianSymm}) = false

isvalidsymm(::Any) = false

nlops(::Any) = 0
"""Return the number of commuting weight (`z`) operators for symmetry type `S`; this is the q-label tuple width used by its representations."""
nzops(::Any) = 0

getsr(::Any) = error("Not implemented")
getszdiag(::Any) = error("Not implemented")

"""
    get_dualq(S, q) -> qdual

Return the q-label of the dual irrep of non-Abelian symmetry family `S`.
For `SU`, duality reverses the Dynkin-label tuple; the currently supported
non-Abelian families are self-dual in this q-label representation.
"""
get_dualq(::Type{S}, q::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ} =
    S<:SU ? reverse(q) : q

include("crystal.jl")
