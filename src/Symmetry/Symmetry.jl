"""Abstract supertype for LurCGT symmetry-family types. Use the concrete family type, not an instance, as API arguments."""
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
"""``\\mathrm{Sp}(N)`` is the compact symplectic non-Abelian family; its Julia type parameter `N` and q-labels use LurCGT's Dynkin convention."""
abstract type Sp{N} <: NonabelianSymm end
"""``\\mathrm{SO}(N)`` is the special-orthogonal non-Abelian family; its Julia type parameter `N` is the defining-group dimension and q-labels use LurCGT's Dynkin convention."""
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
For `SU`, duality reverses the Dynkin-label tuple; the currently supported
non-Abelian families are self-dual in this q-label representation.

`S` selects the non-Abelian symmetry family and `q` is one q-label tuple for
that family.
"""
get_dualq(::Type{S}, q::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ} =
    S<:SU ? reverse(q) : q

include("crystal.jl")
