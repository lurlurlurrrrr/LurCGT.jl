abstract type Symmetry end
abstract type AbelianSymm <: Symmetry end
abstract type NonabelianSymm <: Symmetry end

abstract type Z{N} <: AbelianSymm end
abstract type U1 <: AbelianSymm end
abstract type SU{N} <: NonabelianSymm end
abstract type Sp{N} <: NonabelianSymm end
abstract type SO{N} <: NonabelianSymm end
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

isabelian(::Type{<:AbelianSymm}) = true
isabelian(::Type{<:NonabelianSymm}) = false

isvalidsymm(::Any) = false

nlops(::Any) = 0
nzops(::Any) = 0

getsr(::Any) = error("Not implemented")
getszdiag(::Any) = error("Not implemented")

get_dualq(::Type{S}, q::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ} = 
    S<:SU ? reverse(q) : q

include("crystal.jl")
