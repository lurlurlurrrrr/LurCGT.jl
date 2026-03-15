abstract type Symmetry end
abstract type AbelianSymm <: Symmetry end
abstract type NonabelianSymm <: Symmetry end

abstract type Z{N} <: AbelianSymm end
abstract type U1 <: AbelianSymm end
abstract type SU{N} <: NonabelianSymm end
abstract type Sp{N} <: NonabelianSymm end
abstract type SO{N} <: NonabelianSymm end

include("abelian.jl")
include("SU.jl")
include("Sp.jl")
include("SO.jl")

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