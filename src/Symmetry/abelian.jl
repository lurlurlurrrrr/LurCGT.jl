add_qn(::Type{Z{N}}, q1::Int, q2::Int) where N = mod(q1 + q2, N)
add_qn(::Type{U1}, q1::Int, q2::Int) = q1 + q2
get_dualq(::Type{U1}, q::NTuple{1, Int}) = (-q[1],)
get_dualq(::Type{Z{N}}, q::NTuple{1, Int}) where N = (mod(-q[1], N),)
getqlabel(::Type{<:AbelianSymm}, q::Tuple{Int}) = q
qlab2mwz(::Type{<:AbelianSymm}, q::Tuple{Int}) = q

nzops(::Type{<:AbelianSymm}) = 1
nlops(::Type{<:AbelianSymm}) = 0

totxt(::Type{U1}) = "U1"
totxt(::Type{Z{N}}) where N = "Z$N"