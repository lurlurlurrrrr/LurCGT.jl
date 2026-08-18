add_qn(::Type{Z{N}}, q1::Int, q2::Int) where N = mod(q1 + q2, N)
"""Add Abelian quantum numbers `q1` and `q2` for the specified family. For `U1` this is ordinary integer addition; other methods may impose modular rules."""
add_qn(::Type{U1}, q1::Int, q2::Int) = q1 + q2
"""Return the dual q-label of Abelian label `q`. For `U1`, `q` must be a one-element integer tuple and the returned charge is negated."""
get_dualq(::Type{U1}, q::NTuple{1, Int}) = (-q[1],)
get_dualq(::Type{Z{N}}, q::NTuple{1, Int}) where N = (mod(-q[1], N),)
getqlabel(::Type{<:AbelianSymm}, q::Tuple{Int}) = q
qlab2mwz(::Type{<:AbelianSymm}, q::Tuple{Int}) = q

"""Return the number of weight operators for any Abelian symmetry family (`1`), which is also its q-label tuple width."""
nzops(::Type{<:AbelianSymm}) = 1
nlops(::Type{<:AbelianSymm}) = 0

"""Return the stable filename/cache text representation of symmetry type `S`, such as `"U1"` or `"Z2"`."""
totxt(::Type{U1}) = "U1"
totxt(::Type{Z{N}}) where N = "Z$N"
