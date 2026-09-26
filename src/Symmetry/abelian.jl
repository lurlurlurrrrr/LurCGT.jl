"""
    add_qn(Z{N}, q1, q2) -> Int

Add two `Z{N}` quantum numbers modulo `N`. The result is normalized to the
canonical range `0:N-1`, including when either input is negative.

`Z{N}` selects the cyclic symmetry and its modulus; `q1` and `q2` are the two
integer charges to combine.
"""
add_qn(::Type{Z{N}}, q1::Int, q2::Int) where N = mod(q1 + q2, N)
"""
    add_qn(U1, q1, q2) -> Int

Add two `U1` quantum numbers by ordinary integer addition. Unlike `Z{N}`
labels, `U1` labels are unbounded and are not reduced modulo any period.

`U1` selects the charge symmetry; `q1` and `q2` are the integer charges to
combine.
"""
add_qn(::Type{U1}, q1::Int, q2::Int) = q1 + q2
"""
    get_dualq(S, q) -> qdual

Return the q-label of the dual Abelian charge sector.

`S` selects either `U1` or `Z{N}`. `q` is its one-element q-label tuple. The
returned charge is negated for `U1` and negated modulo `N` for `Z{N}`.
"""
get_dualq(::Type{U1}, q::NTuple{1, Int}) = (-q[1],)
get_dualq(::Type{Z{N}}, q::NTuple{1, Int}) where N = (mod(-q[1], N),)
getqlabel(::Type{<:AbelianSymm}, q::Tuple{Int}) = q
qlab2mwz(::Type{<:AbelianSymm}, q::Tuple{Int}) = q

"""Return the number of weight operators for any Abelian symmetry family (`1`), which is also its q-label tuple width."""
nzops(::Type{<:AbelianSymm}) = 1
nlops(::Type{<:AbelianSymm}) = 0

totxt(::Type{U1}) = "U1"
totxt(::Type{Z{N}}) where N = "Z$N"
