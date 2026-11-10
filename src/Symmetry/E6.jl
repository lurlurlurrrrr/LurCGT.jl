# Basic metadata for the rank-six exceptional family.
isvalidsymm(::Type{E6}) = true
totxt(::Type{E6}) = "E6"
nlops(::Type{E6}) = 6
nzops(::Type{E6}) = 6

"""`fundamental_sequence(E6)` returns the fixed seed/product dependency order."""
fundamental_sequence(::Type{E6}) = [1, 6, 2, 3, 5, 4]

"""`minuscule_nodes(E6)` contains the dual 27-dimensional seeds."""
minuscule_nodes(::Type{E6}) = [1, 6]

"""`fundamental_sources(E6, k)` gives the fixed E6 product recipe."""
fundamental_sources(::Type{E6}, k::Int) = ((0,0), (1,6), (1,1), (3,1), (6,6), (0,0))[k]

"""`_exceptional_weight_rows(E6)` returns the fixed reversed scaled-root conversion."""
_exceptional_weight_rows(::Type{E6}) = ((2,3,4,6,5,4,),(4,6,8,12,10,5,),(6,9,12,18,12,6,),(5,6,10,12,8,4,),(3,6,6,9,6,3,),(4,3,5,6,4,2,),)

"""`_exceptional_cartan(S)` returns the E6 Cartan matrix in book node order."""
function _exceptional_cartan(::Type{S}) where {S<:E6}
    r = nzops(S); A = 2Matrix{Int}(I,r,r)
    for (i,j) in ((1,3),(2,4),(3,4),(4,5),(5,6))
        A[i,j] = A[j,i] = -1
    end
    return A
end

"""`qlab2mwz(S, q)` converts E6 Dynkin labels to reversed scaled-root weights."""
function qlab2mwz(::Type{S}, q::NTuple{R,Int}) where {S<:E6,R}
    R == nzops(S) || throw(DimensionMismatch("Wrong exceptional rank"))
    rows = _exceptional_weight_rows(S)
    return ntuple(i -> foldl(Base.Checked.checked_add,
        (Base.Checked.checked_mul(rows[i][j],q[j]) for j in 1:R); init=0), R)
end

"""`getqlabel(S, z)` exactly converts E6 lattice weights to signed Dynkin labels."""
function getqlabel(::Type{S}, z::NTuple{R,Int}) where {S<:E6,R}
    R == nzops(S) || throw(DimensionMismatch("Wrong exceptional rank"))
    u = reverse(z); d = 3
    neighbors = ((3,),(4,),(1,4),(2,3,5),(4,6),(5,))
    return ntuple(R) do i
        value = Base.Checked.checked_mul(2,u[i])
        for j in neighbors[i]
            value = Base.Checked.checked_sub(value,u[j])
        end
        rem(value,d)==0 || throw(ArgumentError("Weight is outside the integral lattice"))
        div(value,d)
    end
end

"""`getdz(S, i)` returns the E6 positive step at stored slot rank+1−i."""
getdz(::Type{S}, i::Int) where {S<:E6} = ntuple(j -> j == nzops(S)+1-i ? 3 : 0, nzops(S))

"""
    getdzs(S) -> Vector{Vector{Int}}

Return E6 simple-root steps in reversed scaled-root coordinates: color i
has entry d at slot rank+1-i, with d=3. These coordinate
unit vectors are not orthogonal Euclidean roots. Their root-space metric is
`P*A*P/d^2`, where A is the Cartan matrix and P reverses coordinates. Subtract
the returned step from a stored weight when lowering.
"""
getdzs(::Type{S}) where {S<:E6} = [collect(getdz(S,i)) for i in 1:nlops(S)]

"""`getsl_triv(S, RT)` returns empty E6 trivial lowering blocks."""
getsl_triv(::Type{S}, ::Type{RT}) where {S<:E6,RT} = ntuple(_ -> Dict{NTuple{nzops(S),Int},SparseMatrixCSC{RT}}(), nlops(S))

"""`getsz_triv(S)` returns the E6 zero-weight sector."""
getsz_triv(::Type{S}) where {S<:E6} = Dict(ntuple(_ -> 0,nzops(S)) => (1,1))

"""`less_weight(S, a, b)` compares E6 weights in reversed lexicographic order."""
less_weight(::Type{S}, a::Tuple, b::Tuple) where {S<:E6} = rev_less(a,b)

"""`get_dualq(E6, q)` exchanges the two pairs of dual fundamental nodes."""
get_dualq(::Type{E6}, q::NTuple{6,Int}) = (q[6],q[2],q[5],q[4],q[3],q[1])

"""`fundamental_qlabels(E6)` returns unit labels in fixed construction order."""
fundamental_qlabels(::Type{S}) where {S<:E6} = [_unit_label(S,k) for k in fundamental_sequence(S)]

"""`_has_minuscule_recipes(E6)` reports whether Part 3a recipes are available."""
_has_minuscule_recipes(::Type{E6}) = true

"""`fundamental_validation_rules(E6)` returns fixed ordered target/source checks, including family-specific cases."""
function fundamental_validation_rules(::Type{S}) where {S<:E6}
    return _fundamental_validation_rules(S,_unit_label(S,2))
end
