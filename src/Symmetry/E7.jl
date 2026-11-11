# Basic metadata for the rank-seven exceptional family.
isvalidsymm(::Type{E7}) = true
totxt(::Type{E7}) = "E7"
nlops(::Type{E7}) = 7
nzops(::Type{E7}) = 7

"""`_exceptional_weight_rows(E7)` returns the fixed reversed scaled-root conversion."""
_exceptional_weight_rows(::Type{E7}) = ((2,3,4,6,5,4,3,),(4,6,8,12,10,8,4,),(6,9,12,18,15,10,5,),(8,12,16,24,18,12,6,),(6,8,12,16,12,8,4,),(4,7,8,12,9,6,3,),(4,4,6,8,6,4,2,),)

"""`_exceptional_cartan(S)` returns the E7 Cartan matrix in book node order."""
function _exceptional_cartan(::Type{S}) where {S<:E7}
    r = nzops(S); A = 2Matrix{Int}(I,r,r)
    for (i,j) in ((1,3),(2,4),(3,4),(4,5),(5,6))
        A[i,j] = A[j,i] = -1
    end
    A[6,7] = A[7,6] = -1
    return A
end

"""`qlab2mwz(S, q)` converts E7 Dynkin labels to reversed scaled-root weights."""
function qlab2mwz(::Type{S}, q::NTuple{R,Int}) where {S<:E7,R}
    R == nzops(S) || throw(DimensionMismatch("Wrong exceptional rank"))
    rows = _exceptional_weight_rows(S)
    return ntuple(i -> foldl(Base.Checked.checked_add,
        (Base.Checked.checked_mul(rows[i][j],q[j]) for j in 1:R); init=0), R)
end

"""`getqlabel(S, z)` exactly converts E7 lattice weights to signed Dynkin labels."""
function getqlabel(::Type{S}, z::NTuple{R,Int}) where {S<:E7,R}
    R == nzops(S) || throw(DimensionMismatch("Wrong exceptional rank"))
    u = reverse(z); d = 2
    neighbors = ((3,),(4,),(1,4),(2,3,5),(4,6),(5,7),(6,))
    return ntuple(R) do i
        value = Base.Checked.checked_mul(2,u[i])
        for j in neighbors[i]
            value = Base.Checked.checked_sub(value,u[j])
        end
        rem(value,d)==0 || throw(ArgumentError("Weight is outside the integral lattice"))
        div(value,d)
    end
end

"""`getdz(S, i)` returns the E7 positive step at stored slot rank+1−i."""
getdz(::Type{S}, i::Int) where {S<:E7} = ntuple(j -> j == nzops(S)+1-i ? 2 : 0, nzops(S))

"""
    getdzs(S) -> Vector{Vector{Int}}

Return E7 simple-root steps in reversed scaled-root coordinates: color i
has entry d at slot rank+1-i, with d=2. These coordinate
unit vectors are not orthogonal Euclidean roots. Their root-space metric is
`P*A*P/d^2`, where A is the Cartan matrix and P reverses coordinates. Subtract
the returned step from a stored weight when lowering.
"""
getdzs(::Type{S}) where {S<:E7} = [collect(getdz(S,i)) for i in 1:nlops(S)]

"""`getsl_triv(S, RT)` returns empty E7 trivial lowering blocks."""
getsl_triv(::Type{S}, ::Type{RT}) where {S<:E7,RT} = ntuple(_ -> Dict{NTuple{nzops(S),Int},SparseMatrixCSC{RT}}(), nlops(S))

"""`getsz_triv(S)` returns the E7 zero-weight sector."""
getsz_triv(::Type{S}) where {S<:E7} = Dict(ntuple(_ -> 0,nzops(S)) => (1,1))

"""`less_weight(S, a, b)` compares E7 weights in reversed lexicographic order."""
less_weight(::Type{S}, a::Tuple, b::Tuple) where {S<:E7} = rev_less(a,b)

"""`fundamental_qlabels(E7)` preserves historical label enumeration without scheduling construction."""
function fundamental_qlabels(::Type{S}) where {S<:E7}
    recipes = fundamental_recipes(S)
    return [_recipe_label(S,recipes[k].target) for k in (7,1,6,2,5,3,4)]
end

"""`irrep_validation_rules(E7)` returns fixed ordered target/source checks, including family-specific cases."""
function irrep_validation_rules(::Type{S}) where {S<:E7}
    return _irrep_validation_rules(S,_unit_label(S,1))
end

"""
    fundamental_recipes(E7) -> Vector{FundamentalRecipe}

Return owned metadata indexed by fundamental node. Dependencies are recursive;
no numerical representations or crystal vertices are constructed here.
"""
function fundamental_recipes(::Type{S}) where {S<:E7}
    return FundamentalRecipe[
        ProductRecipe(collect(_unit_label(S,1)), (collect(_unit_label(S,7)),collect(_unit_label(S,7))),1,133),
        ProductRecipe(collect(_unit_label(S,2)), (collect(_unit_label(S,1)),collect(_unit_label(S,7))),1,912),
        ProductRecipe(collect(_unit_label(S,3)), (collect(_unit_label(S,1)),collect(_unit_label(S,1))),1,8645),
        ProductRecipe(collect(_unit_label(S,4)), (collect(_unit_label(S,5)),collect(_unit_label(S,7))),1,365750),
        ProductRecipe(collect(_unit_label(S,5)), (collect(_unit_label(S,6)),collect(_unit_label(S,7))),1,27664),
        ProductRecipe(collect(_unit_label(S,6)), (collect(_unit_label(S,7)),collect(_unit_label(S,7))),1,1539),
        MinusculeRecipe(collect(_unit_label(S,7)),56)
    ]
end
