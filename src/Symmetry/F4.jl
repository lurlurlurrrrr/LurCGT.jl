# Basic metadata for the rank-four exceptional family.
isvalidsymm(::Type{F4}) = true
totxt(::Type{F4}) = "F4"
nlops(::Type{F4}) = 4
nzops(::Type{F4}) = 4

"""`_exceptional_cartan(F4)` uses long nodes 1,2 and short nodes 3,4."""
_exceptional_cartan(::Type{F4}) = [2 -1 0 0; -1 2 -1 0; 0 -2 2 -1; 0 0 -1 2]

"""`qlab2mwz(F4, q)` converts signed Dynkin coordinates to reversed simple-root weights."""
function qlab2mwz(::Type{F4}, q::NTuple{4,Int})
    rows = ((2,4,3,2),(4,8,6,3),(3,6,4,2),(2,3,2,1))
    return ntuple(i -> foldl(Base.Checked.checked_add,
        (Base.Checked.checked_mul(rows[i][j],q[j]) for j in 1:4); init=0),4)
end

"""`getqlabel(F4, z)` converts stored weights to signed Dynkin coordinates exactly."""
function getqlabel(::Type{F4}, z::NTuple{4,Int})
    A = _exceptional_cartan(F4); u = reverse(z)
    return ntuple(i -> foldl(Base.Checked.checked_add,
        (Base.Checked.checked_mul(A[i,j],u[j]) for j in 1:4); init=0),4)
end

"""`getdz(F4, i)` returns the positive root step at reversed slot 5-i."""
getdz(::Type{F4}, i::Int) = ntuple(j -> Int(j == 5-i),4)

"""`getdzs(F4)` returns the four integer simple-root coordinate steps."""
getdzs(::Type{F4}) = [collect(getdz(F4,i)) for i in 1:4]

"""`less_weight(F4, a, b)` compares stored weights in reversed lexicographic order."""
less_weight(::Type{F4}, a::Tuple, b::Tuple) = rev_less(a,b)

"""`getsl_triv(F4, RT)` returns empty trivial lowering blocks."""
getsl_triv(::Type{F4}, ::Type{RT}) where RT = ntuple(_ -> Dict{NTuple{4,Int},SparseMatrixCSC{RT}}(),4)

"""`getsz_triv(F4)` returns the one-dimensional zero-weight sector."""
getsz_triv(::Type{F4}) = Dict((0,0,0,0)=>(1,1))

"""`fundamental_recipes(F4)` returns the virtual seed and staged binary product recipes."""
fundamental_recipes(::Type{F4}) = FundamentalRecipe[
    ProductRecipe([1,0,0,0],([0,0,0,1],[0,0,0,1]),1,52),
    ProductRecipe([0,1,0,0],([0,0,1,0],[0,0,0,1]),1,1274),
    ProductRecipe([0,0,1,0],([0,0,0,1],[0,0,0,1]),1,273),
    VirtualRecipe([0,0,0,1],26,E6,[[2,2],[4,4],[3,5],[1,6]])]

"""`fundamental_qlabels(F4)` returns fundamental labels in node order."""
fundamental_qlabels(::Type{F4}) = [_recipe_label(F4,r.target) for r in fundamental_recipes(F4)]

"""`irrep_validation_rules(F4)` checks all adjoint copies in fundamental-times-dual products."""
irrep_validation_rules(::Type{F4}) = _irrep_validation_rules(F4,(1,0,0,0))

"""`crystal_basis_version(F4)` isolates the virtual seed and dependent numerical data."""
crystal_basis_version(::Type{F4}) = "crystal-seed-v1"

"""
    _complete_seed_blocks!(F4, lowering, metric, weights, f, local_index)

Complete the 26-dimensional short-root seed in its stored Demazure order.
The two zero vertices y3,y4 have Gram matrix [2 1; 1 2]. Its overlap requires
cross lowering coefficients from y4 under F3 and y3 under F4, although those
numerical actions are not crystal arrows. Reject other multiplicities.
"""
function _complete_seed_blocks!(::Type{F4}, lowering, metric, weights, f, local_index)
    zero_weight = (0,0,0,0)
    length(weights) == 26 && count(==(zero_weight),weights) == 2 ||
        error("Unexpected F4 direct seed multiplicities")
    all(w == zero_weight || size(G)==(1,1) for (w,G) in metric) ||
        error("Unexpected repeated nonzero F4 weight")
    zero_vertices = Int[]
    for i in (3,4)
        root = Tuple(getdz(F4,i))
        source = only(findall(==(root),weights))
        vertex = f[i][source]
        weights[vertex] == zero_weight || error("Missing F4 zero string for node $i")
        push!(zero_vertices,vertex)
    end
    length(unique(zero_vertices)) == 2 || error("F4 zero vertices are not distinct")
    G = metric[zero_weight]
    G == [2 0; 0 2] || error("Unexpected F4 zero-string norms")
    G[1,2] = G[2,1] = 1
    for (a,i) in enumerate((3,4))
        other = zero_vertices[3-a]
        lowering[i][zero_weight][1,local_index[other]] = 1
    end
    return nothing
end
