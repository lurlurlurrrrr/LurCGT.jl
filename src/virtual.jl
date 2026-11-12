"""
    _virtual_ambient_label(recipe, q, A) -> NTuple

Embed signed target Dynkin coordinates by counting ambient indices in each
operator word. Metadata validation must precede this conversion.
"""
function _virtual_ambient_label(recipe::VirtualRecipe, q, ::Type{A}) where A
    return ntuple(j -> sum(q[i]*count(==(j),recipe.operator_words[i])
        for i in eachindex(q)),nzops(A))
end

"""
    _validate_recipe(S, recipe::VirtualRecipe)

Check the orbit partition, contiguous equal repetitions, commuting ambient
nodes, and compatibility of the induced root map. No crystals are constructed.
"""
function _validate_recipe(::Type{S}, recipe::VirtualRecipe) where S
    invoke(_validate_recipe,Tuple{Type{S},FundamentalRecipe},S,recipe)
    A = recipe.ambient
    A <: NonabelianSymm && applicable(isvalidsymm,A) && isvalidsymm(A) &&
        !isempty(fundamental_recipes(A)) ||
        throw(ArgumentError("Invalid virtual ambient symmetry $A"))
    length(recipe.operator_words) == nlops(S) || throw(ArgumentError("Wrong virtual operator count"))
    used = Int[]
    for word in recipe.operator_words
        isempty(word) && throw(ArgumentError("Empty virtual operator word"))
        all(j -> 1 <= j <= nlops(A),word) || throw(ArgumentError("Invalid ambient node"))
        nodes = unique(word); gamma = count(==(first(word)),word)
        for j in nodes
            count(==(j),word) == gamma || throw(ArgumentError("Unequal virtual repetitions"))
            last(findall(==(j),word))-first(findall(==(j),word))+1 == gamma ||
                throw(ArgumentError("Virtual repetitions must be contiguous"))
            for k in nodes
                j == k && continue
                getqlabel(A,Tuple(getdz(A,j)))[k] == 0 || throw(ArgumentError("Noncommuting virtual nodes"))
            end
        end
        append!(used,nodes)
    end
    sort(used) == collect(1:nlops(A)) || throw(ArgumentError("Virtual orbits must partition ambient nodes"))
    for i in 1:nlops(S)
        embedded = _virtual_ambient_label(recipe,getqlabel(S,Tuple(getdz(S,i))),A)
        ambient_root = ntuple(k -> sum(getqlabel(A,Tuple(getdz(A,j)))[k]
            for j in recipe.operator_words[i]),nzops(A))
        embedded == ambient_root || throw(ArgumentError("Incompatible virtual root map for $S node $i"))
    end
    return nothing
end

"""
    _virtual_weight(S, recipe, A, ambient_weight) -> NTuple

Recover signed target Dynkin coordinates by exact orbit-wise division, then
convert to target stored weights. Reject nonaligned ambient weights.
"""
function _virtual_weight(::Type{S}, recipe::VirtualRecipe, ::Type{A}, w) where {S,A}
    q = getqlabel(A,w)
    target = ntuple(nzops(S)) do i
        word = recipe.operator_words[i]; j = first(word); gamma = count(==(j),word)
        all(k -> q[k] == q[j],word) && rem(q[j],gamma) == 0 ||
            error("Unaligned virtual weight for $S $(Tuple(recipe.target)): $w")
        div(q[j],gamma)
    end
    return qlab2mwz(S,target)
end

"""
    _virtual_step(vertex, cops, word, raising=false)

Execute a grouped crystal operation on owned scratch storage. Raising reverses
the ambient word. Undefined constituent steps return nothing without mutation
of the input, even if earlier steps succeeded.
"""
function _virtual_step(vertex::CrystalVertex{A}, cops::Crystal_ops{A}, word::Vector{Int}, raising::Bool=false) where A
    result = CrystalVertex{A}(copy(vertex.word))
    for j in (raising ? Iterators.reverse(word) : word)
        success = raising ? apply_raising!(result,cops,j) : apply_lowering!(result,cops,j)
        success || return nothing
    end
    return result
end

"""
    _virtual_seed_data(S, recipe, A)

Enumerate an ambient virtual component in target Demazure order. Return its
ambient vertices, target weights, and local lowering maps.
Distinct equal-weight vertices remain
separate. All crystal alignment and dimension checks run before returning.
"""
function _virtual_seed_data(::Type{S}, recipe::VirtualRecipe, ::Type{A}) where {S,A}
    _validate_recipe(S,recipe)
    A === recipe.ambient || throw(ArgumentError("Wrong virtual ambient dispatch"))
    q = _virtual_ambient_label(recipe,recipe.target,A)::NTuple{nzops(A),Int}
    highest = mw_crystal_qlabel(A,q); cops = get_crystal_ops(A)
    vertices = [highest]; index = Dict(highest=>1)
    for i in crystal_lowering_word(S)
        nsource = length(vertices)
        for k in 1:nsource
            v = vertices[k]
            while true
                next = _virtual_step(v,cops,recipe.operator_words[i])
                next === nothing && break
                haskey(index,next) && break
                push!(vertices,next); index[next] = length(vertices)
                v = next
            end
        end
    end
    length(vertices) == recipe.dimension || error("Virtual dimension mismatch for $S $(Tuple(recipe.target)): $(length(vertices))")
    weights = [_virtual_weight(S,recipe,A,get_weight(v)) for v in vertices]
    weights[1] == qlab2mwz(S,_recipe_label(S,recipe.target)) || error("Wrong virtual highest weight")
    f = [Dict{Int,Int}() for _ in 1:nlops(S)]
    for (k,v) in enumerate(vertices), i in 1:nlops(S)
        word = recipe.operator_words[i]; gamma = count(==(first(word)),word)
        eps = epsilon(v,cops,first(word)); ph = phi(v,cops,first(word))
        all(j -> epsilon(v,cops,j)==eps && phi(v,cops,j)==ph,word) &&
            rem(eps,gamma)==0 && rem(ph,gamma)==0 || error("Unaligned virtual string for $S at vertex $k node $i")
        for raising in (false,true)
            next = _virtual_step(v,cops,word,raising)
            next === nothing && continue
            haskey(index,next) || error("Virtual crystal not closed for $S at vertex $k node $i")
            _virtual_step(next,cops,word,!raising) == v || error("Virtual operators not inverse for $S")
            if !raising
                f[i][k] = index[next]
                weights[index[next]] == weights[k] .- getdz(S,i) || error("Incorrect virtual lowering weight for $S")
            end
        end
    end
    e = get_eops_std(f); chars = collect(eachindex(vertices))
    epsilon_table = maximal_apply_num(chars,e); phi_table = maximal_apply_num(chars,f)
    for (k,v) in enumerate(vertices), i in 1:nlops(S)
        word = recipe.operator_words[i]; j = first(word); gamma = count(==(j),word)
        get(epsilon_table,(i,k),0) == div(epsilon(v,cops,j),gamma) &&
            get(phi_table,(i,k),0) == div(phi(v,cops,j),gamma) || error("Virtual string length mismatch for $S")
    end
    return vertices,weights,f
end

"""
    construct_recipe(S, RT, recipe::VirtualRecipe; candidates=nothing, verbose)

Construct an exact seed from target crystal strings and symmetry-specific
metric blocks. Ambient data supplies only the crystal; no ambient numerical
irrep is loaded. Integer storage and conflicting norm assignments are checked.
"""
function construct_recipe(::Type{S}, ::Type{RT}, recipe::VirtualRecipe;
        candidates=nothing, verbose=VERBOSE[]) where {S,RT<:Integer}
    weights,f = _seed_component(S,0,recipe)
    return _irep_from_seed_crystal(S,RT,weights,f)
end

"""
    _seed_component(S, k, recipe::VirtualRecipe) -> (weights, lowering_maps)

Return virtual seed data for shared alphabet assembly. Ambient type dispatch
occurs at the traversal boundary; equal-weight vertices retain distinct indices.
"""
function _seed_component(::Type{S}, k::Int, recipe::VirtualRecipe) where S
    _,weights,f = _virtual_seed_data(S,recipe,recipe.ambient)
    return weights,f
end
