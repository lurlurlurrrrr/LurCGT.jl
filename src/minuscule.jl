"""`minuscule_nodes(S)` derives seed nodes from the sole construction catalogue."""
minuscule_nodes(::Type{S}) where S = findall(recipe -> recipe isa MinusculeRecipe, fundamental_recipes(S))

"""
    minuscule_weights(S, k) -> Vector{NTuple}

Enumerate a minuscule orbit in deterministic first-discovery lowering order.
The highest weight is first; each weight occurs exactly once.
"""
function minuscule_weights(::Type{S}, k::Int) where S
    k in minuscule_nodes(S) || throw(ArgumentError("Node $k is not minuscule for $S"))
    highest = qlab2mwz(S, _unit_label(S, k))::NTuple{nzops(S),Int}
    weights = [highest]
    seen = Set(weights)
    cursor = 1
    while cursor <= length(weights)
        w = weights[cursor]
        a = getqlabel(S, w)
        for i in 1:nlops(S)
            a[i] == 1 || continue
            v = w .- getdz(S, i)
            if !(v in seen)
                push!(seen, v); push!(weights, v)
            end
        end
        cursor += 1
    end
    return weights
end

"""
    get_minuscule_irep(S, RT, k) -> Irep

Construct fundamental k directly from its minuscule orbit. Root edges have
unit coefficients and every one-dimensional weight space has unit metric.
No database entries are read or written.
"""
function get_minuscule_irep(::Type{S}, ::Type{RT}, k::Int) where {S, RT<:Integer}
    ws = minuscule_weights(S, k)::Vector{NTuple{nzops(S),Int}}
    sort!(ws; lt=rev_less, rev=true)
    NZ = nzops(S)
    Sz = Dict(w => (i,i) for (i,w) in enumerate(ws))
    Sl = ntuple(_ -> Dict{NTuple{NZ,Int},SparseMatrixCSC{RT}}(), nlops(S))
    for w in ws, i in 1:nlops(S)
        getqlabel(S, w)[i] == 1 || continue
        Sl[i][w] = sparse(reshape(RT[1], 1, 1))
    end
    metric = Dict(w => ones(RT, 1, 1) for w in ws)
    return Irep(S, RT, Sl, Sz, metric, get_inv_innerprod(S, metric), length(ws))
end

"""`_minuscule_symbol_start(S)` starts a fresh integer alphabet at one by default."""
_minuscule_symbol_start(::Type) = 1

"""`_preserve_vector_symbols(S, k)` defaults to an independently numbered seed."""
_preserve_vector_symbols(::Type, k::Int) = false

"""`_encode_minuscule_symbol(S, w, next)` allocates the next symbol with checked increment."""
_encode_minuscule_symbol(::Type, w::Tuple, next::Int) = (next, Base.Checked.checked_add(next,1))

"""
    _minuscule_alphabet(S) -> (weights, symbol_indices, highest_symbols, operators)

Build disjoint deterministic integer alphabets for all minuscule seeds.
Existing defining symbols are preserved; other seeds occupy successive positive
ranges. All tables are owned, and weights are stored in local z coordinates.
"""
function _minuscule_alphabet(::Type{S}) where S
    weights = Vector{Int}[]
    indices = Dict{Int,Int}()
    highest = Dict{Int,Int}()
    f = [Dict{Int,Int}() for _ in 1:nlops(S)]
    nextsymbol = _minuscule_symbol_start(S)
    for k in minuscule_nodes(S)
        ws = minuscule_weights(S, k)
        symbols = Int[]
        preserve = _preserve_vector_symbols(S,k)
        oldweights = preserve ? getsz_def_vec(S) : Vector{Int}[]
        oldmap = preserve ? crystal_chars_map(S) : Dict{Int,Int}()
        for w in ws
            if preserve
                ix = findfirst(==(collect(w)), oldweights)
                ix === nothing && error("Minuscule vector weight missing from existing alphabet")
                symbol = only(c for (c,j) in oldmap if j == ix)
            else
                symbol, nextsymbol = _encode_minuscule_symbol(S,w,nextsymbol)
            end
            push!(symbols, symbol)
            push!(weights, collect(w)); indices[symbol] = length(weights)
        end
        highest[k] = first(symbols)
        lookup = Dict(w => symbols[j] for (j,w) in enumerate(ws))
        for (j,w) in enumerate(ws), i in 1:nlops(S)
            getqlabel(S,w)[i] == 1 || continue
            f[i][symbols[j]] = lookup[w .- getdz(S,i)]
        end
    end
    e = get_eops_std(f)
    chars = sort!(collect(keys(indices)))
    return weights, indices, highest, Crystal_ops{S}(f,e,maximal_apply_num(chars,f),maximal_apply_num(chars,e))
end

