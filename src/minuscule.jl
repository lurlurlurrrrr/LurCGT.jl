"""`_unit_label(S, k)` returns the rank-length unit Dynkin label at k."""
_unit_label(::Type{S}, k::Int) where S = ntuple(i -> Int(i == k), nzops(S))

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

"""`_has_minuscule_recipes(S)` identifies the families supported by Part 3a."""
_has_minuscule_recipes(::Type) = false

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

"""
    _component_vertices(S, highest) -> Vector{CrystalVertex{S}}

Enumerate a component in the fixed Demazure first-discovery order. Frozen source
counts and stopping at known vertices match numerical basis construction.
"""
function _component_vertices(::Type{S}, highest::CrystalVertex{S}) where S
    vertices = [highest]; seen = Set(vertices); cops = get_crystal_ops(S)
    for i in crystal_lowering_word(S)
        nsource = length(vertices)
        for k in 1:nsource
            v = vertices[k]
            while true
                v = apply_lowering(v,cops,i)
                v === nothing && break
                v in seen && break
                push!(seen,v); push!(vertices,v)
            end
        end
    end
    return vertices
end

const _FUNDAMENTAL_WORDS = Dict{DataType,Dict{Int,Vector{Int}}}()
const _FUNDAMENTAL_WORDS_LOCK = ReentrantLock()

"""
    _fundamental_word(S, k) -> Vector{Int}

Return an owned copy of the deterministic fundamental highest word. Cache only
integer words, never ambient numerical representations or mutable vertices.
"""
function _fundamental_word(::Type{S}, k::Int) where S
    lock(_FUNDAMENTAL_WORDS_LOCK) do
        words = get!(() -> Dict{Int,Vector{Int}}(), _FUNDAMENTAL_WORDS,S)
        word = get!(words,k) do
            _build_fundamental_word(S,k)
        end
        return copy(word)
    end
end

"""
    _build_fundamental_word(S, k) -> Vector{Int}

Return a fundamental highest word using its fixed binary recipe. Minuscule
fundamentals are atomic symbols; other fundamentals retain source words.
"""
function _build_fundamental_word(::Type{S}, k::Int) where S
    if k in minuscule_nodes(S)
        return [_minuscule_alphabet(S)[3][k]]
    end
    a,b = fundamental_sources(S,k)
    left = _component_vertices(S, CrystalVertex{S}(_fundamental_word(S,a)))
    right = _component_vertices(S, CrystalVertex{S}(_fundamental_word(S,b)))
    target = qlab2mwz(S,_unit_label(S,k)); cops = get_crystal_ops(S)
    for x in left, y in right
        get_weight(x) .+ get_weight(y) == target || continue
        v = CrystalVertex{S}(vcat(x.word,y.word))
        all(i -> epsilon(v,cops,i)==0, 1:nlops(S)) && return v.word
    end
    error("Fundamental recipe has no highest vertex for $S node $k")
end

"""`_exact_divide(x, d)` divides integer storage values, rejecting nonzero remainder."""
function _exact_divide(x::Integer, d::Integer)
    q,r = divrem(x,d)
    iszero(r) || error("Nonintegral representation coefficient: $x / $d")
    return q
end

"""
    _orthogonal_remainder(v, basis, metric) -> Vector

Project an integer candidate off mutually metric-orthogonal stored columns.
Clear denominators by integer cross multiplication and primitive rescaling;
no rational vector or nonexact division is introduced.
"""
function _orthogonal_remainder(v::Vector{RT}, basis::Vector{Vector{RT}}, metric::Matrix{RT}) where RT
    v = copy(v)
    for u in basis
        overlap = dot(u,metric*v)
        iszero(overlap) && continue
        norm = dot(u,metric*u)
        v = norm*v - overlap*u
        iszero(v) && return v
        v = div.(v,gcd(v))
    end
    return v
end

"""
    get_irep_frommw(S, Sl, Sz, metric1, metric2, mw, weight_info; copy_index=1, candidates=nothing)

Extract a target component from an ordered product, constructing all preceding
multiplets first. Select copies by exact projection of standard basis vectors
in index order. Optional `candidates` collects nonfundamental normalized irreps
without saving them. Reject absent targets and unavailable copies.
"""
function get_irep_frommw(::Type{S},
        Sl::NTuple{NL,Dict{NTuple{NZ,Int},SparseMatrixCSC{RT}}},
        Sz::Dict{NTuple{NZ,Int},Tuple{Int,Int}},
        metric1::Dict{NTuple{NZ,Int},Matrix{RT}},
        metric2::Dict{NTuple{NZ,Int},Matrix{RT}},
        mw::Union{Nothing,NTuple{NZ,Int}},
        winfo::Dict{NTuple{NZ,Int},WeightInfo{NZ}};
        copy_index::Int=1, candidates=nothing, verbose=VERBOSE[]) where {S,NL,NZ,RT<:Integer}
    copy_index > 0 || throw(ArgumentError("Copy index must be positive"))
    weights = sorted_zvals(Sz)
    target = isnothing(mw) ? first(weights) : mw
    haskey(Sz,target) || throw(ArgumentError("Target weight absent from product"))
    # The maximal Cartan component has a one-dimensional highest sector.
    if target == first(weights) && copy_index == 1
        rep, _ = _irep_component(S,Sl,Sz,metric1,metric2,target,winfo; verbose)
        if candidates !== nothing && sum(rep.qlabel) > 1; push!(candidates,rep); end
        return rep
    end
    metrics = Dict{NTuple{NZ,Int},Matrix{RT}}()
    found = Dict(w => Vector{RT}[] for w in weights)
    target_copies = 0
    for w in weights
        rev_less(w,target) && break
        d = Sz[w][2]-Sz[w][1]+1
        if !haskey(metrics,w)
            metrics[w] = get_innerprodmat(S,RT,w,Matrix{RT}(I,d,d),d,metric1,metric2,winfo; verbose)
        end
        for index in 1:d
            length(found[w]) == d && break
            v = zeros(RT,d); v[index] = 1
            v = _orthogonal_remainder(v,found[w],metrics[w])
            iszero(v) && continue
            all(>=(0),getqlabel(S,w)) || error("Residual product state is not highest weight")
            rep, bases = _irep_component(S,Sl,Sz,metric1,metric2,w,winfo; mwstate=v, verbose)
            if candidates !== nothing && sum(rep.qlabel)>1; push!(candidates,rep); end
            if w == target
                target_copies += 1
                target_copies == copy_index && return rep
            end
            for z in sorted_zvals(bases)
                B = bases[z]; zd = size(B,1)
                if !haskey(metrics,z)
                    metrics[z] = get_innerprodmat(S,RT,z,Matrix{RT}(I,zd,zd),zd,metric1,metric2,winfo; verbose)
                end
                for j in axes(B,2)
                    u = _orthogonal_remainder(B[:,j],found[z],metrics[z])
                    iszero(u) && error("Dependent multiplet column")
                    push!(found[z],u)
                end
            end
        end
    end
    throw(ArgumentError("Target irrep or requested copy absent from product"))
end

"""
    construct_fundamental(S, RT, k; candidates=nothing) -> Irep

Execute the requested hard-coded fundamental recipe, recursively loading only
its required sources. Optional candidates receives incidental nonfundamental
components; the caller saves the returned fundamental.
"""
function construct_fundamental(::Type{S}, ::Type{RT}, k::Int; candidates=nothing, verbose=VERBOSE[]) where {S,RT<:Integer}
    k in fundamental_sequence(S) || throw(ArgumentError("Unknown fundamental node"))
    k in minuscule_nodes(S) && return get_minuscule_irep(S,RT,k)
    a,b = fundamental_sources(S,k)
    x = getNsave_irep(S,RT,_unit_label(S,a); verbose)
    y = getNsave_irep(S,RT,_unit_label(S,b); verbose)
    Sl,Sz,_,wi = tensor_prod_reps(x,y; verbose)
    return get_irep_frommw(S,Sl,Sz,x.innerprod,y.innerprod,qlab2mwz(S,_unit_label(S,k)),wi; candidates,verbose)
end

"""
    same_irep_fields(a, b) -> Bool

Compare every mathematical Irep field exactly, including primitive inverse
metric factors. Ignore only the allocation-dependent `size_byte` field.
"""
function same_irep_fields(a::Irep, b::Irep)
    typeof(a) === typeof(b) || return false
    return a.Sl == b.Sl && a.Sz == b.Sz && a.innerprod == b.innerprod &&
        a.inv_innerprod == b.inv_innerprod && a.qlabel == b.qlabel && a.dimension == b.dimension
end

"""
    validate_fundamental_ireps(S, RT=BigInt; verbose=0) -> Int

Construct and save fundamentals in hard-coded dependency order, collecting all
incidental nonfundamental copies outside the cache. After every fundamental is
saved, compare each candidate with its canonical construction. Return the count
of validated copies; throw on the first unequal field set. Run with a fresh
database so canonical references cannot come from a previous backend.
"""
function validate_fundamental_ireps(::Type{S}, ::Type{RT}=BigInt; verbose=0) where {S,RT<:Integer}
    initialize_symmetry(S,RT;verbose)
    candidates = Irep{S,nlops(S),nzops(S),RT}[]
    for k in fundamental_sequence(S)
        rep = construct_fundamental(S,RT,k; candidates,verbose)
        save_irep_sqlite(rep;verbose)
    end
    for candidate in candidates
        canonical = getNsave_irep(S,RT,candidate.qlabel;verbose)
        same_irep_fields(candidate,canonical) || error("Canonical Irep mismatch for $S $(candidate.qlabel)")
    end
    return length(candidates)
end

"""
    _fundamental_validation_rules(S, adjoint) -> Vector{Tuple{label,label,label}}

Build the shared adjoint-in-fundamental-times-dual checks in recipe order.
Family methods supply the adjoint label and append their special cases.
"""
function _fundamental_validation_rules(::Type{S}, adjoint::NTuple{R,Int}) where {S,R}
    rules = Tuple{NTuple{R,Int},NTuple{R,Int},NTuple{R,Int}}[]
    for k in fundamental_sequence(S)
        left = _unit_label(S,k); right = get_dualq(S,left)
        left > right && continue
        push!(rules,(adjoint,left,right))
    end
    return rules
end

"""
    validate_fundamental_products(S, RT=BigInt; rules=fundamental_validation_rules(S)) -> Int

Execute the supplied finite cross-product rules after fundamental construction.
Compare each independently extracted target with its canonical cached irrep;
return the number of passing rules. Large exceptional rules are opt-in work.
"""
function validate_fundamental_products(::Type{S}, ::Type{RT}=BigInt;
        rules=fundamental_validation_rules(S), verbose=0) where {S,RT<:Integer}
    for (target,left,right) in rules
        x = getNsave_irep(S,RT,left;verbose); y = getNsave_irep(S,RT,right;verbose)
        Sl,Sz,_,wi = tensor_prod_reps(x,y;verbose)
        result = get_irep_frommw(S,Sl,Sz,x.innerprod,y.innerprod,qlab2mwz(S,target),wi;verbose)
        canonical = getNsave_irep(S,RT,target;verbose)
        same_irep_fields(result,canonical) || error("Cross-product mismatch for $S $target from $left × $right")
    end
    return length(rules)
end

