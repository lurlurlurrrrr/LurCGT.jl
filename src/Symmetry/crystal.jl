"""
    CrystalVertex{S}(word::Vector{Int})

A crystal vertex of symmetry `S`, represented by a left-to-right tensor word of
integer symbols. Each symbol determines its weight and local crystal operators.
The constructor takes ownership of `word`; never mutate it after hashing the
vertex or inserting it into a dictionary/set. Use nonmutating crystal operators
on published vertices and in-place operators only on owned scratch words.
"""
struct CrystalVertex{S<:NonabelianSymm}
    word::Vector{Int}
end

"""
    Crystal_ops{S}(f, e, ϕ, ϵ)

Local crystal operator data for symmetry `S`. `f[i]` and `e[i]` map integer
symbols under lowering and raising of color `i`; missing entries mean undefined
action. `ϕ[(i, symbol)]` and `ϵ[(i, symbol)]` store nonzero string lengths.
Zero is a valid symbol in B and G2, not an undefined-action sentinel.
"""
struct Crystal_ops{S<:NonabelianSymm}
    f::Vector{Dict{Int, Int}}
    e::Vector{Dict{Int, Int}}
    ϕ::Dict{Tuple{Int, Int}, Int}
    ϵ::Dict{Tuple{Int, Int}, Int}
end

const _CRYSTAL_WEIGHT_DATA = Dict{DataType, Tuple{Vector{Vector{Int}}, Dict{Int, Int}}}()
const _CRYSTAL_WEIGHT_DATA_LOCK = ReentrantLock()

"""
    _crystal_weight_data(::Type{S}) -> (weights, character_map)

Return cached defining-representation weights and character indices for `S`.
The cache owns these mutable results; callers must treat them as read-only.
Initialization is synchronized because crystal generation may run concurrently.
"""
function _crystal_weight_data(::Type{S}) where {S<:NonabelianSymm}
    lock(_CRYSTAL_WEIGHT_DATA_LOCK) do
        return get!(_CRYSTAL_WEIGHT_DATA, S) do
            weights = getsz_def_vec(S)::Vector{Vector{Int}}
            charmap = crystal_chars_map(S)::Dict{Int, Int}
            (weights, charmap)
        end
    end
end

"""
    mw_crystal_qlabel(S, qlabel) -> CrystalVertex{S}

Construct the highest tensor word in ascending fundamental-node order by
concatenating `mw_column(S, i, qlabel)` for each fundamental node `i`. The
symmetry method handles its label convention and returns all columns for `i`.
The returned word is owned.
"""
function mw_crystal_qlabel(::Type{S}, qlabel::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ}
    @assert NZ == nlops(S)
    @assert all(>=(0), qlabel)
    word = Int[]
    for i in 1:NZ
        append!(word, mw_column(S, i, qlabel))
    end
    return CrystalVertex{S}(word)
end

"""
    show(io, vertex::CrystalVertex)

Display the symmetry and tensor word of a vertex without reconstructing a tableau.
"""
function Base.show(io::IO, vertex::CrystalVertex{S}) where S<:NonabelianSymm
    print(io, "CrystalVertex{", S, "}(")
    show(io, vertex.word)
    print(io, ")")
end

"""
    get_crystal_ops(S) -> Crystal_ops{S}

Build local lowering/raising maps and their string lengths from the existing
symbol alphabet of `S`. No tensor-product vertices are enumerated.
"""
function get_crystal_ops(::Type{S}) where S <: NonabelianSymm
    char_list = charlist(S)
    f = get_fops_std(S)
    e = get_eops_std(f)
    ϕ = maximal_apply_num(char_list, f)
    ϵ = maximal_apply_num(char_list, e)
    return Crystal_ops{S}(f, e, ϕ, ϵ)
end

"""
    maximal_apply_num(char_list, ops) -> Dict{Tuple{Int,Int},Int}

Count successive local operations for each `(color, symbol)` in the finite
crystal alphabet. Store only nonzero lengths. `ops` must consist of finite strings.
"""
function maximal_apply_num(char_list::Vector{Int}, ops::Vector{Dict{Int, Int}})
    max_apply = Dict{Tuple{Int, Int}, Int}()
    for (i, op) in enumerate(ops)
        for char in char_list
            count = 0
            current = char
            while haskey(op, current)
                current = op[current]
                count += 1
            end
            if count > 0 max_apply[(i, char)] = count end
        end
    end
    return max_apply
end

"""
    get_eops_std(fops_std) -> Vector{Dict{Int,Int}}

Reverse the local lowering maps to obtain raising maps. Each color must be
injective on its defined domain, as required for disjoint crystal strings.
"""
function get_eops_std(fops_std::Vector{Dict{Int, Int}})
    eops_std = Vector{Dict{Int, Int}}()
    for fop in fops_std
        eop = Dict{Int, Int}()
        for (k, v) in fop
            @assert !haskey(eop, v) "Crystal lowering edges must have unique inverses"
            eop[v] = k
        end
        push!(eops_std, eop)
    end
    return eops_std
end

"""
    _crystal_signature(vertex, cops, i) -> (epsilon, phi, raising_index, lowering_index)

Evaluate the chapter 2 tensor rule without allocating signs. Each factor supplies
`phi` minus signs then `epsilon` plus signs; cancel `+-` pairs. Return the first
surviving plus and last surviving minus positions, or zero if absent. The input
word and operator data are never mutated.
"""
function _crystal_signature(vertex::CrystalVertex{S}, cops::Crystal_ops{S}, i::Int) where S<:NonabelianSymm
    @boundscheck checkbounds(cops.f, i)
    eps, ph, raising, lowering = 0, 0, 0, 0
    for k in eachindex(vertex.word)
        symbol = vertex.word[k]
        local_phi = get(cops.ϕ, (i, symbol), 0)
        local_eps = get(cops.ϵ, (i, symbol), 0)
        # New minus signs cancel the previous plus suffix, never earlier minuses.
        if local_phi > eps
            lowering = k
        end
        if local_phi >= eps
            raising = local_eps > 0 ? k : 0
        end
        ph += max(0, local_phi - eps)
        eps = local_eps + max(0, eps - local_phi)
    end
    return eps, ph, raising, lowering
end

"""
    epsilon(vertex::CrystalVertex, cops, i) -> Int

Return the maximal number of successive color-`i` raising operations on the
product vertex, using the chapter 2 tensor rule. Does not mutate or enumerate it.
"""
epsilon(vertex::CrystalVertex{S}, cops::Crystal_ops{S}, i::Int) where S<:NonabelianSymm =
    _crystal_signature(vertex, cops, i)[1]

"""
    phi(vertex::CrystalVertex, cops, i) -> Int

Return the maximal number of successive color-`i` lowering operations on the
product vertex, using the chapter 2 tensor rule. Does not mutate or enumerate it.
"""
phi(vertex::CrystalVertex{S}, cops::Crystal_ops{S}, i::Int) where S<:NonabelianSymm =
    _crystal_signature(vertex, cops, i)[2]

"""
    apply_lowering!(vertex::CrystalVertex, cops, i) -> Bool

Lower an owned scratch word at the rightmost surviving minus of its `i`-signature.
Return false without mutation if undefined. Never use on a hashed/published vertex.
"""
function apply_lowering!(vertex::CrystalVertex{S}, cops::Crystal_ops{S}, i::Int) where S<:NonabelianSymm
    k = _crystal_signature(vertex, cops, i)[4]
    iszero(k) && return false
    vertex.word[k] = cops.f[i][vertex.word[k]]
    return true
end

"""
    apply_lowering(vertex::CrystalVertex, cops, i) -> Union{CrystalVertex,Nothing}

Return a lowered vertex with a new owned word, or `nothing` if undefined. The
input word is never mutated, including after it has been hashed or stored as a key.
"""
function apply_lowering(vertex::CrystalVertex{S}, cops::Crystal_ops{S}, i::Int) where S<:NonabelianSymm
    k = _crystal_signature(vertex, cops, i)[4]
    iszero(k) && return nothing
    word = copy(vertex.word)
    word[k] = cops.f[i][word[k]]
    return CrystalVertex{S}(word)
end

"""
    apply_raising!(vertex::CrystalVertex, cops, i) -> Bool

Raise an owned scratch word at the leftmost surviving plus of its `i`-signature.
Return false without mutation if undefined. Never use on a hashed/published vertex.
"""
function apply_raising!(vertex::CrystalVertex{S}, cops::Crystal_ops{S}, i::Int) where S<:NonabelianSymm
    k = _crystal_signature(vertex, cops, i)[3]
    iszero(k) && return false
    vertex.word[k] = cops.e[i][vertex.word[k]]
    return true
end

"""
    apply_raising(vertex::CrystalVertex, cops, i) -> Union{CrystalVertex,Nothing}

Return a raised vertex with a new owned word, or `nothing` if undefined. The
input word is never mutated, including after it has been hashed or stored as a key.
"""
function apply_raising(vertex::CrystalVertex{S}, cops::Crystal_ops{S}, i::Int) where S<:NonabelianSymm
    k = _crystal_signature(vertex, cops, i)[3]
    iszero(k) && return nothing
    word = copy(vertex.word)
    word[k] = cops.e[i][word[k]]
    return CrystalVertex{S}(word)
end

"""
    get_weight(vertex::CrystalVertex{S}) -> NTuple{nzops(S),Int}

Sum the existing symbol weights in the package's z-coordinate convention.
The word is read-only; equal weights do not imply equal vertices.
"""
function get_weight(vertex::CrystalVertex{S}) where S <: NonabelianSymm
    cnt = get_count(vertex)
    NZ = nzops(S)
    # Supported alphabets store concrete integer z-coordinate lists.
    szs, chardict = _crystal_weight_data(S)
    @assert length(szs[1]) == NZ
    weight = zeros(Int, NZ)
    for char in keys(cnt)
        weight += cnt[char] * szs[chardict[char]]
    end
    return NTuple{NZ, Int}(weight)
end

"""
    det_getvector(vertex::CrystalVertex, i, cops, verticesNvecs) -> Bool

Return whether the lowering string below `vertex` contains an undiscovered word.
Preserve the numerical traversal's existing look-ahead rule. Only new words are
created; dictionary keys and their words remain unchanged.
"""
function det_getvector(vertex::CrystalVertex{S}, i::Int, cops::Crystal_ops{S},
    verticesNvecs::Dict{CrystalVertex{S}, Vector{RT}}) where {S<:NonabelianSymm, RT}
    while true
        vertex = apply_lowering(vertex, cops, i)
        vertex === nothing && break
        !haskey(verticesNvecs, vertex) && return true
    end
    return false
end

"""
    get_count(vertex::CrystalVertex) -> Dict{Int,Int}

Count occurrences of each symbol without changing the tensor word.
"""
function get_count(vertex::CrystalVertex)
    cnt = Dict{Int, Int}()
    for val in vertex.word
        cnt[val] = get(cnt, val, 0) + 1
    end
    return cnt
end

"""
    hash(vertex::CrystalVertex, h::UInt) -> UInt

Hash the word by contents using the existing sequential symbol hash. The caller
must keep the word unchanged after hashing; no hidden freezing or mutation occurs.
"""
function Base.hash(vertex::CrystalVertex, h::UInt)
    for symbol in vertex.word
        h = hash(symbol, h)
    end
    return h
end

"""
    ==(a::CrystalVertex{S}, b::CrystalVertex{S}) -> Bool

Compare tensor words by contents within the same symmetry. Different weights or
embeddings are not identified by this operation; it does not mutate either word.
"""
Base.:(==)(a::CrystalVertex{S}, b::CrystalVertex{S}) where S<:NonabelianSymm = a.word == b.word
