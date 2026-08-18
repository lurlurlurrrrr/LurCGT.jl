"""
Outer-multiplicity permutation needed to conjugate a CGT.

# Fields

- `perm`: one-based permutation of flattened outer-multiplicity states.
- `upsp`: sorted upper-leg q-label tuple defining that OM basis.
- `size_byte`: cached memory footprint for LRU eviction.
"""
struct Conjperm{S<:NonabelianSymm, U, NZ}
    perm::Vector{Int}
    upsp::NTuple{U, NTuple{NZ, Int}}
    size_byte::Int

    function Conjperm{S, U, NZ}(perm, upsp, size_byte::Int=0) where {S<:NonabelianSymm, U, NZ}
        if size_byte == 0
            obj = new{S, U, NZ}(perm, upsp, 0)
            size_byte = Base.summarysize(obj)
        end
        new{S, U, NZ}(perm, upsp, size_byte)
    end
end

"""
    getNsave_Conjperm(::Type{S}, upsp; save=true)

Load or compute the OM-basis permutation needed for CGT conjugation.

`upsp` is the sorted tuple of upper/input qlabels defining the CGT OM basis.
Abelian symmetries return `nothing` because their OM space is trivial. For
non-Abelian symmetries, leading trivial qlabels are removed before the
standardized cache lookup. `save` controls SQLite persistence for generated
permutations.
"""
getNsave_Conjperm(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}};
    save=true) where {S<:AbelianSymm, U, NZ} = nothing

function getNsave_Conjperm(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}};
    save=true) where {S<:NonabelianSymm, U, NZ}

    @assert issorted(upsp)
    upsp_, _ = remove_zeros(S, upsp, ntuple(identity, Val(U)))
    getNsave_Conjperm_std(S, upsp_; save=save)
end

"""
    getNsave_Conjperm_std(::Type{S}, upsp; save=true) -> Conjperm

Load or compute a standardized non-Abelian conjugation permutation.

`upsp` must already be sorted and have qlabel width `nzops(S)`. The method
checks SQLite first and computes the permutation from the CGT OM metadata on a
cache miss.
"""
function getNsave_Conjperm_std(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}};
    save=true) where {S<:NonabelianSymm, U, NZ}

    @assert NZ == nzops(S)
    @assert issorted(upsp)

    loaded = load_Conjperm_sqlite(S, upsp)
    !isnothing(loaded) && return loaded
    return computeNsave_Conjperm_std(S, upsp; save)
end

"""
    computeNsave_Conjperm_std(::Type{S}, upsp; save=true) -> Conjperm

Compute the conjugation permutation for one non-Abelian CGT OM basis.

`upsp` defines the input and output qlabel tuple of the self-conjugation
problem. The permutation is obtained from `get_conj_perm(get_CGTom(S, upsp,
upsp))`, wrapped in `Conjperm`, and optionally saved to SQLite.
"""
function computeNsave_Conjperm_std(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}};
    save=true) where {S<:NonabelianSymm, U, NZ}

    cgt_oms = get_CGTom(S, upsp, upsp)
    obj = Conjperm{S, U, NZ}(get_conj_perm(cgt_oms), upsp)
    save && save_Conjperm_sqlite(S, obj)
    return obj
end
