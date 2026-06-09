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

function getNsave_Conjperm_std(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}};
    save=true) where {S<:NonabelianSymm, U, NZ}

    @assert NZ == nzops(S)
    @assert issorted(upsp)

    loaded = load_Conjperm_sqlite(S, upsp)
    !isnothing(loaded) && return loaded

    cgt_oms = get_CGTom(S, upsp, upsp)
    obj = Conjperm{S, U, NZ}(get_conj_perm(cgt_oms), upsp)
    save && save_Conjperm_sqlite(S, obj)
    return obj
end
