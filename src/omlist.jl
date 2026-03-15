# S: Non-abelian symmetry type, N: number of incoming spaces
# NZ: length of Dynkin label, M1=N-1, M2=N-2
struct OMList{S, N, NZ, M1, M2} 
    # List of incoming spaces. Tuple of Dynkin labels. Sorted
    incom_spaces::NTuple{N, NTuple{NZ, Int}}
    out_space::NTuple{NZ, Int}
    # [(e1, e2, e3), ...] where e1, e2, e3 are intermediate spaces
    # The vector is sorted by dictionary order
    interm_spaces::Vector{NTuple{M2, NTuple{NZ, Int}}}
    # [(u1, u2, u3, u4), ...] where u1, ..., u4 are outer multiplicities of CG3s
    cg3_oms::Vector{NTuple{M1, Int}}
    # [1, a1, a2, a3, ...], 1 ~ a1-1th entry of outer multiplicity have
    # intermediate spaces (e1, e2, e3), and so on
    cumul::Vector{Int}
    # Total outer multiplicity
    totalOM::Int
    # Cached memory footprint in bytes, computed once at construction for LRU cache eviction by size
    size_byte::Int
    
    function OMList{S, N, NZ, M1, M2}(incom_spaces, out_space, interm_spaces, cg3_oms, cumul, totalOM, size_byte::Int=0) where {S, N, NZ, M1, M2}
        if size_byte == 0
            obj = new{S, N, NZ, M1, M2}(incom_spaces, out_space, interm_spaces, cg3_oms, cumul, totalOM, 0)
            size_byte = Base.summarysize(obj)
        end
        new{S, N, NZ, M1, M2}(incom_spaces, out_space, interm_spaces, cg3_oms, cumul, totalOM, size_byte)
    end
end

struct ValidOuts{S, N, NZ}
    # Incoming spaces are sorted from lowest to highest
    incom_spaces::NTuple{N, NTuple{NZ, Int}}
    # Dynkin labels of possible outgoing spaces. This is sorted.
    out_spaces::Vector{NTuple{NZ, Int}}
    # Cached memory footprint in bytes, computed once at construction for LRU cache eviction by size
    size_byte::Int
    
    function ValidOuts{S, N, NZ}(incom_spaces, out_spaces, size_byte::Int=0) where {S, N, NZ}
        if size_byte == 0
            obj = new{S, N, NZ}(incom_spaces, out_spaces, 0)
            size_byte = Base.summarysize(obj)
        end
        new{S, N, NZ}(incom_spaces, out_spaces, size_byte)
    end
end

function OMList(::Type{S},  
    incom_spaces::NTuple{N, NTuple{NZ, Int}},
    out_space, 
    interm_spaces, 
    cg3_oms) where {S<:NonabelianSymm, N, NZ}
    M1, M2 = N-1, max(N-2, 0)
    @assert length(out_space) == NZ
    @assert all(length.(interm_spaces) .== M2)
    @assert all(length.(cg3_oms) .== M1)
    cumul = Int[1]
    for t in cg3_oms
        push!(cumul, cumul[end] + prod(t))
    end
    totalOM = cumul[end] - 1
    return OMList{S, N, NZ, M1, M2}(incom_spaces, out_space, 
    interm_spaces, cg3_oms, cumul, totalOM)
end

function ValidOuts(::Type{S}, 
    incom_spaces::NTuple{N, NTuple{NZ, Int}},
    out_spaces::Vector{NTuple{NZ, Int}}) where {S<:NonabelianSymm, N, NZ}
    @assert NZ == nzops(S)
    return ValidOuts{S, N, NZ}(incom_spaces, out_spaces)
end


function getNsave_omlist(::Type{S},
        incom_spaces::NTuple{N, NTuple{NZ, Int}},
        out_space::NTuple{NZ, Int}) where {S<:NonabelianSymm, N, NZ}
    @assert NZ == nzops(S)
    @assert collect(incom_spaces) == sort(collect(incom_spaces))
    
    # If the file exists, just load it
    loaded = load_omlist_sqlite(S, incom_spaces, out_space)
    if !isnothing(loaded) return loaded end

    if N >= 2
        # If not, first load the vout file
        vo_obj = getNsave_validout(S, incom_spaces)
        # If there out_space is not possible from incoming spaces, return nothing
        if out_space ∉ vo_obj.out_spaces return nothing end
    end
    
    # If N==1, input and output spaces are the same
    if N == 1
        # If out_space is impossible, return nothing
        if incom_spaces[1] != out_space return nothing end
        return OMList(S, incom_spaces, out_space, Tuple{}[()], [()])
    elseif N == 2
        # If N == 2, this is just CG3 case.
        q1_, q2_ = incom_spaces
        # TODO: Always BigInt?
        omlookup(S, BigInt, q1_, q2_; verbose=0)
        loaded = load_omlist_sqlite(S, incom_spaces, out_space)
        @assert !isnothing(loaded); return loaded
    else
        # For N > 2, generate from N-1 rank OMList
        incom_spaces_prev = incom_spaces[1:(N-1)]
        last_incoming = incom_spaces[N]
        validout_prev = load_validout_sqlite(S, incom_spaces_prev)
        # Since it is created when getNsave_validout function call
        @assert !isnothing(validout_prev)
        # Initialize the field of OMList
        interm_spaces = Vector{NTuple{N-2, NTuple{NZ, Int}}}()
        cg3_oms = Vector{NTuple{N-1, Int}}()

        # validout_prev.out_spaces are all possible outgoing spaces
        # This is sorted
        for out_space_prev in validout_prev.out_spaces
            q1_, q2_ = minmax(out_space_prev, last_incoming)
            lastcg3_om = getNsave_validout(S, (q1_, q2_))
            idx = searchsortedfirst(lastcg3_om.out_spaces, out_space)

            # If out_space is possible outcome
            if idx <= length(lastcg3_om.out_spaces) && lastcg3_om.out_spaces[idx] == out_space
                # Load the low-rank OMList
                omlist_prev = getNsave_omlist(S, incom_spaces_prev, out_space_prev)
                lastcg3_omlist = getNsave_omlist(S, (q1_, q2_), out_space)

                interm_spaces_appended = N==3 ?  [(out_space_prev,)] : 
                map(x -> (out_space_prev, x...), omlist_prev.interm_spaces)
                append!(interm_spaces, interm_spaces_appended)
                cg3_oms_appended =
                    map(x -> (lastcg3_omlist.totalOM, x...), omlist_prev.cg3_oms)
                append!(cg3_oms, cg3_oms_appended)
            end
        end
        final_omlist = OMList(S, Tuple(incom_spaces), out_space,
            interm_spaces, cg3_oms)
        save_omlist_sqlite(S, final_omlist)
        return final_omlist
    end
end

function getNsave_validout(::Type{S},
        incom_spaces::NTuple{N, NTuple{NZ, Int}}) where {S<:NonabelianSymm, N, NZ}

    # Sanity checks
    if N == 1 return ValidOuts(S, incom_spaces, [incom_spaces[1]]) end
    @assert N >= 2; @assert NZ == nzops(S)
    @assert collect(incom_spaces) == sort(collect(incom_spaces))

    # If it exists, load and return
    loaded = load_validout_sqlite(S, incom_spaces)
    if !isnothing(loaded) return loaded end

    # If not, generate the new validout object
    if N == 2
        # If N == 2, this is just CG3 case.
        q1, q2 = incom_spaces
        # TODO: Always BigInt?
        omlookup(S, BigInt, q1, q2; verbose=0)
        loaded = load_validout_sqlite(S, incom_spaces)
        @assert !isnothing(loaded); return loaded
    end
    
    # If not, generate it from validout object which has N-1 incoming spaces
    incom_spaces_prev = incom_spaces[1:(N-1)]
    last_incoming = incom_spaces[N]
    vo_obj_prev = getNsave_validout(S, incom_spaces_prev)
    
    out_space_set = Set{NTuple{NZ, Int}}()
    for out_space_prev in vo_obj_prev.out_spaces
        q1_, q2_ = minmax(out_space_prev, last_incoming)
        lastcg3_om = getNsave_validout(S, (q1_, q2_))
        out_space_set = union(out_space_set, Set(lastcg3_om.out_spaces))
    end
    vo_obj = ValidOuts(S, Tuple(incom_spaces), sort(collect(out_space_set)))
    save_validout_sqlite(S, vo_obj)
    return vo_obj
end
