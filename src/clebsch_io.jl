# NZ : length of qlabels
# N : Dimension of CGT 
# I defined it to cover arbitrary CGT dimensions
# but in my calculation, it is 3 for almost every case
# 2 for rare cases (1j-symbol)
# Clebsch-Gordan tensor in block-sparse format
# Each block is a dense array, and the blocks are indexed by weights
# RT : Number type of irreps, CT : Number type of CGT
struct CGT{S<:NonabelianSymm, CT<:Number, NZ, N}
    qlabels::NTuple{N, NTuple{NZ, Int}} # qlabels for each leg
    # Nonzero blocks of CGT. 
    # Key is a tuple of qlabels, value is the corresponding array
    # Each array has rank N+1, the last dimension is outer multiplicity
    blocks::Dict{NTuple{N, NTuple{NZ, Int}}, Array{CT}}

    nfactor::Vector{Rational{BigInt}} 
    dir::NTuple{N, Char}
    # Cached memory footprint in bytes, computed once at construction for LRU cache eviction by size
    size_byte::Int
    
    function CGT{S, CT, NZ, N}(qlabels, blocks, nfactor, dir, size_byte::Int=0) where {S, CT, NZ, N}
        if size_byte == 0
            obj = new{S, CT, NZ, N}(qlabels, blocks, nfactor, dir, 0)
            size_byte = Base.summarysize(obj)
        end
        new{S, CT, NZ, N}(qlabels, blocks, nfactor, dir, size_byte)
    end
end


function load_cg3blk(::Type{S},
    ::Type{CT},
    qlabels_in::NTuple{2, NTuple{NZ, Int}},
    qlabels_out::Vector{NTuple{NZ, Int}}) where {S<:NonabelianSymm, NZ, CT<:Number}

    cg3s = getNsave_cg3(S, CT, qlabels_in, qlabels_out)
    return Dict(q => (cg3s[q].blocks, cg3s[q].nfactor) for q in qlabels_out)
end

function getNsave_cg3(::Type{S}, 
    ::Type{CT}, 
    qlabels_in::NTuple{2, NTuple{NZ, Int}},
    qlabels_out::Vector{NTuple{NZ, Int}}) where {S<:NonabelianSymm, CT<:Number, NZ}

    @assert NZ == nzops(S)
    dir = ['+', '+', '-']
    qlabels_in_sorted = Tuple([minmax(qlabels_in[1], qlabels_in[2])...])
    permute = qlabels_in[1] > qlabels_in[2]

    cg3s = Dict{NTuple{NZ, Int}, CGT{S, CT, NZ, 3}}()
    remain_qs = Vector{NTuple{NZ, Int}}()
    # Find as many CG3s as possible from saved files
    for qout in qlabels_out
        qlabels_sorted = Tuple([qlabels_in_sorted..., qout])
        cg3 = load_cgt_sqlite(S, CT, qlabels_sorted)
        if !isnothing(cg3) cg3s[qout] = cg3
        else push!(remain_qs, qout) end
    end

    # Generate missing CG3s
    if !isempty(remain_qs)
        new_cg3s = generate_every_CGT(S, CT, CT, qlabels_in_sorted, remain_qs)
        for q in remain_qs 
            @assert haskey(new_cg3s, q)
            cg3s[q] = new_cg3s[q]
        end
    end
    
    for q in remain_qs @assert haskey(cg3s, q) end

    if permute
        for (q, cg3) in cg3s
            new_blocks = Dict{NTuple{3, NTuple{NZ, Int}}, Array{CT}}()
            for k in keys(cg3.blocks)
                new_key = (k[2], k[1], k[3])
                new_blocks[new_key] = permutedims(cg3.blocks[k], (2, 1, 3, 4))
            end
            cg3s[q] = CGT{S, CT, NZ, 3}(Tuple([qlabels_in..., q]), 
            new_blocks, cg3.nfactor, Tuple(dir), cg3.size_byte)
        end
    end
    return cg3s
end
