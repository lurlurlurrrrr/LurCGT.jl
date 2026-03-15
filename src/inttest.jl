using Serialization
using DataStructures

totxt(::Type{BigInt}) = "BigInt"
totxt(::Type{Int128}) = "Int128"
totxt(::Type{Int64}) = "Int64"

struct irep_maxnums{S<:NonabelianSymm, NZ}
    qlabel::NTuple{NZ, Int}
    # Dimension of the irrep. Sorted by dimension
    dim::Integer
    # The largest number appearing when getting matrix elements
    matelem_max::Integer
    innerprod_max::Integer
    inv_innerprod_max::Integer
    inv_innerprod_fac_max::Integer
end

function irep_maxnums(::Type{S},
    qlabel::NTuple{NZ, Int},
    dim::Integer,
    matelem_max::Integer,
    innerprod_max::Integer,
    inv_innerprod_max::Integer,
    inv_innerprod_fac_max::Integer) where {S<:NonabelianSymm, NZ}

    return irep_maxnums{S, NZ}(qlabel, dim, matelem_max, 
    innerprod_max, inv_innerprod_max, inv_innerprod_fac_max)
end

function list_combinations(n::Int, N::Int)
    if n == 1
        return [[N]]
    end

    result = Vector{Vector{Int}}()
    for k in 0:N
        for tail in list_combinations(n - 1, N - k)
            push!(result, [k; tail])
        end
    end
    return result
end

function analyze_irep(irep::Irep{S, NL, NZ, RT}) where {S<:NonabelianSymm, NL, NZ, RT}
    matelem_max = 0
    for sl in irep.Sl
        for (_, mat) in sl
            matelem_max = max(matelem_max, maximum(abs.(mat.nzval)))
        end
    end

    innerprod_max = 0
    for (_, mat) in irep.innerprod
        innerprod_max = max(innerprod_max, maximum(abs.(mat)))
    end

    inv_innerprod_max = 0
    inv_innerprod_fac_max = 0
    for (_, (mat, fac)) in irep.inv_innerprod
        inv_innerprod_max = max(inv_innerprod_max, maximum(abs.(mat)))
        inv_innerprod_fac_max = max(inv_innerprod_fac_max, abs(numerator(fac)))
        inv_innerprod_fac_max = max(inv_innerprod_fac_max, abs(denominator(fac)))
    end

    return irep_maxnums{S, NZ}(irep.qlabel,
        irep.dimension,
        matelem_max,
        innerprod_max,
        inv_innerprod_max,
        inv_innerprod_fac_max)
end

Base.isless(a::irep_maxnums{S, NZ}, 
    b::irep_maxnums{S, NZ}) where {S<:NonabelianSymm, NZ} = a.dim < b.dim

function get_maximum(a::Integer, b::AbstractArray{<:Integer}) 
    bmax = maximum(abs.(b))
    return a > bmax ? a : bmax
end

function add_next_qlabels!(::Type{S}, 
    tobe_processed::BinaryMinHeap{irep_maxnums{S, NZ}},
    tobe_processed_set::Set{NTuple{NZ, Int}},
    processed::Set{NTuple{NZ, Int}}, 
    q::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ}
    
    funda_lst = fundamental_qlabels(S)
    for f in funda_lst
        nq = q .+ f
        if !in(nq, processed) && !in(nq, tobe_processed_set)
            irep = getNsave_irep(S, BigInt, nq)
            imn = analyze_irep(irep)
            push!(tobe_processed, imn)
            push!(tobe_processed_set, nq)
        end
    end
end

function test_irep_numtype(::Type{S}) where {S<:NonabelianSymm}
    fname = "$(homedir())/LurCGT_inttest/irep_maxnums_$(totxt(S))"
    NZ = nzops(S)
    processed = Set{NTuple{NZ, Int}}()
    tobe_processed = BinaryMinHeap{irep_maxnums{S, NZ}}()
    tobe_processed_set = Set{NTuple{NZ, Int}}()
    v = Vector{irep_maxnums{S, NZ}}()
    if isfile(fname) 
        v::Vector{irep_maxnums{S, NZ}} = deserialize(fname) 
        for imn in v push!(processed, imn.qlabel) end
    else 
        push!(processed, Tuple(0 for _=1:NZ)) 
    end
    println(processed)

    for q in processed
        add_next_qlabels!(S, tobe_processed, tobe_processed_set, processed, q)
    end
    
    try
        while true
            imn = pop!(tobe_processed)
            q = imn.qlabel
            delete!(tobe_processed_set, q)
            @info "Processing irep $(totxt(S)), qlabel=$(q), dim=$(imn.dim)"
            println(imn)
            push!(processed, q)
            push!(v, imn)
            add_next_qlabels!(S, tobe_processed, tobe_processed_set, processed, q)
        end
    catch
        serialize(fname, v)
    end
end

function test_irep_numtype(::Type{S},
    ::Type{RT},
    smin::Int,
    smax::Int) where {S<:NonabelianSymm, RT<:Number}

    NZ = nzops(S)
    for i in smin:smax
        for comb in list_combinations(NZ, i)
            qlabel = Tuple(comb)
            @info "Testing irep number type $(totxt(S)), RT=$(totxt(RT)), qlabel=$(qlabel)"
            rep = getNsave_irep(S, RT, qlabel)
            imn = analyze_irep(rep)
            println(imn)
            if imn.matelem_max > typemax(RT) || 
               imn.innerprod_max > typemax(RT) || 
               imn.inv_innerprod_max > typemax(RT) || 
               imn.inv_innerprod_fac_max > typemax(RT)
                @error "Number type $(totxt(RT)) is insufficient for irep $(totxt(S)), qlabel=$(qlabel)"
            else
                @info "Number type $(totxt(RT)) is sufficient for irep $(totxt(S)), qlabel=$(qlabel)"
            end
            
        end
    end
end

function test_cg3_numtype(::Type{S},
    ::Type{CT},
    smin::Int,
    smax::Int) where {S<:NonabelianSymm, CT<:Number}

    NZ = nzops(S)
    for j in smin:smax
        for i in 0:j
            for comb1 in list_combinations(NZ, i)
                for comb2 in list_combinations(NZ, j)
                    q1, q2 = Tuple(comb1), Tuple(comb2)
                    @info "Testing CG3 number type $(totxt(S)), CT=$(totxt(CT)), qlabels=($(q1), $(q2))"
                    generate_every_CGT(S, CT, CT, minmax(q1, q2), nothing; assertlev=1, save=false)
                end
            end
        end
    end
end

Base.show(io::IO, imn::irep_maxnums) =
    print(io, "irep_maxnums(qlabel:$(imn.qlabel), dim:$(imn.dim), 
    Sl:$(imn.matelem_max), inprod:$(imn.innerprod_max), 
    inv_inprod:$(imn.inv_innerprod_max), inv_inprod_fac:$(imn.inv_innerprod_fac_max))")