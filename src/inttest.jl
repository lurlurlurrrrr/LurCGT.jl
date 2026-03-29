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
    #println(processed)

    for q in processed
        add_next_qlabels!(S, tobe_processed, tobe_processed_set, processed, q)
    end
    
    try
        while true
            imn = pop!(tobe_processed)
            q = imn.qlabel
            delete!(tobe_processed_set, q)
            @info "Processing irep $(totxt(S)), qlabel=$(q), dim=$(imn.dim)"
            #println(imn)
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
            #println(imn)
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

all_nonnegative_qlabel(q::NTuple{N, Int}) where {N} = all(x -> x >= 0, q)

function is_valid_int128_sum_slice_qlabel(::Type{S},
    q::NTuple{NZ, Int}) where {S<:Union{SU, Sp}, NZ}
    return all_nonnegative_qlabel(q)
end

function is_valid_int128_sum_slice_qlabel(::Type{SO{N}},
    q::NTuple{NZ, Int}) where {N, NZ}
    all_nonnegative_qlabel(q) || return false
    if isodd(N)
        return iseven(q[end])
    end
    return iseven(q[end - 1] + q[end])
end

function int128_sum_slice_qlabels(::Type{S},
    n::Int) where {S<:NonabelianSymm}

    n >= 0 || throw(ArgumentError("sum of qlabel entries must be nonnegative"))
    NZ = nzops(S)
    result = NTuple{NZ, Int}[]
    for comb in list_combinations(NZ, n)
        q = Tuple(comb)
        is_valid_int128_sum_slice_qlabel(S, q) || continue
        push!(result, q)
    end
    return sort!(unique!(result))
end

int128_sum_slice_normalize(n1::Int, n2::Int) = max(n1, n2), min(n1, n2)

function int128_sum_slice_pairs(::Type{S},
    n1::Int,
    n2::Int) where {S<:NonabelianSymm}

    n1_, n2_ = int128_sum_slice_normalize(n1, n2)
    qlabels1 = int128_sum_slice_qlabels(S, n1_)
    qlabels2 = int128_sum_slice_qlabels(S, n2_)
    pairs = Set{Tuple{NTuple{nzops(S), Int}, NTuple{nzops(S), Int}}}()
    for q1 in qlabels1, q2 in qlabels2
        push!(pairs, minmax(q1, q2))
    end
    return sort!(collect(pairs))
end

function int128_irep_equal_except_size_byte(rep1::Irep, rep2::Irep)
    return rep1.Sl == rep2.Sl &&
           rep1.Sz == rep2.Sz &&
           rep1.innerprod == rep2.innerprod &&
           rep1.inv_innerprod == rep2.inv_innerprod &&
           rep1.qlabel == rep2.qlabel &&
           rep1.dimension == rep2.dimension
end

function int128_sum_slice_irep_status(::Type{S},
    q::NTuple{NZ, Int};
    verbose=0) where {S<:NonabelianSymm, NZ}

    try
        rep_big = getNsave_irep(S, BigInt, q)
        rep_i128 = getNsave_irep(S, Int128, q)
        if int128_irep_equal_except_size_byte(rep_big, rep_i128)
            return (status=:match, reason="match")
        end
        verbose > 0 && println("IREP MISMATCH $(q)")
        return (status=:mismatch, reason="Int128 and BigInt irreps differ")
    catch err
        if err isa AssertionError || err isa OverflowError || err isa InexactError
            verbose > 0 && println("IREP FAIL $(q): $(typeof(err))")
            return (status=:error, reason=string(typeof(err)))
        end
        rethrow()
    end
end

function int128_sum_slice_print_summary(summary)
    println("Int128 sum-slice summary for $(summary.symmetry) with (n1, n2)=($(summary.n1), $(summary.n2)):")
    println("  total pairs: $(summary.total_pairs)")
    println("  passed pairs: $(summary.passed_pairs)")
    println("  failed pairs: $(summary.failed_pairs)")
    println("  irep mismatch pairs: $(summary.irep_mismatch_pairs)")
    println("  CGT failed pairs: $(summary.cgt_failed_pairs)")
    return nothing
end

function run_int128_cgt_sum_slice(::Type{S},
    n1::Int,
    n2::Int;
    verbose=1) where {S<:NonabelianSymm}

    n1_, n2_ = int128_sum_slice_normalize(n1, n2)
    pairs = int128_sum_slice_pairs(S, n1_, n2_)
    qlabels = sort!(unique!([q for pair in pairs for q in pair]))
    irep_status = Dict{NTuple{nzops(S), Int}, NamedTuple{(:status, :reason), Tuple{Symbol, String}}}()
    for q in qlabels
        irep_status[q] = int128_sum_slice_irep_status(S, q; verbose=max(verbose - 1, 0))
    end

    results = NamedTuple[]
    for pair in pairs
        q1, q2 = pair
        q1_status = irep_status[q1]
        q2_status = irep_status[q2]
        if q1_status.status != :match
            result = (q1=q1, q2=q2, status=:irep_mismatch, reason="q1: $(q1_status.reason)")
        elseif q2_status.status != :match
            result = (q1=q1, q2=q2, status=:irep_mismatch, reason="q2: $(q2_status.reason)")
        else
            try
                generate_every_CGT(S, Int128, Int128, pair, nothing; assertlev=1, save=false)
                result = (q1=q1, q2=q2, status=:passed, reason="passed")
            catch err
                if err isa AssertionError || err isa OverflowError || err isa InexactError
                    result = (q1=q1, q2=q2, status=:cgt_failed, reason=string(typeof(err)))
                else
                    rethrow()
                end
            end
        end
        push!(results, result)
        println("($(result.q1), $(result.q2)) => $(result.status): $(result.reason)")
        GC.gc()
    end

    failures = [result for result in results if result.status != :passed]
    summary = (
        symmetry=totxt(S),
        n1=n1_,
        n2=n2_,
        total_pairs=length(results),
        passed_pairs=count(result -> result.status == :passed, results),
        failed_pairs=length(failures),
        irep_mismatch_pairs=count(result -> result.status == :irep_mismatch, results),
        cgt_failed_pairs=count(result -> result.status == :cgt_failed, results),
        results=results,
        failures=failures,
    )
    verbose > 0 && int128_sum_slice_print_summary(summary)
    return summary
end

Base.show(io::IO, imn::irep_maxnums) =
    print(io, "irep_maxnums(qlabel:$(imn.qlabel), dim:$(imn.dim), 
    Sl:$(imn.matelem_max), inprod:$(imn.innerprod_max), 
    inv_inprod:$(imn.inv_innerprod_max), inv_inprod_fac:$(imn.inv_innerprod_fac_max))")
