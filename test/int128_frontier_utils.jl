int128_frontier_cover_steps(::Type{SU{N}}) where {N} =
    [ntuple(i -> i == j ? 1 : 0, N - 1) for j in 1:N-1]

int128_frontier_cover_steps(::Type{Sp{N}}) where {N} =
    [ntuple(i -> i == j ? 1 : 0, div(N, 2)) for j in 1:div(N, 2)]

function int128_frontier_cover_steps(::Type{SO{N}}) where {N}
    NZ = div(N, 2)
    steps = NTuple{NZ, Int}[]
    upper = isodd(N) ? NZ - 1 : NZ - 2
    for j in 1:max(upper, 0)
        push!(steps, ntuple(i -> i == j ? 1 : 0, NZ))
    end
    if isodd(N)
        push!(steps, ntuple(i -> i == NZ ? 2 : 0, NZ))
    else
        push!(steps, ntuple(i -> i == NZ - 1 ? 2 : 0, NZ))
        push!(steps, ntuple(i -> i == NZ ? 2 : 0, NZ))
        push!(steps, ntuple(i -> i >= NZ - 1 ? 1 : 0, NZ))
    end
    return sort!(unique!(steps))
end

is_nonzero_qlabel(q) = any(!iszero, q)
all_nonnegative_qlabel(q) = all(>=(0), q)

function is_valid_int128_frontier_qlabel(::Type{S}, q::NTuple{NZ, Int}) where {S<:Union{SU, Sp}, NZ}
    return all_nonnegative_qlabel(q) && is_nonzero_qlabel(q)
end

function is_valid_int128_frontier_qlabel(::Type{SO{N}}, q::NTuple{NZ, Int}) where {N, NZ}
    all_nonnegative_qlabel(q) || return false
    is_nonzero_qlabel(q) || return false
    if isodd(N)
        return iseven(q[end])
    end
    return iseven(q[end - 1] + q[end])
end

canonicalize_int128_frontier_pair(q1, q2) = minmax(q1, q2)

add_qlabels(q::NTuple{NZ, Int}, step::NTuple{NZ, Int}) where {NZ} = q .+ step
sub_qlabels(q::NTuple{NZ, Int}, step::NTuple{NZ, Int}) where {NZ} = q .- step

function int128_frontier_successors(::Type{S},
    pair::Tuple{NTuple{NZ, Int}, NTuple{NZ, Int}}) where {S<:NonabelianSymm, NZ}

    q1, q2 = pair
    successors = Set{typeof(pair)}()
    for step in int128_frontier_cover_steps(S)
        next1 = add_qlabels(q1, step)
        if is_valid_int128_frontier_qlabel(S, next1)
            push!(successors, canonicalize_int128_frontier_pair(next1, q2))
        end

        next2 = add_qlabels(q2, step)
        if is_valid_int128_frontier_qlabel(S, next2)
            push!(successors, canonicalize_int128_frontier_pair(q1, next2))
        end
    end
    return sort!(collect(successors); by=int128_frontier_sort_key)
end

function int128_frontier_predecessors(::Type{S},
    pair::Tuple{NTuple{NZ, Int}, NTuple{NZ, Int}}) where {S<:NonabelianSymm, NZ}

    q1, q2 = pair
    predecessors = Set{typeof(pair)}()
    for step in int128_frontier_cover_steps(S)
        prev1 = sub_qlabels(q1, step)
        if all_nonnegative_qlabel(prev1) && is_valid_int128_frontier_qlabel(S, prev1)
            push!(predecessors, canonicalize_int128_frontier_pair(prev1, q2))
        end

        prev2 = sub_qlabels(q2, step)
        if all_nonnegative_qlabel(prev2) && is_valid_int128_frontier_qlabel(S, prev2)
            push!(predecessors, canonicalize_int128_frontier_pair(q1, prev2))
        end
    end
    return sort!(collect(predecessors); by=int128_frontier_sort_key)
end

function int128_frontier_seed_pairs(::Type{S}) where {S<:NonabelianSymm}
    steps = int128_frontier_cover_steps(S)
    seeds = Set{Tuple{eltype(steps), eltype(steps)}}()
    for q1 in steps, q2 in steps
        push!(seeds, canonicalize_int128_frontier_pair(q1, q2))
    end
    return sort!(collect(seeds); by=int128_frontier_sort_key)
end

function int128_frontier_sort_key(pair)
    q1, q2 = pair
    return (sum(q1) + sum(q2), q1, q2)
end

function int128_frontier_irep_equal_except_size_byte(rep1::Irep, rep2::Irep)
    return rep1.Sl == rep2.Sl &&
           rep1.Sz == rep2.Sz &&
           rep1.innerprod == rep2.innerprod &&
           rep1.inv_innerprod == rep2.inv_innerprod &&
           rep1.qlabel == rep2.qlabel &&
           rep1.dimension == rep2.dimension
end

function int128_frontier_irep_matches(::Type{S},
    q::NTuple{NZ, Int};
    verbose=0) where {S<:NonabelianSymm, NZ}

    try
        rep_big = getNsave_irep(S, BigInt, q)
        rep_i128 = getNsave_irep(S, Int128, q)
        same = int128_frontier_irep_equal_except_size_byte(rep_big, rep_i128)
        if verbose > 0 && !same
            println("IREP MISMATCH $(q)")
        end
        return same
    catch err
        if err isa AssertionError || err isa OverflowError || err isa InexactError
            if verbose > 0
                println("IREP FAIL $(q): $(typeof(err))")
            end
            return false
        end
        rethrow()
    end
end

function int128_frontier_passes(::Type{S},
    q1::NTuple{NZ, Int},
    q2::NTuple{NZ, Int};
    verbose=0) where {S<:NonabelianSymm, NZ}

    pair = canonicalize_int128_frontier_pair(q1, q2)
    int128_frontier_irep_matches(S, pair[1]; verbose) || return false
    pair[1] == pair[2] || int128_frontier_irep_matches(S, pair[2]; verbose) || return false
    try
        LurCGT.generate_every_CGT(S, Int128, Int128, pair, nothing; assertlev=1, save=false)
        verbose > 0 && println("PASS $(pair)")
        return true
    catch err
        if err isa AssertionError || err isa OverflowError || err isa InexactError
            verbose > 0 && println("FAIL $(pair): $(typeof(err))")
            return false
        end
        rethrow()
    end
end

function int128_frontier_candidate_status(::Type{S},
    pair,
    passed,
    failed,
    blocked) where {S<:NonabelianSymm}

    predecessors = int128_frontier_predecessors(S, pair)
    any(pred -> pred in failed || pred in blocked, predecessors) && return :blocked
    all(pred -> pred in passed, predecessors) && return :ready
    return :wait
end

function find_int128_cgt_frontier(::Type{S};
    pass_predicate=nothing,
    max_tested::Int=typemax(Int),
    verbose=0) where {S<:NonabelianSymm}

    tested = 0
    pending = Set(int128_frontier_seed_pairs(S))
    passed = Set{eltype(pending)}()
    failed = Set{eltype(pending)}()
    blocked = Set{eltype(pending)}()
    pass_fn = isnothing(pass_predicate) ? ((q1, q2) -> int128_frontier_passes(S, q1, q2; verbose)) : pass_predicate

    while !isempty(pending)
        blocked_now = eltype(pending)[]
        ready_now = eltype(pending)[]
        for pair in pending
            status = int128_frontier_candidate_status(S, pair, passed, failed, blocked)
            if status == :blocked
                push!(blocked_now, pair)
            elseif status == :ready
                push!(ready_now, pair)
            end
        end

        for pair in blocked_now
            delete!(pending, pair)
            push!(blocked, pair)
        end

        isempty(pending) && break
        isempty(ready_now) || sort!(ready_now; by=int128_frontier_sort_key)
        if isempty(ready_now)
            error("No eligible Int128 frontier candidates remain; search is stuck with $(length(pending)) pending states.")
        end

        pair = first(ready_now)
        delete!(pending, pair)
        tested += 1
        tested <= max_tested || error("Int128 frontier search exceeded max_tested=$(max_tested)")

        q1, q2 = pair
        if pass_fn(q1, q2)
            push!(passed, pair)
            verbose > 0 && println("tested $(tested): pass $(pair)")
            for successor in int128_frontier_successors(S, pair)
                if !(successor in passed || successor in failed || successor in blocked)
                    push!(pending, successor)
                end
            end
        else
            push!(failed, pair)
            verbose > 0 && println("tested $(tested): fail $(pair)")
        end
    end

    frontier = [pair for pair in passed if !any(successor -> successor in passed, int128_frontier_successors(S, pair))]
    sort!(frontier; by=int128_frontier_sort_key)
    return frontier
end
