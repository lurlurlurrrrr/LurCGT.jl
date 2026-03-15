lt_multisymm_weights(t1::NTuple{N, Tuple{Vararg{Int}}},
    t2::NTuple{N, Tuple{Vararg{Int}}}) where N =
    rev_less_symms(reverse(t1), reverse(t2))


function rev_less_symms(t1::NTuple{N, Tuple{Vararg{Int}}},
    t2::NTuple{N, Tuple{Vararg{Int}}}) where N 

    for i in 1:N
        if reverse(t1[i]) != reverse(t2[i])
            return reverse(t1[i]) < reverse(t2[i])
        end
    end
    return false
end

rev_less(t1::NTuple{N, Int}, t2::NTuple{N, Int}) where N = reverse(t1) < reverse(t2)

sorted_zvals(Sz::Dict{NTuple{NZ, Int}}) where NZ =
    sort!(collect(keys(Sz)); lt=rev_less, rev=true)

function sortperm_sz(Sz::NTuple{N, Vector{Int}}) where N
    tuples = collect(zip(Sz...))
    return sortperm(tuples; lt=rev_less, rev=true)
end

# Input : permutation vector
# Return 1 if it is an even permutation, -1 otherwise (odd permutation)
function permutation_sign(perm)
    n = length(perm)
    visited = falses(n)
    swaps = 0

    for i in 1:n
        if visited[i]
            continue
        end
        cycle_size = 0
        j = i
        while !visited[j]
            visited[j] = true
            j = perm[j]
            cycle_size += 1
        end
        swaps += cycle_size - 1
    end

    return (-1)^swaps  
end

