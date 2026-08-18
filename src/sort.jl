"""
    lt_multisymm_weights(t1, t2) -> Bool

Compare product-symmetry weight tuples in the package's deterministic order.

`t1` and `t2` contain one weight tuple per symmetry. The symmetry order is
reversed first, then each individual weight is compared with reversed coordinate
order. This matches the ordering expected by CGT/TLArray sector metadata.
"""
lt_multisymm_weights(t1::NTuple{N, Tuple{Vararg{Int}}},
    t2::NTuple{N, Tuple{Vararg{Int}}}) where N =
    rev_less_symms(reverse(t1), reverse(t2))


"""
    rev_less_symms(t1, t2) -> Bool

Compare same-length tuples of symmetry weights by reversed coordinates.

The first component whose reversed coordinate tuple differs determines the
ordering. Equal tuples return `false`, following Julia `lt` predicate
conventions.
"""
function rev_less_symms(t1::NTuple{N, Tuple{Vararg{Int}}},
    t2::NTuple{N, Tuple{Vararg{Int}}}) where N 

    for i in 1:N
        if reverse(t1[i]) != reverse(t2[i])
            return reverse(t1[i]) < reverse(t2[i])
        end
    end
    return false
end

"""
    rev_less(t1, t2) -> Bool

Compare two weight tuples after reversing coordinate order.

This is the scalar-weight ordering used for sorted z-sector traversal in irrep
and Clebsch generation.
"""
rev_less(t1::NTuple{N, Int}, t2::NTuple{N, Int}) where N = reverse(t1) < reverse(t2)

"""
    sorted_zvals(Sz::Dict) -> Vector

Return z-weight keys from `Sz` in descending package weight order.

`Sz` is an irrep sector dictionary mapping z-weight tuples to basis ranges.
Only keys are used; the returned vector is sorted by `rev_less` with
`rev=true`.
"""
sorted_zvals(Sz::Dict{NTuple{NZ, Int}}) where NZ =
    sort!(collect(keys(Sz)); lt=rev_less, rev=true)

"""
    sortperm_sz(Sz::NTuple{N, Vector{Int}}) -> Vector{Int}

Return the permutation that sorts zipped z-weight coordinate vectors.

`Sz` contains one coordinate vector per Cartan direction. The function zips
them into weight tuples and returns the sort permutation in descending
`rev_less` order.
"""
function sortperm_sz(Sz::NTuple{N, Vector{Int}}) where N
    tuples = collect(zip(Sz...))
    return sortperm(tuples; lt=rev_less, rev=true)
end

"""
    permutation_sign(perm) -> Int

Return `1` for an even permutation and `-1` for an odd permutation.

`perm` is a one-based permutation vector. The implementation decomposes it into
cycles and counts the number of swaps implied by those cycles.
"""
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
