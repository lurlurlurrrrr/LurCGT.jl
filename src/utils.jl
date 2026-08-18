"""
    combination_index(m, n, target) -> Int

Return the lexicographic rank of one `n`-combination drawn from `1:m`.

`target` is the sorted selected combination. The result is one-based and counts
how many combinations would appear before `target` by summing binomial blocks at
each selected position.
"""
function combination_index(m, n, target)
    index = 1
    for (i, v) in enumerate(target)
        if i == 1
            start = 1
        else
            start = target[i-1] + 1
        end

        for x in start:v-1
            index += binomial(m - x, n - i)
        end
    end
    return index
end

"""
    inner_prod(v1, v2, basisnorms)

Compute a diagonal-basis inner product.

`v1` and `v2` are coefficient vectors in the same basis. `basisnorms` contains
the norm factor for each basis vector, so the result is
`dot(basisnorms, v1 .* v2)`.
"""
inner_prod(v1::AbstractVector, v2::AbstractVector, basisnorms) = 
dot(basisnorms, v1 .* v2)

"""
    orthog_2vecs(v1, v2, v1normsq, inprod)

Orthogonalize `v2` against `v1` using exact rational arithmetic.

`v1normsq` is `<v1,v1>` and `inprod` is `<v1,v2>`. The returned tuple
`(p, q, r, v)` represents an integer residual proportional to
`p * v2 - q * v1`, with common factor `r` removed from `v`.
"""
function orthog_2vecs(v1::AbstractVector, v2::AbstractVector, v1normsq, inprod) 
    rat = inprod // v1normsq; p, q = rat.den, rat.num
    v = p * v2 - q * v1
    r, v = divcfac(v)
    return p, q, r, v
end

"""
    contract_ith(arr::Array{CT,N}, mat::AbstractMatrix{RT}, ::Val{I})

Contract matrix `mat` into axis `I` of `arr`.

`arr` is an N-dimensional CGT-like coefficient array. `mat` is interpreted with
its first index contracted against axis `I`, replacing that axis by the matrix's
second index. `Val{I}` keeps the generated tensor expression type-stable.
"""
# mat: Matrix, inner product matrix or its inverse.
# arr: N-dimensional array. Ith index of arr and 1st index of mat is contracted.
@generated function contract_ith(arr::Array{CT, N}, 
    mat::AbstractMatrix{RT}, 
    ::Val{I}) where {RT<:Number, CT<:Number, N, I}
    @assert I <= N
    inds = [Symbol(:i, i) for i in 1:N]
    out_inds = copy(inds)
    out_inds[I] = :a  # Replace ith index with 'a'
    
    contraction = :(mat[$(inds[I]), a])
    
    quote
        @tensor cgt[$(out_inds...)] := arr[$(inds...)] * $contraction
        return cgt
    end
end
