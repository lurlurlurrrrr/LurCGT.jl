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

inner_prod(v1::AbstractVector, v2::AbstractVector, basisnorms) = 
dot(basisnorms, v1 .* v2)

function orthog_2vecs(v1::AbstractVector, v2::AbstractVector, v1normsq, inprod) 
    rat = inprod // v1normsq; p, q = rat.den, rat.num
    v = p * v2 - q * v1
    r, v = divcfac(v)
    return p, q, r, v
end

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