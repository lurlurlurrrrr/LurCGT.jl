# Test permutation using F/R-symbol works well
function test_permutation(::Type{S}, ins, out, perm, ::Type{FT}; verbose=0) where {S<:NonabelianSymm, FT<:AbstractFloat}
    FTree = LurCGT.random_FTree(S, Tuple(ins), (out,))

    # First, get an Float array and permute
    arr_perm = FTree2arr(FTree, FT)
    arr_perm = permutedims(arr_perm, perm)

    # Second, get a permuted FTree and convert it to array
    switch_list = decompose_perm(perm)
    for (i, j) in switch_list
        FTree = LurCGT.permute_adjacent!(FTree, i, j; verbose)
    end
    arr_perm2 = FTree2arr(FTree, FT)

    # The two results should be the same
    return arr_perm, arr_perm2, norm(arr_perm - arr_perm2)
end

function decompose_perm(perm)
    len = length(perm)
    vec = collect(1:len)
    switch_list = Vector{Tuple{Int, Int}}()

    for i=1:len
        target_elem = perm[i]
        start_idx = findfirst(x -> x == target_elem, vec)
        @assert start_idx !== nothing "Invalid permutation: $perm"
        if start_idx == i continue end  # Already in place
        # Swap the elements
        for k=start_idx-1:-1:i
            vec[k], vec[k+1] = vec[k+1], vec[k]
            push!(switch_list, (k, k+1))
        end
        if Tuple(vec) == perm return switch_list end
    end
end


function test_permutation_randinput(::Type{S}, dim_limit=100000, outlimit=100, test_inputs=20; verbose=0) where S<:NonabelianSymm
    for i in 1:test_inputs
        ins = sort(generate_incom(S, dim_limit))
        N = length(ins)
        if N < 2 continue end
        vo = getNsave_validout(S, Tuple(ins))
        out = select_out(vo, outlimit)
        perm = (get_rand_perm(N)..., N+1)

        _, _, norm_diff = test_permutation(S, ins, out, perm, BigFloat; verbose)
        if norm_diff < 1e-50
            println("Test #$i: Permutation test passed for: $(ins) -> $(out) with permutation $(perm) with error $(norm_diff)\n")
        else
            error("Test #$i: Permutation test failed for: $(ins) -> $(out) with permutation $(perm) with error $(norm_diff)")
        end
    end
end

function get_rand_perm(N)
    while true
        s = Set(1:N)
        perm = Int[]
        for _=1:N
            i = rand(s)
            push!(perm, i)
            delete!(s, i)
        end
        if perm != collect(1:N) return Tuple(perm) end # not identity permutation
    end
end

