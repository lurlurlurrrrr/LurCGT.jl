⊗(a::Matrix, b::Matrix) = kron(b, a)
⊗(a::Vector, b::Vector) = kron(b, a)

comm(A, B) = A * B - B * A 

# Test about F-symbol.

# Get F_symbol and convert it into real matrix form
function Fsymbol_toreal(::Type{S}, 
    ::Type{CT}, 
    ::Type{FT}, 
    in1::NTuple{NZ, Int}, 
    in2::NTuple{NZ, Int}, 
    in3::NTuple{NZ, Int}, 
    out::NTuple{NZ, Int}) where {S<:NonabelianSymm, CT<:Number, FT<:Number, NZ}
    fsym = getNsave_Fsymbol(S, CT, in1, in2, in3, out)

    mat = fsym.fsym_mat
    sz = size(mat)[1]
    realmat = zeros(FT, sz, sz)

    for i=1:sz
        for j=1:sz
            realmat[i, j] = mat[i, j] * sqrt(fsym.fs_nfactor[j]*fsym.es_nfactor[i])
        end
    end
    return realmat
end

function Fsymbol_fullarr(::Type{S}, 
    in1::NTuple{NZ, Int}, 
    in2::NTuple{NZ, Int}, 
    in3::NTuple{NZ, Int}, 
    out::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ}

    # in1 * in2 -> e, e * in3 -> out
    es_list = LurCGT.find_espaces(S, in1, in2, in3, out, BigInt)
    # in2 * in3 -> e, e * in1 -> out
    fs_list = LurCGT.find_espaces(S, in2, in3, in1, out, BigInt)

    # Get the outer multiplicity
    mat_size = 0
    for (_, (in12_om, qin3_om)) in es_list
        mat_size += in12_om * qin3_om
    end
    
    # Get the outer multiplicity in other way
    mat_size_f = 0
    for (_, (qin2_om, in3_om)) in fs_list
        mat_size_f += qin2_om * in3_om
    end

    @assert mat_size == mat_size_f
    # Initialize the F-symbol matrix
    fsym_mat = zeros(Float64, mat_size, mat_size)

    es_contract_res = Array{Float64, 6}[]
    out_dim = 0
    # Fill in the es_contract_res
    for (e, _) in es_list
        # Load CG3 in1 ⊗ in2 -> e 
        μblk, _, _ = load_cg3_float(S, BigInt, (in1, in2, e))
        # Load CG3 e ⊗ in3 -> out
        νblk, _, _ = load_cg3_float(S, BigInt, (e, in3, out))

        @tensor blk_contract[in1, in2, in3, out, ν, μ] := 
            μblk[in1, in2, e, μ] * νblk[e, in3, out, ν]
        push!(es_contract_res, blk_contract)
        out_dim = size(blk_contract, 4)
    end

    fs_contract_res = Array{Float64, 6}[]
    for (f, _) in fs_list
        # Load CG3 in2 ⊗ in3 -> f
        κblk, _, _ = load_cg3_float(S, BigInt, (in2, in3, f))
        # Load CG3 f ⊗ in1 -> out
        λblk, _, _ = load_cg3_float(S, BigInt, (f, in1, out))

        @tensor blk_contract[in1, in2, in3, out, λ, κ] := 
            κblk[in2, in3, f, κ] * λblk[f, in1, out, λ]
        push!(fs_contract_res, blk_contract)
        @assert out_dim == size(blk_contract, 4)
    end

    i, j = 1, 1
    for (k1, eblk) in enumerate(es_contract_res)
        _, (m1, m2) = es_list[k1]
        for (k2, fblk) in enumerate(fs_contract_res)
            _, (n1, n2) = fs_list[k2]
            irange = i:(i + m1*m2 - 1)
            jrange = j:(j + n1*n2 - 1)

            @tensor fsym[ν, μ, λ, κ] := 
                eblk[i1, i2, i3, out, ν, μ] * 
                fblk[i1, i2, i3, out, λ, κ]
            fsym_mat[irange, jrange] += reshape(fsym, m1*m2, n1*n2)
            j += n1 * n2
        end
        i += m1 * m2
        j = 1
    end

    return fsym_mat / out_dim
end

function compare_two_ways_fsym()
    in1, in2, in3, out = (2, 2), (2, 2), (4, 4), (4, 4)
    @time fsym_float1 = Fsymbol_toreal(SU{3}, BigInt, in1, in2, in3, out)
    @time fsym_float2 = Fsymbol_fullarr(SU{3}, in1, in2, in3, out)
    display(fsym_float1)
    display(fsym_float2)
    println(norm(fsym_float1 - fsym_float2))
    return fsym_float1, fsym_float2
end

function Rsymbol_toreal(::Type{S}, 
    ::Type{CT}, 
    in::NTuple{NZ, Int}, 
    out::NTuple{NZ, Int}) where {S<:NonabelianSymm, CT<:Number, NZ}
    rsym = getNsave_Rsymbol(S, CT, in, out)
    mat = rsym.rsym_mat
    
    sz = size(mat)[1]
    realmat = zeros(sz, sz)
    
    for i=1:sz
        for j=1:sz
            realmat[i, j] = mat[i, j] * sqrt(rsym.nfactor[i] * rsym.nfactor[j])
        end
    end
    return realmat
end

function Rsymbol_fullarr(::Type{S},
    in::NTuple{NZ, Int},
    out::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ}

    # in ⊗ in -> out
    blk, _, _ = load_cg3_float(S, BigInt, (in, in, out))
    out_dim = size(blk, 3)
    @tensor rsym[ν, μ] := blk[i2, i1, out, μ] * blk[i1, i2, out, ν]
    return Matrix{Float64}(rsym / out_dim)
end

function compare_two_ways_rsym()
    #in, out = (4, 4), (4, 4)
    in, out = (1,), (0,)
    @time rsym_float1 = Rsymbol_toreal(SU{2}, BigInt, in, out)
    @time rsym_float2 = Rsymbol_fullarr(SU{2}, in, out)
    display(rsym_float1)
    display(rsym_float2)
    println(norm(rsym_float1 - rsym_float2))
    return rsym_float1, rsym_float2
end

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

# Test permutation with random inputs
function generate_incom(::Type{S}, dim_limit, qlimit=4) where S<:NonabelianSymm
    NZ = nzops(S)
    ins = Vector{NTuple{NZ, Int}}()
    dim = 1
    while true
        qlabel = Tuple(rand(0:qlimit) for _=1:NZ)
        irep_dim = LurCGT.dimension(getNsave_irep(S, BigInt, qlabel))
        prod_dim = dim * irep_dim
        if prod_dim > dim_limit break end
        push!(ins, qlabel)
        dim = prod_dim
    end
    return ins
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

function dimension(::Type{S}, qlabel::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ}
    irep = getNsave_irep(S, BigInt, qlabel)
    return LurCGT.dimension(irep)
end

function select_out(vo::ValidOuts{S, N, NZ},
    outlimit=100) where {S<:NonabelianSymm, N, NZ}
    outlist = vo.out_spaces
    outdims = [dimension(S, out) for out in outlist]
    if minimum(outdims) > outlimit
        i = argmin(outdims)
        return outlist[i]
    end
    while true
        out = rand(vo.out_spaces)
        if dimension(S, out) <= outlimit return out end
    end
end

function select_out(::Type{S},
    vec::Vector{NTuple{NZ, Int}},
    outlimit=100) where {S<:NonabelianSymm, NZ}

    outdims = [dimension(S, out) for out in vec]
    if minimum(outdims) > outlimit
        i = argmin(outdims)
        return vec[i]
    end
    while true
        out = rand(vec)
        if dimension(S, out) <= outlimit return out end
    end
end

# Getting random F-symbol
function generate_incom_3spaces(::Type{S}, qlimit=4) where S<:NonabelianSymm
    NZ = nzops(S)
    ins = Vector{NTuple{NZ, Int}}()
    for _=1:3
        qlabel = Tuple(rand(0:qlimit) for _=1:NZ)
        push!(ins, qlabel)
    end

    ins_sorted = sort(ins)
    vo = getNsave_validout(S, Tuple(ins_sorted))
    out = select_out(vo, 100)

    return ins, out
end

function get_random_Fsymbol(::Type{S}, ::Type{FT}, qlimit=4) where {S<:NonabelianSymm, FT<:Integer}
    ins, out = generate_incom_3spaces(S, qlimit)
    return ins, out, getNsave_Fsymbol(S, BigInt, ins[1], ins[2], ins[3], out)
end

function get_random_Fsymbol_real(::Type{S}, ::Type{FT}, qlimit=4) where {S<:NonabelianSymm, FT<:AbstractFloat}
    ins, out = generate_incom_3spaces(S, qlimit)
    return ins, out, Fsymbol_toreal(S, BigInt, BigFloat, ins[1], ins[2], ins[3], out)
end

function get_random_Rsymbol_real(::Type{S}, ::Type{FT}, qlimit=4) where {S<:NonabelianSymm, FT<:AbstractFloat}
    NZ = nzops(S)
    qlabel = Tuple(rand(0:qlimit) for _=1:NZ)
    vo = getNsave_validout(S, (qlabel, qlabel))
    out = select_out(vo, 100)
    return qlabel, out, Rsymbol_toreal(S, BigInt, qlabel, out)
end

@generated function contract_om_step3(arr_up, 
    arr_dn, 
    om_arr, 
    ::Val{NI}, 
    ::Val{NO}) where {NI, NO}

    result_inds = [[Symbol(:in, i) for i=1:NI]..., [Symbol(:out, i) for i=1:NO]...]
    up_inds = [[Symbol(:in, i) for i=1:NI]..., :c, [Symbol(:oup, i) for i=1:NI-1]...]
    dn_inds = [[Symbol(:out, i) for i=1:NO]..., :c, [Symbol(:odn, i) for i=1:NO-1]...]
    om_inds = [[Symbol(:oup, i) for i=1:NI-1]..., [Symbol(:odn, i) for i=1:NO-1]...]

    quote
        @tensor result[$(result_inds...)] :=
            arr_up[$(up_inds...)] * arr_dn[$(dn_inds...)] * om_arr[$(om_inds...)]
        return result
    end
end

function get_spaces(::Type{S}, 
    dim_limit::Int, 
    qlimit::Int, 
    outlimit::Int; 
    verbose) where S<:NonabelianSymm

    NZ = nzops(S)
    # Get input and output of the first CGT
    incom = sort(generate_incom(S, dim_limit, qlimit))
    len = length(incom)
    if len < 3 return nothing end
    # The number of outgoing legs of CGT1
    cgt1_up = rand(1:len-1)
    # The number of incoming legs of CGT1
    cgt1_dn = len - cgt1_up

    shuffle_len = shuffle(1:len)
    
    cgt1_outsp = incom[sort(shuffle_len[1:cgt1_up])]
    cgt1_insp = incom[sort(shuffle_len[cgt1_up+1:end])]

    # If there is no common factor of product of input and product of output, skip
    if cgt1_up > 1
        vo_out = getNsave_validout(S, Tuple(cgt1_outsp))
        possible_out = Set(vo_out.out_spaces)
    else possible_out = Set(cgt1_outsp[1]) end

    if cgt1_dn > 1
        vo_in = getNsave_validout(S, Tuple(cgt1_insp))
        possible_in = Set(vo_in.out_spaces)
    else possible_in = Set(cgt1_insp[1]) end

    intersect_outsps = intersect(possible_out, possible_in)
    if isempty(intersect_outsps) return nothing end
    cgt1_internalsp = rand(collect(intersect_outsps))
    
    # Number of contracted legs
    nctlegs = rand(1:len)
    # The number of contracted output spaces from CGT1
    nctlegs_out1 = rand(max(nctlegs-cgt1_dn,0):min(nctlegs, cgt1_up))
    # The number of contracted input spaces from CGT1
    nctlegs_in1 = nctlegs - nctlegs_out1

    # Pick 'nctlegs_out1' legs from the output spaces of CGT1
    contracted_1out2in = cgt1_outsp[sort(shuffle(1:cgt1_up)[1:nctlegs_out1])]
    # Pick 'nctlegs_in1' legs from the input spaces of CGT1
    contracted_1in2out = cgt1_insp[sort(shuffle(1:cgt1_dn)[1:nctlegs_in1])]

    # List of contracted spaces
    contracted_spaces = vcat(contracted_1out2in, contracted_1in2out)
    # Product of dimensions of contracted spaces
    contracted_dim = prod([LurCGT.dimension(getNsave_irep(S, BigInt, sp)) 
                            for sp in contracted_spaces])


    # Get the input and output of the second CGT
    # The number of contracted input spaces from CGT2 (=nctlegs_out1)
    nctlegs_in2 = nctlegs_out1
    # The number of contracted output spaces from CGT2 (=nctlegs_in1)
    nctlegs_out2 = nctlegs_in1
    remain_dim = 10 * div(dim_limit, contracted_dim)
    found = false
    find_cnt = 0

    cgt2_insp, cgt2_outsp = [], []
    cgt2_internalsp = nothing

    # Possible minimum length of incomnew2
    len2_min = max(3-nctlegs_in2-nctlegs_out2,
                    nctlegs_in2==0 || nctlegs_out2==0 ? 1 : 0)

    while true
        incomnew2 = sort(generate_incom(S, remain_dim, qlimit))
        len2 = length(incomnew2)
        if len2 < len2_min continue end

        nopenlegs_in2_min = nctlegs_in2>0 ? 0 : 1
        nopenlegs_in2_max = nctlegs_out2>0 ? len2 : len2-1
        # The number of open input legs of CGT2
        nopenlegs_in2 = rand(nopenlegs_in2_min:nopenlegs_in2_max)
        # The number of open output legs of CGT2
        nopenlegs_out2 = len2 - nopenlegs_in2
        @assert nctlegs_in2 + nopenlegs_in2 >= 1 && nctlegs_out2 + nopenlegs_out2 >= 1
        @assert nctlegs_in2 + nopenlegs_in2 + nctlegs_out2 + nopenlegs_out2 >= 3

        shuffle_len2 = shuffle(1:len2)
        # Open input spaces of CGT2
        opensp_in2 = incomnew2[sort(shuffle_len2[1:nopenlegs_in2])]
        # Open output spaces of CGT2
        opensp_out2 = incomnew2[sort(shuffle_len2[nopenlegs_in2+1:end])]

        # Open input, output spaces of CGT2, respectively
        cgt2_insp = sort(vcat(opensp_in2, contracted_1out2in))
        cgt2_outsp = sort(vcat(opensp_out2, contracted_1in2out))

        # Check whether the output spaces are valid
        if length(cgt2_insp) > 1
            vo_in2 = getNsave_validout(S, Tuple(cgt2_insp))
            possible_in = Set(vo_in2.out_spaces)
        else possible_in = Set(cgt2_insp[1]) end

        if length(cgt2_outsp) > 1
            vo_out2 = getNsave_validout(S, Tuple(cgt2_outsp))
            possible_out = Set(vo_out2.out_spaces)
        else possible_out = Set(cgt2_outsp[1]) end

        intersect_outsps = intersect(possible_out, possible_in)
        if !isempty(intersect_outsps) 
            found = true
            cgt2_internalsp = rand(collect(intersect_outsps))
            break 
        end

        if find_cnt > 100 break end
        find_cnt += 1
    end
    if !found return nothing end

    # cgt1_insp, cgt1_outsp, cgt2_insp, cgt2_outsp are now ready
    # also contracted_1in2out and contracted_1out2in
    # also cgt1_internalsp and cgt2_internal
    @assert length(cgt1_insp) >= 1 && length(cgt1_outsp) >= 1
    @assert length(cgt2_insp) >= 1 && length(cgt2_outsp) >= 1
    @assert length(cgt1_insp) + length(cgt1_outsp) >= 3
    @assert length(cgt2_insp) + length(cgt2_outsp) >= 3

    cgt1legs = Int[]
    cgt2legs = Int[]
    contracted_1in2out_dict = Dict{NTuple{NZ, Int}, Int}()
    contracted_1out2in_dict = Dict{NTuple{NZ, Int}, Int}()

    for sp in contracted_1in2out
        if verbose > 1 println(sp) end
        if haskey(contracted_1in2out_dict, sp)
            contracted_1in2out_dict[sp] += 1
        else
            contracted_1in2out_dict[sp] = 1
        end
    end

    for sp in contracted_1out2in
        if verbose > 1 println(sp) end
        if haskey(contracted_1out2in_dict, sp)
            contracted_1out2in_dict[sp] += 1
        else
            contracted_1out2in_dict[sp] = 1
        end
    end

    if verbose > 1
        println("contracted_1in2out: $(contracted_1in2out)")
        println("contracted_1out2in: $(contracted_1out2in)")
        println("contracted_1in2out_dict: $(contracted_1in2out_dict)")
        println("contracted_1out2in_dict: $(contracted_1out2in_dict)")
    end

    for (sp, cnt) in contracted_1in2out_dict
        idx1_in = sample(findall(==(sp), cgt1_insp), cnt; replace=false)
        idx2_out = sample(findall(==(sp), cgt2_outsp), cnt; replace=false)
        append!(cgt1legs, idx1_in.+length(cgt1_outsp))
        append!(cgt2legs, idx2_out)
    end
    for (sp, cnt) in contracted_1out2in_dict
        idx1_out = sample(findall(==(sp), cgt1_outsp), cnt; replace=false)
        idx2_in = sample(findall(==(sp), cgt2_insp), cnt; replace=false)
        append!(cgt1legs, idx1_out)
        append!(cgt2legs, idx2_in.+length(cgt2_outsp))
    end

    if verbose > 1
        println(cgt2_insp)
        println(cgt2_outsp)
    end
    return cgt1_outsp, cgt1_insp, cgt2_outsp, cgt2_insp, cgt1_internalsp, cgt2_internalsp, 
    cgt1legs, cgt2legs, contracted_1out2in, contracted_1in2out
end

# Ultimate test for generating X-symbol
function test_Xsymbol_randinput(::Type{S},
    ::Type{FT},
    dim_limit=100000,
    qlimit=4,
    outlimit=100,
    ninput=40;
    verbose=0) where {S<:NonabelianSymm, FT<:AbstractFloat}

    NZ = nzops(S)
    for i in 1:ninput
        println("Test input #$i")
        spaces = nothing
        while isnothing(spaces) 
            spaces = get_spaces(S, dim_limit, qlimit, outlimit; verbose)
        end

        cgt1_outsp, cgt1_insp, cgt2_outsp, cgt2_insp, cgt1_internalsp, cgt2_internalsp, 
            cgt1legs, cgt2legs, contracted_1out2in, contracted_1in2out = spaces

        #cgt1_outsp = [(2, 1), (2, 3)]
        #cgt1_insp = [(1, 1), (1, 1)]
        #cgt2_outsp = [(1, 1), (3, 2)]
        #cgt2_insp = [(2, 3), (3, 1)]
        #cgt1_internalsp = (1, 1)
        #cgt2_internalsp = (4, 3)
        #cgt1legs = [4, 2]
        #cgt2legs = [1, 3]
        #contracted_1out2in = [(2, 3)]
        #contracted_1in2out = [(1, 1)]


        CGT1up = LurCGT.random_FTree(S, Tuple(cgt1_outsp), (cgt1_internalsp,))
        CGT1down = LurCGT.random_FTree(S, Tuple(cgt1_insp), (cgt1_internalsp,))
        CGT2up = LurCGT.random_FTree(S, Tuple(cgt2_outsp), (cgt2_internalsp,))
        CGT2down = LurCGT.random_FTree(S, Tuple(cgt2_insp), (cgt2_internalsp,))

        CGT1up_arr = FTree2arr(CGT1up, FT, false)
        CGT1down_arr = FTree2arr(CGT1down, FT, false)
        CGT2up_arr = FTree2arr(CGT2up, FT, false)
        CGT2down_arr = FTree2arr(CGT2down, FT, false)

        # Written in old convention except this part
        println("CGT1: $(cgt1_outsp)->$(cgt1_internalsp)->$(cgt1_insp)")
        println("CGT2: $(cgt2_outsp)->$(cgt2_internalsp)->$(cgt2_insp)")
        println("contracted legs: $(cgt1legs) of CGT1, $(cgt2legs) of CGT2")

        CGT1 = contract_arrs(CGT1up_arr, CGT1down_arr, 
                             (length(cgt1_outsp)+1,), (length(cgt1_insp)+1,))
        CGT2 = contract_arrs(CGT2up_arr, CGT2down_arr, 
                             (length(cgt2_outsp)+1,), (length(cgt2_insp)+1,))
        arr1 = contract_arrs(CGT1, CGT2, Tuple(cgt1legs), Tuple(cgt2legs))

        # permute the contraction result
        nout1, nin1 = length(cgt1_outsp), length(cgt1_insp)
        nout2, nin2 = length(cgt2_outsp), length(cgt2_insp)
        if verbose > 1 println("nout1: $nout1, nin1: $nin1, nout2: $nout2, nin2: $nin2") end
        ncont_1out2in = length(contracted_1out2in)
        ncont_1in2out = length(contracted_1in2out)

        if verbose > 1
            println("ncont_1out2in: $ncont_1out2in, ncont_1in2out: $ncont_1in2out")
            println(contracted_1out2in)
            println(contracted_1in2out)
        end

        nout1_open, nin2_open = nout1 - ncont_1out2in, nin2 - ncont_1out2in
        nin1_open, nout2_open = nin1 - ncont_1in2out, nout2 - ncont_1in2out
        nout1, nin1, nout2, nin2 = 
        length(CGT1up.ins), length(CGT1down.ins), 
        length(CGT2up.ins), length(CGT2down.ins)

        if verbose > 1
            println("nout1_open: $nout1_open, nin2_open: $nin2_open")
            println("nin1_open: $nin1_open, nout2_open: $nout2_open")
        end

        perm1 = vcat(collect(1:nout1_open),
                     collect(nout1_open+nin1_open+1:nout1_open+nin1_open+nout2_open),
                     collect(nout1_open+1:nout1_open+nin1_open),
                     collect(nout1_open+nin1_open+nout2_open+1:nout1_open+nin1_open+nout2_open+nin2_open))
        if verbose > 1 println("perm1: $perm1") end
        arr1 = permutedims(arr1, perm1)

        # Finished to get arr1, which is obtained X-symbol from direct contraction of sparse CGTs

        out1_opensp = [CGT1up.ins[i] for i in 1:nout1 if !(i in cgt1legs)]
        in1_opensp = [CGT1down.ins[i-nout1] for i in nout1+1:nout1+nin1 if !(i in cgt1legs)]
        out2_opensp = [CGT2up.ins[i] for i in 1:nout2 if !(i in cgt2legs)]
        in2_opensp = [CGT2down.ins[i-nout2] for i in nout2+1:nout2+nin2 if !(i in cgt2legs)]
        outsp = vcat(out1_opensp, out2_opensp)
        insp = vcat(in1_opensp, in2_opensp)

        out_sortperm, int_sortperm = sortperm(outsp), sortperm(insp)
        perm2 = vcat(out_sortperm, int_sortperm .+ length(outsp))
        arr1 = permutedims(arr1, perm2)
        final_sz = size(arr1)
        if isempty(outsp) final_sz = (1, final_sz...) end
        if isempty(insp) final_sz = (final_sz..., 1) end
        arr1 = reshape(arr1, final_sz...)


        by_frsymbols = LurCGT.contN2canonical(CGT1up, CGT1down, 
        CGT2up, CGT2down, Tuple(cgt1legs), Tuple(cgt2legs); verbose=0)

        arr2 = contract(by_frsymbols, FT; verbose=0)

        ndiff = norm(arr1 - arr2) / norm(arr1)
        cr1 = FT == BigInt ? 1e-68 : 1e-2
        cr2 = FT == BigInt ? 1e-70 : 1e-10
        if norm(arr1) < cr1 && arr1 + arr2 == arr1
            println("Success: arr1 and arr2 are both zero.")
        elseif ndiff < cr2
            println("Success: arr1 and arr2 are the same. with error $(ndiff)")
        else
            error("Failure: arr1 and arr2 are not the same. with error $(ndiff)")
        end
    end
end

# Get the sparse array from final result of getting X-symbol
function contract(step2_result::LurCGT.Step2_result{S, NI, NO, NZ, M},
    ::Type{FT};
    verbose=0) where {S<:NonabelianSymm, NI, NO, NZ, M, FT<:AbstractFloat}
    szs = Int[]
    for insp in step2_result.ins
        irep = getNsave_irep(S, BigInt, insp)
        push!(szs, LurCGT.dimension(irep))
    end
    for outsp in step2_result.outs
        irep = getNsave_irep(S, BigInt, outsp)
        push!(szs, LurCGT.dimension(irep))
    end
    arr = SparseArray(zeros(FT, szs...))

    if NI == 1 && NO == 1
        for intsps in keys(step2_result.coeff)
            coeff, nfac = step2_result.coeff[intsps], step2_result.coeff_nfac[intsps]
            arr += Matrix(I, szs...) * coeff[] / nfac
        end
    elseif NO == 1
        @assert NI >= 2
        for intsps in keys(step2_result.coeff)
            cont_res = contract_CG3s(S, collect(step2_result.ins), step2_result.outs[1], intsps, Val(NI), FT)
            om_arr = SparseArray{FT}(step2_result.coeff[intsps] * step2_result.coeff_nfac[intsps])
            arr += contract_om(cont_res, om_arr, Val(NI))
        end
    elseif NI == 1
        @assert NO >= 2
        for intsps in keys(step2_result.coeff)
            cont_res = contract_CG3s(S, collect(step2_result.outs), step2_result.ins[1], intsps, Val(NO), FT)
            om_arr = SparseArray{FT}(step2_result.coeff[intsps] * step2_result.coeff_nfac[intsps])
            contracted = contract_om(cont_res, om_arr, Val(NO))
            arr += permutedims(contracted, (NO+1, collect(1:NO)...))
        end
    else # Both NI and NO are greater than 1
        @assert NO >= 2 && NI >= 2
        for (i, intsps) in enumerate(keys(step2_result.coeff))
            if verbose > 1 println("$(i)/$(length(keys(step2_result.coeff)))") end
            center_sp = intsps[NI-1]
            cont_res_up = contract_CG3s(S, collect(step2_result.ins), center_sp, intsps[1:NI-2], Val(NI), FT)
            cont_res_dn = contract_CG3s(S, collect(step2_result.outs), center_sp, intsps[NI:end], Val(NO), FT)
            om_arr = SparseArray{FT}(step2_result.coeff[intsps] * step2_result.coeff_nfac[intsps])
            newentry = contract_om_step3(cont_res_up, cont_res_dn, om_arr, Val(NI), Val(NO))
            arr += newentry
        end
    end
    return arr
end

# Test 1j-symbol == CG3 with qlabel = (qlabel, reverse(qlabel), 0)
function test_1j(::Type{S}, qlimit=6, test_inputs=10) where {S<:NonabelianSymm}
    NZ = nzops(S)
    for i in 1:test_inputs
        qlabel = Tuple(rand(0:qlimit) for _=1:NZ)
        dualq = get_dualq(S, qlabel)
        zq = Tuple(0 for _=1:NZ)
        blks_1j, nfac_1j = LurCGT.load_1jblk(S, BigInt, BigInt, qlabel)
        blks_cg3, nfac_cg3 = LurCGT.load_cg3blk(S, BigInt, (qlabel, dualq), [zq])[zq]

        # outer multiplicity = 1
        @assert length(nfac_cg3) == 1
        @assert length(nfac_1j) == 1
        for ((w1, w2, _), blk) in blks_cg3
            #println("for block with weight ($(w1), $(w2)):")
            #println("for cg3:")
            #println(blk[:, :, :, 1])
            #println("normalization factor: $(nfac_cg3[1])")

            #println("for 1j:")
            #println(blks_1j[(w1, w2)])
            #println("normalization factor: $(nfac_1j)")

            if !(blk[:, :, :, 1] == blks_1j[(w1, w2)] && nfac_cg3[1] == nfac_1j[1])
                error("1j-symbol test failed for weight ($(w1), $(w2))")
            end
        end
        println("Test #$i: $(totxt(S)) 1j-symbol test passed for qlabel: $(qlabel)")
    end
end

function test_rsym_from_1j(::Type{S},
    qlimit=4,
    test_inputs=10) where {S<:NonabelianSymm}

    NZ = nzops(S)
    for i in 1:test_inputs
        println("Test input #$i")
        out = Tuple(0 for _=1:NZ)
        if S <: SU
            q = rand(0:qlimit)
            in = Tuple(q for _=1:NZ)
        else
            in = Tuple(rand(0:qlimit) for _=1:NZ)
        end

        blk, fac = LurCGT.load_cg3blk(S, BigInt, (in, in), [out])[out]
        blk_mw = LurCGT.mwpartof(S, blk, (in, in, out), 3)
        om = length(fac)
        rsym_mat_ = zeros(Float64, om, om)

        blk_permuted = LurCGT.permute12(blk_mw)
        LurCGT.conjugate_rsym!(S, blk_mw, in)
        for key in keys(blk_mw)
            if haskey(blk_permuted, key)
                @tensor mat2add[ν, μ] := 
                blk_permuted[key][in1, in2, μ] * blk_mw[key][in1, in2, ν]
                rsym_mat_ += mat2add
            end
        end
        @assert size(rsym_mat_) == (1, 1)
        @assert abs(rsym_mat_[1, 1]) == fac[1].den

        r = getNsave_Rsymbol(S, BigInt, in, out; verbose=0)
        @assert r.nfactor == fac
        @assert r.rsym_mat == rsym_mat_
        println("R-symbol test passed for input: $(in), output: $(out) with value $(rsym_mat_[1, 1]*fac[1])")

    end
end

# This function should be tested when there is no HDF5 file existing.
# Test HDF5 I/O functionality
function test_hdf5_io(::Type{S}) where {S<:NonabelianSymm}
    NZ = nzops(S)
    
    # Test Irep save/load
    println("Testing Irep save/load...")
    test_qlabel = ntuple(i -> i == 1 ? 2 : (i == 2 && NZ >= 2 ? 1 : 0), NZ)
    # Use getNsave to generate and save
    irep_orig = getNsave_irep(S, BigInt, test_qlabel)
    
    # Load directly to verify save worked
    irep_loaded = getNsave_irep(S, BigInt, test_qlabel)
    @assert !isnothing(irep_loaded)
    @assert irep_loaded.qlabel == irep_orig.qlabel
    @assert irep_loaded.dimension == irep_orig.dimension
    @assert length(irep_loaded.Sz) == length(irep_orig.Sz)
    println("  Irep save/load: PASSED")
    
    # Test CG3 save/load
    println("Testing CG3 save/load...")
    q1 = ntuple(i -> i == 2 && NZ >= 2 ? 1 : (i == 1 && NZ == 1 ? 1 : 0), NZ)
    q2 = ntuple(i -> i == 1 ? 1 : 0, NZ)
    # Use getNsave_cg3 from clebsch_io which handles generation
    LurCGT.generate_every_CGT(S, BigInt, BigInt, (q1, q2), nothing; verbose=0)
    
    # Get valid outputs using getNsave
    vo = getNsave_validout(S, (q1, q2))
    @assert !isnothing(vo) && length(vo.out_spaces) > 0
    possible_out = vo.out_spaces[1]
    
    # Use getNsave_cg3 to verify data accessibility
    cg3s_loaded = getNsave_cg3(S, BigInt, (q1, q2), [possible_out])
    @assert !isnothing(cg3s_loaded)
    @assert !isempty(cg3s_loaded[possible_out].blocks)
    @assert length(cg3s_loaded[possible_out].nfactor) > 0
    println("  CG3 save/load: PASSED")
    
    # Test F-symbol save/load
    println("Testing F-symbol save/load...")
    in1, in2, in3 = q1, q2, q1
    # Find a valid output for F-symbol
    vo_f = getNsave_validout(S, sort((in1, in2, in3)))
    @assert !isnothing(vo_f) && length(vo_f.out_spaces) > 0
    out = vo_f.out_spaces[1]
    
    # Use getNsave to generate and save
    fsym = getNsave_Fsymbol(S, BigInt, in1, in2, in3, out; verbose=0)
    @assert !isnothing(fsym)
    @assert fsym.in1 == in1
    @assert fsym.in2 == in2
    @assert fsym.in3 == in3
    @assert fsym.out == out
    
    # Call again to verify load works
    fsym_loaded = getNsave_Fsymbol(S, BigInt, in1, in2, in3, out; verbose=0)
    @assert !isnothing(fsym_loaded)
    @assert size(fsym_loaded.fsym_mat) == size(fsym.fsym_mat)
    println("  F-symbol save/load: PASSED")
    
    # Test R-symbol save/load
    println("Testing R-symbol save/load...")
    in_r = q1
    out_r = q1 .* 2
    # Use getNsave to generate and save
    rsym = getNsave_Rsymbol(S, BigInt, in_r, out_r; verbose=0)
    @assert !isnothing(rsym)
    @assert rsym.in == in_r
    @assert rsym.out == out_r
    
    # Call again to verify load works
    rsym_loaded = getNsave_Rsymbol(S, BigInt, in_r, out_r; verbose=0)
    @assert !isnothing(rsym_loaded)
    @assert size(rsym_loaded.rsym_mat) == size(rsym.rsym_mat)
    println("  R-symbol save/load: PASSED")
    
    # Test OMList save/load
    println("Testing OMList save/load...")
    incom_test = (q1, q2)
    out_test = possible_out
    # Use getNsave to generate and save
    omlist = getNsave_omlist(S, incom_test, out_test)
    @assert !isnothing(omlist)
    @assert omlist.incom_spaces == incom_test
    @assert omlist.out_space == out_test
    
    # Call again to verify load works
    omlist_loaded = getNsave_omlist(S, incom_test, out_test)
    @assert !isnothing(omlist_loaded)
    @assert omlist_loaded.totalOM == omlist.totalOM
    println("  OMList save/load: PASSED")
    
    # Test ValidOuts save/load
    println("Testing ValidOuts save/load...")
    # Use getNsave to generate and save
    vo = getNsave_validout(S, incom_test)
    @assert !isnothing(vo)
    @assert vo.incom_spaces == incom_test
    
    # Call again to verify load works
    vo_loaded = getNsave_validout(S, incom_test)
    @assert !isnothing(vo_loaded)
    @assert length(vo_loaded.out_spaces) == length(vo.out_spaces)
    @assert vo_loaded.out_spaces == vo.out_spaces
    println("  ValidOuts save/load: PASSED")
    
    println("All HDF5 I/O tests passed!")
end

# Test thread-safe HDF5 access
function test_hdf5_threadsafety(::Type{S}) where {S<:NonabelianSymm}
    NZ = nzops(S)
    println("Testing thread-safe concurrent access...")
    
    # Generate some test data first using getNsave
    test_qlabels = [ntuple(i -> i <= j ? 1 : 0, NZ) for j in 0:min(3, NZ)]
    for q in test_qlabels
        getNsave_irep(S, BigInt, q)
    end
    
    # Concurrent reads using getNsave (which will load from cache)
    results = Vector{Any}(undef, length(test_qlabels))
    println(Threads.nthreads(), " threads will be used for concurrent reads.")
    Threads.@threads for i in 1:length(test_qlabels)
        results[i] = getNsave_irep(S, BigInt, test_qlabels[i])
    end
    
    # Verify all reads succeeded
    for (i, irep) in enumerate(results)
        @assert !isnothing(irep)
        @assert irep.qlabel == test_qlabels[i]
    end
    
    println("  Thread-safe concurrent reads: PASSED")
    println("Thread-safety test completed!")
end

function test_pentagon_randinput(::Type{S},
    ::Type{CT},
    qlimit=3,
    test_inputs=10;
    verbose=0) where {S<:NonabelianSymm, CT<:AbstractFloat}

    for i in 1:test_inputs
        println("Pentagon test #$i")
        ins_ = generate_incom(S, 1e10, qlimit)
        @assert length(ins_) >= 4; ins = sort(ins_[1:4])
        vo = getNsave_validout(S, Tuple(ins))
        out = select_out(vo, 300)
        println("Input spaces: $(ins), Output space: $(out)")
        test_pentagon(S, CT, ins, out; verbose)
    end
end

function test_pentagon(::Type{S},
    ::Type{CT},
    ins::Vector{NTuple{NZ, Int}},
    out::NTuple{NZ, Int};
    verbose=0) where {S<:NonabelianSymm, CT<:AbstractFloat, NZ}

    @assert NZ == nzops(S)
    @assert length(ins) == 4
    FTree_path1 = LurCGT.random_FTree(S, Tuple(ins), (out,))
    FTree_path2 = copy(FTree_path1)
    outs = (out,); ins = Tuple(ins)

    # Complete FTree_path1
    FTree_path1_ = LurCGT.apply_Fsymbol(FTree_path1, ins, outs, 1, false, path1s1_iofunc; verbose=0)
    FTree_path1 = change_internal_order(FTree_path1_)
    FTree_path1 = LurCGT.apply_Fsymbol(FTree_path1, ins, outs, 1, false, path1s2_iofunc; verbose=0)

    # Complete FTree_path2
    FTree_path2 = LurCGT.apply_Fsymbol(FTree_path2, ins, outs, 2, false, path2s1_iofunc; verbose=0)
    FTree_path2 = LurCGT.apply_Fsymbol(FTree_path2, ins, outs, 1, false, path2s2_iofunc; verbose=0)
    FTree_path2 = LurCGT.apply_Fsymbol(FTree_path2, ins, outs, 2, false, path2s3_iofunc; verbose=0)

    # Compare them
    compare_FTrees(FTree_path1, FTree_path2; verbose)
    println("Pentagon test passed.")
end

function compare_FTrees(FTree1::FTree{S, 4},
    FTree2::FTree{S, 4};
    verbose=0) where S<:NonabelianSymm

    @assert FTree1.ins == FTree2.ins
    @assert FTree1.outs == FTree2.outs

    for intsps in keys(FTree1.coeff)
        coeff1 = FTree1.coeff[intsps]
        nfac1 = FTree1.coeff_nfac[intsps]
        coeff2 = FTree2.coeff[intsps]
        nfac2 = FTree2.coeff_nfac[intsps]

        @assert size(coeff1) == size(coeff2)
        @assert norm(coeff1 - coeff2) == 0
        @assert nfac1 == nfac2
    end
    for intsps in keys(FTree2.coeff)
        @assert haskey(FTree1.coeff, intsps)
    end
end

function change_internal_order(FTree::FTree{S, 4}) where S<:NonabelianSymm
    NZ = nzops(S)
    ncoeff = Dict{NTuple{2, NTuple{NZ, Int}}, Array{BigInt, 3}}()
    ncoeff_nfac = Dict{NTuple{2, NTuple{NZ, Int}}, Rational{BigInt}}()

    for intsps in keys(FTree.coeff)
        coeff = FTree.coeff[intsps]
        nfac = FTree.coeff_nfac[intsps]
        nkey = (intsps[2], intsps[1])
        ncoeff[nkey] = permutedims(coeff, (1, 3, 2))
        ncoeff_nfac[nkey] = nfac
    end

    new_FTree = LurCGT.create_FTree(S, FTree.ins, FTree.outs, ncoeff, ncoeff_nfac, false)
    return new_FTree
end

function path1s1_iofunc(FTree::FTree{S, 4},
    rem::NTuple{1, NTuple{NZ, Int}},
    i::Int) where {S<:NonabelianSymm, NZ}

    @assert NZ == nzops(S) && i == 1
    return (rem[1], FTree.ins[3], FTree.ins[4], FTree.outs[1])
end

function path1s2_iofunc(FTree::FTree{S, 4},
    rem::NTuple{1, NTuple{NZ, Int}},
    i::Int) where {S<:NonabelianSymm, NZ}

    @assert NZ == nzops(S) && i == 1
    return (FTree.ins[1], FTree.ins[2], rem[1], FTree.outs[1])
end

function path2s1_iofunc(FTree::FTree{S, 4},
    rem::NTuple{1, NTuple{NZ, Int}},
    i::Int) where {S<:NonabelianSymm, NZ}

    @assert NZ == nzops(S) && i == 2
    return (FTree.ins[1], FTree.ins[2], FTree.ins[3], rem[1])
end

function path2s2_iofunc(FTree::FTree{S, 4},
    rem::NTuple{1, NTuple{NZ, Int}},
    i::Int) where {S<:NonabelianSymm, NZ}

    @assert NZ == nzops(S) && i == 1
    return (FTree.ins[1], rem[1], FTree.ins[4], FTree.outs[1])
end

function path2s3_iofunc(FTree::FTree{S, 4},
    rem::NTuple{1, NTuple{NZ, Int}},
    i::Int) where {S<:NonabelianSymm, NZ}

    @assert NZ == nzops(S) && i == 2
    return (FTree.ins[2], FTree.ins[3], FTree.ins[4], rem[1])
end

function test_Xsym_1j(::Type{S},
    dim_limit=100000,
    qlimit=4,
    outlimit=200,
    ninput=20;
    verbose=0) where {S<:NonabelianSymm}

    NZ = nzops(S)
    for i in 1:ninput
        println("Test input #$i")
        spaces = nothing
        while isnothing(spaces)
            spaces = get_spaces(S, dim_limit, qlimit, outlimit; verbose)
        end
        # TODO: currently, get_spaces is written in old convention. Need to fix it later
        # This function is written in new convention
        cgt_in, cgt_out, _, _, _, _, _, _, _, _ = spaces
        NI = length(cgt_in); NO = length(cgt_out)
        cgt_in = Tuple(cgt_in); cgt_out = Tuple(cgt_out)

        cgt_1st = rand(Bool)
        ctleg = (rand(1:NI+NO),); out_contracted = ctleg[1] > NI

        # Contracted space
        ctsp = out_contracted ? cgt_out[ctleg[1]-NI] : cgt_in[ctleg[1]]
        dualq = get_dualq(S, ctsp)

        if ctsp == dualq
            legs_1j = (ctsp, ctsp); ctleg_1j = (rand(1:2),)
        else
            legs_1j = minmax(ctsp, dualq)
            ctleg_1j = ((ctsp == legs_1j[1] ? 1 : 2),)
        end

        println("CGT $(cgt_1st ? 1 : 2) legs: $(cgt_in)->$(cgt_out), $(ctleg[1])th leg $(ctsp) is contracted, 1j: $(legs_1j), $(ctleg_1j[1])th leg is contracted")

        ins_1j = out_contracted ? legs_1j : ()
        outs_1j = out_contracted ? () : legs_1j

        if cgt_1st
            cgt1up, cgt1dn = cgt_in, cgt_out
            cgt2up, cgt2dn = ins_1j, outs_1j
            ctlegs1, ctlegs2 = ctleg, ctleg_1j
        else
            cgt1up, cgt1dn = ins_1j, outs_1j
            cgt2up, cgt2dn = cgt_in, cgt_out
            ctlegs1, ctlegs2 = ctleg_1j, ctleg
        end

        Xsym_1j = getNsave_Xsymbol(S, cgt1up, cgt1dn, cgt2up, cgt2dn, 
        ctlegs1, ctlegs2; verbose=verbose, use1j=true, save=false)

        Xsym_no1j = getNsave_Xsymbol(S, cgt1up, cgt1dn, cgt2up, cgt2dn, 
        ctlegs1, ctlegs2; verbose=verbose, use1j=false, save=false)

        ndiff = norm(Xsym_no1j.xsym_arr - Xsym_1j.xsym_arr) / norm(Xsym_no1j.xsym_arr)
        if ndiff < 1e-10
            println("Success: X-symbol with and without 1j-symbol are the same. with error $(ndiff), norm is $(norm(Xsym_no1j.xsym_arr))")
        else
            error("Failure: X-symbol with and without 1j-symbol are not the same. with error $(ndiff)")
        end
    end
end

function get_final_CGTinfo(::Type{S},
    cgt1_insp::NTuple{NI1, NTuple{NZ, Int}},
    cgt1_outsp::NTuple{NO1, NTuple{NZ, Int}},
    cgt2_insp::NTuple{NI2, NTuple{NZ, Int}},
    cgt2_outsp::NTuple{NO2, NTuple{NZ, Int}},
    cgt1legs::NTuple{M, Int},
    cgt2legs::NTuple{M, Int},
    arr3_cont::SparseArray{Float64}) where {S<:NonabelianSymm, NI1, NO1, NI2, NO2, M, NZ}

    cgt1_inopen = [cgt1_insp[i] for i in 1:NI1 if !(i in cgt1legs)]
    cgt1_outopen = [cgt1_outsp[i] for i in 1:NO1 if !(i+NI1 in cgt1legs)]

    cgt2_inopen = [cgt2_insp[i] for i in 1:NI2 if !(i in cgt2legs)]
    cgt2_outopen = [cgt2_outsp[i] for i in 1:NO2 if !(i+NI2 in cgt2legs)]

    cgt3_insp = vcat(cgt1_inopen, cgt2_inopen)
    cgt3_outsp = vcat(cgt1_outopen, cgt2_outopen)

    n1ino, n1outo = length(cgt1_inopen), length(cgt1_outopen)
    n2ino, n2outo = length(cgt2_inopen), length(cgt2_outopen)

    perm1 = vcat(collect(1:n1ino),
        collect(n1ino+n1outo+1:n1ino+n1outo+n2ino),
        collect(n1ino+1:n1ino+n1outo),
        collect(n1ino+n1outo+n2ino+1:n1ino+n1outo+n2ino+n2outo))

    arr3_cont = permutedims(arr3_cont, perm1)
    in_sortperm, out_sortperm = sortperm(cgt3_insp), sortperm(cgt3_outsp)
    perm2 = vcat(in_sortperm, out_sortperm .+ length(cgt3_insp))
    arr3_cont = permutedims(arr3_cont, perm2)

    cgt3_insp = sort(cgt3_insp); cgt3_outsp = sort(cgt3_outsp)
    if isempty(cgt3_insp) cgt3_insp = [Tuple(0 for _=1:NZ)] end
    if isempty(cgt3_outsp) cgt3_outsp = [Tuple(0 for _=1:NZ)] end
    cgt3_insp_, _ = remove_zeros(S, Tuple(cgt3_insp))
    cgt3_outsp_, _ = remove_zeros(S, Tuple(cgt3_outsp))
    CGT3om = get_CGTom(S, Tuple(cgt3_insp_), Tuple(cgt3_outsp_))

    return cgt3_insp, cgt3_outsp, CGT3om, arr3_cont
end

function test_Xsym(::Type{S},
    dim_limit=100000,
    qlimit=4,
    outlimit=200,
    ninput=20;
    verbose=0) where {S<:NonabelianSymm}

    NZ = nzops(S)
    for i in 1:ninput
        println("Test input #$i")
        spaces = nothing
        while isnothing(spaces)
            spaces = get_spaces(S, dim_limit, qlimit, outlimit; verbose)
        end

        cgt1_insp, cgt1_outsp, cgt2_insp, cgt2_outsp, _, _, 
            cgt1legs, cgt2legs, _, _ = spaces

        # Written in old convention except this part
        println("CGT1: $(cgt1_insp)->$(cgt1_outsp)")
        println("CGT2: $(cgt2_insp)->$(cgt2_outsp)")
        println("contracted legs: $(cgt1legs) of CGT1, $(cgt2legs) of CGT2")

        cgt1_insp = Tuple(cgt1_insp); cgt1_outsp = Tuple(cgt1_outsp)
        cgt2_insp = Tuple(cgt2_insp); cgt2_outsp = Tuple(cgt2_outsp)
        cgt1legs = Tuple(cgt1legs); cgt2legs = Tuple(cgt2legs)

        cgt1_insp_, _ = remove_zeros(S, cgt1_insp)
        cgt1_outsp_, _ = remove_zeros(S, cgt1_outsp)
        cgt2_insp_, _ = remove_zeros(S, cgt2_insp)
        cgt2_outsp_, _ = remove_zeros(S, cgt2_outsp)
        CGTom1 = get_CGTom(S, cgt1_insp_, cgt1_outsp_)
        CGTom2 = get_CGTom(S, cgt2_insp_, cgt2_outsp_)

        om1, om2 = CGTom1.totalOM, CGTom2.totalOM
        canbasis1 = get_canonical_basis(S, cgt1_insp, cgt1_outsp, CGTom1; verbose)
        canbasis2 = get_canonical_basis(S, cgt2_insp, cgt2_outsp, CGTom2; verbose)
        @assert !isempty(canbasis1) && !isempty(canbasis2)
        coeff1, coeff2 = rand(Float64, om1), rand(Float64, om2)

        arr1 = SparseArray(zeros(Float64, size(canbasis1[1])))
        arr2 = SparseArray(zeros(Float64, size(canbasis2[1])))
        for i in 1:om1 arr1 += coeff1[i] * canbasis1[i] end
        for i in 1:om2 arr2 += coeff2[i] * canbasis2[i] end

        # Contract arr1 and arr2
        arr3_cont = contract_arrs(arr1, arr2, cgt1legs, cgt2legs)
        
        Xsym = getNsave_Xsymbol(S, cgt1_insp, cgt1_outsp,
        cgt2_insp, cgt2_outsp, cgt1legs, cgt2legs; verbose, save=true)
        if isnothing(Xsym) || iszero(Xsym.xsym_arr)
            # Check if arr3_cont is zero. If yes, test passed
            if norm(arr3_cont) < 1e-10
                println("Success: both direct contraction and X-symbol result are zero. (norm: $(norm(arr3_cont)))")
                continue
            else
                error("Failure: X-symbol could not be generated, but direct contraction is non-zero. (norm: $(norm(arr3_cont)))")
            end
        end

        cgt3_insp, cgt3_outsp, CGTom3, arr3_cont = get_final_CGTinfo(S, 
        cgt1_insp, cgt1_outsp, cgt2_insp, cgt2_outsp, cgt1legs, cgt2legs, arr3_cont)

        canbasis3 = get_canonical_basis(S, Tuple(cgt3_insp), Tuple(cgt3_outsp), CGTom3)
        om3 = CGTom3.totalOM
        @tensor coeff3[i3] := coeff1[i1] * coeff2[i2] * Xsym.xsym_arr[i1, i2, i3]
        arr3_Xsym = SparseArray(zeros(Float64, size(canbasis3[1])))
        for i in 1:om3 arr3_Xsym += coeff3[i] * canbasis3[i] end

        arr3_cont = reshape(arr3_cont, size(arr3_Xsym)...)

        # Compare arr3_cont and arr3_Xsym
        ndiff = norm(arr3_cont - arr3_Xsym) / norm(arr3_cont)
        if ndiff < 1e-10
            println("Success: direct contraction and X-symbol result are the same. with error $(ndiff), norm is $(norm(arr3_cont))")
        else
            error("Failure: direct contraction and X-symbol result are not the same. with error $(ndiff)")
        end
    end
end

function test_CGTperm(::Type{S},
    dim_limit=500000,
    qlimit=3,
    ninput=10;
    verbose=0) where {S<:NonabelianSymm}

    for i in 1:ninput
        NZ = nzops(S)

        nrepeat, qrep = nothing, nothing
        cgt_insp, cgt_outsp = nothing, nothing
        repeat_in = rand(Bool)
        CGTom = nothing

        # Get incoming & outgoing spaces for CGT until found
        while true
            nrepeat = rand(2:4)
            qrep = Tuple(rand(0:qlimit) for i in 1:NZ)
            rep_dim = LurCGT.dimension(getNsave_irep(S, BigInt, qrep))

            rem_dim = div(dim_limit, rep_dim^nrepeat)
            if rem_dim < 100 continue end

            cnt = 0; found = false
            while cnt < 20
                cnt += 1
                spaces = generate_incom(S, rem_dim, qlimit)
                nsp = length(spaces)
                if length(spaces) < 1 continue end

                shuffle_len = shuffle(1:nsp)
                nsp_otherside = rand(1:nsp)
                repside = spaces[shuffle_len[nsp_otherside+1:end]]
                othside = spaces[shuffle_len[1:nsp_otherside]]

                repside = sort(vcat(repside, [qrep for _=1:nrepeat]))
                othside = sort(othside)

                if repeat_in
                    cgt_insp, cgt_outsp = Tuple(repside), Tuple(othside)
                else
                    cgt_insp, cgt_outsp = Tuple(othside), Tuple(repside)
                end
                
                cgt_insp_, _ = remove_zeros(S, cgt_insp)
                cgt_outsp_, _ = remove_zeros(S, cgt_outsp)
                CGTom = get_CGTom(S, cgt_insp_, cgt_outsp_)
                if CGTom.totalOM > 0 found = true; break end
            end
            if found break end
        end
        println("Test input #$i: CGT $(cgt_insp)->$(cgt_outsp) with repeated qlabel $(qrep)")
        NI = length(cgt_insp); NO = length(cgt_outsp)

        repeat_sps = repeat_in ? cgt_insp : cgt_outsp
        i1 = findfirst(x -> x == qrep, repeat_sps)
        i2 = findlast(x -> x == qrep, repeat_sps)
        @assert i2 >= i1 + 1
        perm_repeat = nothing
        while true
            perm_repeat = shuffle(i1:i2)
            if !issorted(perm_repeat) break end
        end
        perm_repeat = vcat(collect(1:i1-1), perm_repeat, collect(i2+1:length(repeat_sps)))
        perm_other = collect(1:(repeat_in ? NO : NI))

        total_perm = repeat_in ? (perm_repeat..., (perm_other .+ NI)...) : (perm_other..., (perm_repeat .+ NI)...)
        canbasis = get_canonical_basis(S, cgt_insp, cgt_outsp, CGTom; verbose)
        @assert !isempty(canbasis); om = CGTom.totalOM

        coeff_before = rand(Float64, om)
        arr_before = SparseArray(zeros(Float64, size(canbasis[1])))
        for i in 1:om arr_before += coeff_before[i] * canbasis[i] end

        arr_perm = permutedims(arr_before, total_perm)

        CGTperm = getNsave_CGTperm(S, cgt_insp, cgt_outsp, total_perm; save=true)
        println(total_perm)
        if isnothing(CGTperm)
            println("Permutation becomes identity")
            arr_CGTperm = arr_before
        else
            coeff_after = CGTperm.perm_arr * coeff_before
            arr_CGTperm = SparseArray(zeros(Float64, size(canbasis[1])))
            for i in 1:om arr_CGTperm += coeff_after[i] * canbasis[i] end
        end

        # Check whether arr_perm and arr_CGTperm are the same
        ndiff = norm(arr_perm - arr_CGTperm) / norm(arr_perm)
        if ndiff < 1e-10
            println("Success: direct permutation and CGTperm result are the same. with error $(ndiff), norm is $(norm(arr_perm))")
        else
            error("Failure: direct permutation and CGTperm result are not the same. with error $(ndiff), norm of correct result is $(norm(arr_perm))")
        end
    end
end

function test_CGT_conj(::Type{S},
    dim_limit=100000,
    qlimit=4,
    outlimit=200,
    ninput=20;
    verbose=0) where {S<:NonabelianSymm}

    NZ = nzops(S)
    for i in 1:ninput
        spaces = nothing
        while isnothing(spaces) || spaces[1] == spaces[2]
            spaces = get_spaces(S, dim_limit, qlimit, outlimit; verbose)
        end
        # TODO: currently, get_spaces is written in old convention. Need to fix it later
        # This function is written in new convention
        cgt_in, cgt_out, _, _, _, _, _, _, _, _ = spaces
        NI = length(cgt_in); NO = length(cgt_out)
        cgt_in = Tuple(cgt_in); cgt_out = Tuple(cgt_out)

        println("Test input #$i: CGT $(cgt_in)->$(cgt_out)")

        CGTom = get_CGTom(S, cgt_in, cgt_out)
        om = CGTom.totalOM
        CGTom_conj = get_CGTom(S, cgt_out, cgt_in)
        @assert CGTom_conj.totalOM == om

        canbasis = get_canonical_basis(S, cgt_in, cgt_out, CGTom; verbose)
        canbasis_conj = get_canonical_basis(S, cgt_out, cgt_in, CGTom_conj; verbose)

        perm = (NI+1:NI+NO..., 1:NI...)
        for j in 1:om
            diffnorm = norm(permutedims(canbasis[j], perm) - canbasis_conj[j])
            @assert diffnorm < 1e-10 "Difference in conjugate basis is too large: $diffnorm"
        end
        println("Test passed")
    end
end

function test_CGT_conj_sameq(::Type{S},
    dim_limit=1000,
    qlimit=4,
    outlimit=200,
    ninput=20;
    verbose=0) where {S<:NonabelianSymm}

    NZ = nzops(S)
    for i in 1:ninput
        sps = Tuple(sort(generate_incom(S, dim_limit, qlimit)))
        N = length(sps)
        println("Test input #$i: CGT $(sps) -> $(sps)")

        CGTom = get_CGTom(S, sps, sps)
        om = CGTom.totalOM

        canbasis = get_canonical_basis(S, sps, sps, CGTom; verbose)
        perm_CGT = (N+1:2N..., 1:N...)
        perm_om = LurCGT.get_conj_perm(CGTom)

        for j in 1:om
            diffnorm = norm(permutedims(canbasis[j], perm_CGT) - canbasis[perm_om[j]])
            @assert diffnorm < 1e-10 "Difference in conjugate basis is too large: $diffnorm"
        end
        println("test passed")
    end
end

