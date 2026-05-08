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

