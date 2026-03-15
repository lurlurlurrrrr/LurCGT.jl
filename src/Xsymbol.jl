# The file containing the procedure to get the X-symbol
# U1: The number of upper legs of the first CGT
# D1: The number of lower legs of the first CGT, similarly for U2, D2
# This is expressed as 64-bit float array. Every element is # obtained 
# from sqrt(p/q), so it has enough precision regardless of its size.
struct Xsymbol{S<:NonabelianSymm, U1, D1, U2, D2, NZ, M}
    xsym_arr::Array{Float64, 3}
    up1sp::NTuple{U1, NTuple{NZ, Int}}
    dn1sp::NTuple{D1, NTuple{NZ, Int}}
    up2sp::NTuple{U2, NTuple{NZ, Int}}
    dn2sp::NTuple{D2, NTuple{NZ, Int}}
    legs1::NTuple{M, Int}
    legs2::NTuple{M, Int}
    size_byte::Int

    function Xsymbol{S, U1, D1, U2, D2, NZ, M}(xsym_arr, up1sp, dn1sp, up2sp, dn2sp, legs1, legs2, size_byte::Int=0) where {S<:NonabelianSymm, U1, D1, U2, D2, NZ, M}
        if size_byte == 0
            obj = new{S, U1, D1, U2, D2, NZ, M}(xsym_arr, up1sp, dn1sp, up2sp, dn2sp, legs1, legs2, 0)
            size_byte = Base.summarysize(obj)
        end
        new{S, U1, D1, U2, D2, NZ, M}(xsym_arr, up1sp, dn1sp, up2sp, dn2sp, legs1, legs2, size_byte)
    end
end

struct CGTom{S<:NonabelianSymm, NZ}
    totalOM::Int
    # The list of possible central spaces
    spaces::Vector{NTuple{NZ, Int}}
    # Accumulative sum of outer multiplicities
    om_accumul::Vector{Int}
    # Possible outer multiplicities of FTrees for each central space
    FTree_oms::Vector{Tuple{Int, Int}}
    # True if the outgoing leg is 'larger', which means that
    # # of outgoing legs > # of incoming legs or 
    # outgoing == incoming and qlabels of outgoing is larger 
    outgoing_larger::Bool
end

# mfac: If there is an integer multiplicative factor, it is stored here.
# NI: The number of open legs of upper side FTree. (except upper/lower_contracted) 
# NO: The number of open legs of lower side FTree.
# upper_contracted: The legs to be contracted with the other Step1_result struct
# similar for lower_contracted. upper_contracted is the output leg of upper FTree
# M: The total number of intermediate legs
struct Step1_result{S<:NonabelianSymm, NI, NO, NZ, M} <: AbstractCG3contract{S, NI, NO, NZ}
    # Incoming arrows of the upper FTree.
    ins::NTuple{NI, NTuple{NZ, Int}}
    # Outgoing arrows of the lower FTree.
    outs::NTuple{NO, NTuple{NZ, Int}}
    # Legs to be contracted with the other Step1_result struct
    upper_contracted::NTuple{NZ, Int}
    lower_contracted::NTuple{NZ, Int}
    coeff::Dict{NTuple{M, NTuple{NZ, Int}}, Array{BigInt}}
    coeff_nfac::Dict{NTuple{M, NTuple{NZ, Int}}, Rational{BigInt}}
end

# Contraction of 4 CG3s appear in step 2. 
# It is eventually reduced to contraction of two CG3s.
struct FourCG3s{S<:NonabelianSymm, NZ} <: AbstractCG3contract{S, 2, 2, NZ}
    # It has only two incoming and two outgoing arrows
    ins::NTuple{2, NTuple{NZ, Int}}
    outs::NTuple{2, NTuple{NZ, Int}}
    coeff::Dict{NTuple{4, NTuple{NZ, Int}}, Array{BigInt, 4}}
    coeff_nfac::Dict{NTuple{4, NTuple{NZ, Int}}, Rational{BigInt}}
end

# The result of step 2. With appropriate action of F-symbol, it becomes canonical form
struct Step2_result{S<:NonabelianSymm, NI, NO, NZ, M, D} <: AbstractCG3contract{S, NI, NO, NZ}
    # Incoming arrows at the lower side
    ins::NTuple{NI, NTuple{NZ, Int}}
    # Outgoing arrows at the upper side
    outs::NTuple{NO, NTuple{NZ, Int}}
    # M = NI + NO - 3
    coeff::Dict{NTuple{M, NTuple{NZ, Int}}, Array{BigInt, D}}
    coeff_nfac::Dict{NTuple{M, NTuple{NZ, Int}}, Rational{BigInt}}
end

get_Xsym_path(::Type{S}, ::Type{CT}) where {S<:NonabelianSymm, CT<:Integer} =
    "$(homedir())/LurCGT/Xsymbols/$(totxt(S))_$(CT)"

to_txt(tup::NTuple{N, NTuple{NZ, Int}}) where {N, NZ} =
    join(["$(to_txt(sp))" for sp in tup], "_")

function get_Xsym_fname(::Type{S},
    up1sp::NTuple{U1, NTuple{NZ, Int}},
    dn1sp::NTuple{D1, NTuple{NZ, Int}},
    up2sp::NTuple{U2, NTuple{NZ, Int}},
    dn2sp::NTuple{D2, NTuple{NZ, Int}},
    ctlegs1::NTuple{M, Int},
    ctlegs2::NTuple{M, Int}) where {S<:NonabelianSymm, U1, D1, U2, D2, NZ, M}

    @assert NZ == nzops(S)
    fpath = get_Xsym_path(S, BigInt)
    mkpath(fpath)
    qlabels_txt = "$(to_txt(dn1sp))→$(to_txt(up1sp));__$(to_txt(ctlegs1))__;$(to_txt(dn2sp))→$(to_txt(up2sp));__$(to_txt(ctlegs2))__"
    return joinpath(fpath, "$(qlabels_txt)")
end

remove_zeros(::Type{S}, tup::Tuple{}) where {S<:Symmetry} = (Tuple(0 for _ in 1:nzops(S)),), 0

# Remove all zero spaces, and if there is no non-zero space, add one zero space
function remove_zeros(::Type{S}, 
    tup::NTuple{N, NTuple{NZ, Int}}) where {S<:Symmetry, N, NZ}

    zerotup = Tuple(0 for _ in 1:NZ)
    zcnt = 0
    res = NTuple{NZ, Int}[]
    for sp in tup if sp != zerotup push!(res, sp) else zcnt += 1 end end
    if isempty(res) push!(res, zerotup) end
    return Tuple(res), zcnt
end

function standardize_spaces_and_legs(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    legs::NTuple{M, Int},
    upperfirst::Bool) where {S<:NonabelianSymm, U, D, NZ, M}
    # 'legs' should be sorted in ascending order
    up_contracted = Tuple(i for i in legs if i <= U)
    dn_contracted = Tuple(i-U for i in legs if i > U)

    @assert issorted(upsp) && issorted(dnsp)
    upsp, nzup = remove_zeros(S, upsp)
    dnsp, nzdn = remove_zeros(S, dnsp)

    upc = Tuple(i-nzup for i in up_contracted if i > nzup)
    dnc = Tuple(i-nzdn+length(upsp) for i in dn_contracted if i > nzdn)
    legs = upperfirst ? (upc..., dnc...) : (dnc..., upc...)
    return upsp, dnsp, legs
end

# use1j: Keyword argument for test purpose only
function getNsave_Xsymbol(::Type{S},
    up1sp::NTuple{U1, NTuple{NZ, Int}},
    dn1sp::NTuple{D1, NTuple{NZ, Int}},
    up2sp::NTuple{U2, NTuple{NZ, Int}},
    dn2sp::NTuple{D2, NTuple{NZ, Int}},
    ctlegs1::NTuple{M, Int},
    ctlegs2::NTuple{M, Int};
    verbose=0,
    use1j=true,
    save=true) where {S<:NonabelianSymm, U1, D1, U2, D2, NZ, M}

    up1sp, dn1sp, ctlegs1 = standardize_spaces_and_legs(S, up1sp, dn1sp, ctlegs1, true)
    up2sp, dn2sp, ctlegs2 = standardize_spaces_and_legs(S, up2sp, dn2sp, ctlegs2, false)
    perm = sortperm(collect(ctlegs1))
    ctlegs1 = Tuple(ctlegs1[i] for i in perm); ctlegs2 = Tuple(ctlegs2[i] for i in perm)
    getNsave_Xsymbol_zeroadded(S, up1sp, dn1sp, up2sp, dn2sp, ctlegs1, ctlegs2; verbose, use1j, save)
end

function detect_1j(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    legs::NTuple{M, Int}=(2,)) where {S<:NonabelianSymm, U, D, NZ, M}

    @assert NZ == nzops(S)
    zerotup = Tuple(0 for _ in 1:NZ)
    if M != 1 return false end
    if U == 1 && D == 2 && upsp[1] == zerotup
        if dnsp[2] != get_dualq(S, dnsp[1]) return false end
        if legs[1] == 2 || legs[1] == 3 return true end
    elseif U == 2 && D == 1 && dnsp[1] == zerotup
        if upsp[1] != get_dualq(S, upsp[2]) return false end
        if legs[1] == 1 || legs[1] == 2 return true end
    end
    return false
end

# Get X-symbol from given space configuration
# TODO: Consider two 1j-symbols are contracted
function getNsave_Xsymbol_zeroadded(::Type{S},
    up1sp::NTuple{U1, NTuple{NZ, Int}},
    dn1sp::NTuple{D1, NTuple{NZ, Int}},
    up2sp::NTuple{U2, NTuple{NZ, Int}},
    dn2sp::NTuple{D2, NTuple{NZ, Int}},
    ctlegs1::NTuple{M, Int},
    ctlegs2::NTuple{M, Int};
    verbose,
    use1j,
    save) where {S<:NonabelianSymm, U1, D1, U2, D2, NZ, M}

    @assert NZ == nzops(S)
    # Try to load from HDF5, use it if exists
    loaded = load_Xsymbol_sqlite(S, up1sp, dn1sp, up2sp, dn2sp, ctlegs1, ctlegs2)
    if !isnothing(loaded) return loaded end

    # First, we need to get outer multiplicities of CGTs involved
    is1j_1 = detect_1j(S, up1sp, dn1sp, ctlegs1)
    CGT1_oms = get_CGTom(S, up1sp, dn1sp, is1j_1)
    @assert CGT1_oms.totalOM > 0 "CGT1 has zero outer multiplicity"
    is1j_2 = detect_1j(S, up2sp, dn2sp, ctlegs2)
    CGT2_oms = get_CGTom(S, up2sp, dn2sp, is1j_2)
    @assert CGT2_oms.totalOM > 0 "CGT2 has zero outer multiplicity"
    
    # Upper / lower spaces for the resulting CGT
    up3sp, dn3sp = get_resulting_spaces(
        up1sp, dn1sp, up2sp, dn2sp, ctlegs1, ctlegs2)
    up3sp, _ = remove_zeros(S, up3sp); dn3sp, _ = remove_zeros(S, dn3sp)
    # TODO: Test when the resulting CGT is 1j-symbol
    is1j_res = detect_1j(S, up3sp, dn3sp)
    CGT3_oms = get_CGTom(S, up3sp, dn3sp, is1j_res)
    if CGT3_oms.totalOM == 0 return end

    is1j_2 = detect_1j(S, up2sp, dn2sp, ctlegs2)
    if (is1j_1 || is1j_2) && use1j
        Xsym_arr = getNsave_Xsymbol_1j(S, up1sp, dn1sp, up2sp, dn2sp, ctlegs1, ctlegs2, is1j_1, is1j_2, CGT3_oms; verbose)
    else
        # bp: before preprocessing
        CGT1_FTrees_up_bp = Dict{NTuple{NZ, Int}, Vector{FTree{S, U1, NZ}}}()
        CGT1_FTrees_dn_bp = Dict{NTuple{NZ, Int}, Vector{FTree{S, D1, NZ}}}()
        CGT2_FTrees_up_bp = Dict{NTuple{NZ, Int}, Vector{FTree{S, U2, NZ}}}()
        CGT2_FTrees_dn_bp = Dict{NTuple{NZ, Int}, Vector{FTree{S, D2, NZ}}}()
        # Prepare unit FTrees, csp: central space
        for csp1 in CGT1_oms.spaces
            fill_FTrees!(CGT1_FTrees_up_bp, up1sp, csp1)
            fill_FTrees!(CGT1_FTrees_dn_bp, dn1sp, csp1)
        end

        for csp2 in CGT2_oms.spaces
            fill_FTrees!(CGT2_FTrees_up_bp, up2sp, csp2)
            fill_FTrees!(CGT2_FTrees_dn_bp, dn2sp, csp2)
        end

        # Step 1-1, first permute legs & apply F-symbols
        # l1up, l1dn, l2up, l2dn: which legs are cocntracted for each FTree
        l1up, l1dn, l2up, l2dn = separate_legs(ctlegs1, ctlegs2, U1, D1, U2, D2)
        @assert length(l1up) == length(l2dn); MU = length(l1up)
        @assert length(l1dn) == length(l2up); MD = length(l1dn)
        for i in 1:MU @assert up1sp[l1up[i]] == dn2sp[l2dn[i]] end
        for i in 1:MD @assert dn1sp[l1dn[i]] == up2sp[l2up[i]] end

        CGT1_FTrees_up, add0_up1 = preprocess_FTrees(CGT1_FTrees_up_bp, l1up)
        CGT1_FTrees_dn, add0_dn1 = preprocess_FTrees(CGT1_FTrees_dn_bp, l1dn)
        CGT2_FTrees_up, add0_up2 = preprocess_FTrees(CGT2_FTrees_up_bp, l2up)
        CGT2_FTrees_dn, add0_dn2 = preprocess_FTrees(CGT2_FTrees_dn_bp, l2dn)
        nz_in, nz_out = add0_up1 + add0_up2, add0_dn1 + add0_dn2

        # Step 1-2, merge the contracted legs + apply R-symbol
        step1_up = combine_FTrees(CGT2_FTrees_dn, CGT1_FTrees_up, true, Val(max(MU, 1)))
        step1_dn = combine_FTrees(CGT1_FTrees_dn, CGT2_FTrees_up, false, Val(max(MD, 1)))

        # Do Step 2 for all combinations of outer multiplicities
        step2_res, nr_input, nr_output = do_step2_forall(step1_up, step1_dn, CGT1_oms, CGT2_oms; verbose=0)

        # Return the X-symbol array before normalized
        Xsym_arr = zeros(Float64, CGT1_oms.totalOM, CGT2_oms.totalOM, CGT3_oms.totalOM)
        # Do Step 3 for all step2 results obtained before
        # This can be parallelized easily
        do_step3_forall!(Xsym_arr, step2_res, nr_input, nr_output, nz_in, nz_out, CGT3_oms; verbose)
    end
    
    if !is1j_1 || !is1j_2
        # The last step: normalize everything
        # 1. Normalize all CGCs appear in the calculation
        # 2. Normalize the canonical basis elements
        # They are form of 1/sqrt(N)
        CGT1_norms = get_CGT_norms(S, up1sp, dn1sp, CGT1_oms, is1j_1; use1j)
        CGT2_norms = get_CGT_norms(S, up2sp, dn2sp, CGT2_oms, is1j_2; use1j)
        CGT3_norms = get_CGT_norms(S, up3sp, dn3sp, CGT3_oms, is1j_res; use1j)

        # multiply by CGT1_norms & CGT2_norms, divide by CGT3_norms
        mul_along_dim!(Xsym_arr, CGT1_norms, 1)
        mul_along_dim!(Xsym_arr, CGT2_norms, 2)
        div_along_dim!(Xsym_arr, CGT3_norms, 3)
    end

    xsym_obj = Xsymbol{S, U1, D1, U2, D2, NZ, M}(Xsym_arr, up1sp, dn1sp, up2sp, dn2sp, ctlegs1, ctlegs2)
    if save save_Xsymbol_sqlite(xsym_obj) end
    return xsym_obj
end

function get_CGT_norms(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    CGT_oms::CGTom{S, NZ},
    is1j::Bool;
    use1j) where {S<:NonabelianSymm, U, D, NZ}

    if is1j && use1j 
        q = U==2 ? upsp[1] : dnsp[1]
        return sqrt.(getNsave_1jsym(S, BigInt, BigInt, q).nfactor)
    end
    om = CGT_oms.totalOM
    norms = Vector{Float64}(undef, om)
    for (i, csp) in enumerate(CGT_oms.spaces)
        csp_dim = dimension(getNsave_irep(S, BigInt, csp))
        start_idx = CGT_oms.om_accumul[i]
        omup, omdn = CGT_oms.FTree_oms[i]

        factors_up = get_CGC_factors(S, upsp, csp)
        factors_dn = get_CGC_factors(S, dnsp, csp)
        partial_mat = fill(1/sqrt(csp_dim), omup, omdn)
        mul_along_dim!(partial_mat, factors_up, 1)
        mul_along_dim!(partial_mat, factors_dn, 2)
        norms[start_idx:start_idx+omup*omdn-1] = 
        CGT_oms.outgoing_larger ? partial_mat[:] : partial_mat'[:]
    end
    return norms
end

function get_CGC_factors(::Type{S},
    sp_list::NTuple{N, NTuple{NZ, Int}},
    csp::NTuple{NZ, Int}) where {S<:NonabelianSymm, N, NZ}

    if N == 1 return [1.0] end
    omlist = getNsave_omlist(S, sp_list, csp)
    factors = Vector{Float64}(undef, omlist.totalOM)
    for (i, isps) in enumerate(omlist.interm_spaces)
        sidx = omlist.cumul[i]
        partial_fac = ones(Float64, omlist.cg3_oms[i]...)
        @assert length(size(partial_fac)) == N-1
        ncg3 = length(omlist.cg3_oms[i])
        for j in 1:ncg3
            cg3in1 = j == ncg3 ? sp_list[1] : isps[j]
            cg3in2 = sp_list[N+1-j]
            cg3in3 = j == 1 ? csp : isps[j-1]
            _, fac = load_cg3blk(S, BigInt, (cg3in1, cg3in2), [cg3in3])[cg3in3]
            mul_along_dim!(partial_fac, [Float64(sqrt(f)) for f in fac], j)
        end
        l = length(partial_fac)
        factors[sidx:sidx+l-1] = partial_fac[:]
    end
    return factors
end

function mul_along_dim!(A, v, dim)
    shape = ntuple(d -> d == dim ? length(v) : 1, ndims(A))
    A .*= reshape(v, shape)
end

function div_along_dim!(A, v, dim)
    shape = ntuple(d -> d == dim ? length(v) : 1, ndims(A))
    A ./= reshape(v, shape)
end



function fill_FTrees!(CGT_FTrees::Dict{NTuple{NZ, Int}, Vector{FTree{S, N, NZ}}},
    incom_sps::NTuple{N, NTuple{NZ, Int}},
    csp::NTuple{NZ, Int}) where {S<:NonabelianSymm, N, NZ}

    omlist = getNsave_omlist(S, incom_sps, csp)
    om = omlist.totalOM
    CGT_FTrees[csp] = Vector{FTree{S, N, NZ}}(undef, om)
    for ii in 1:om CGT_FTrees[csp][ii] = create_unit_FTree(omlist, ii) end
end

# ci: contracted legs indices
function preprocess_FTrees(FTree_dict::Dict{NTuple{NZ, Int}, Vector{FTree{S, N, NZ}}},
    ci::NTuple{M, Int}) where {S<:NonabelianSymm, N, M, NZ}

    Nn = (M == N || M == 0) ? N+1 : N # New number of legs after adding singlet dims
    FTrees = Dict{NTuple{NZ, Int}, Vector{FTree{S, Nn, NZ}}}()
    perm = get_perm_step1(N, ci)
    for (sp, ftree_vec) in FTree_dict
        FTrees[sp] = Vector{FTree{S, Nn, NZ}}(undef, length(ftree_vec))
        for i in 1:length(ftree_vec)
            # First, permute the legs so that the contracted legs are at rightmost
            result = permute!(ftree_vec[i], perm)
            if M == 0 result = add_singdim(result, false; verbose=0) end
            if M == N result = add_singdim(result, true; verbose=0) end

            # Next, apply F-symbols to separate the contracted /legs
            for i=M-1:-1:1
                ft = (a...) -> fsymiofunc_step1(Val(M), a...)
                result = apply_Fsymbol(result, result.ins, result.outs, i, false, ft)
            end
            FTrees[sp][i] = result
        end
    end
    return FTrees, M == N ? 1 : 0
end

const Dict_s1res{S, NI, NO, NZ, Mn} = Dict{Tuple{NTuple{NZ, Int}, NTuple{NZ, Int}}, 
                                       Matrix{Step1_result{S, NI, NO, NZ, Mn}}}

# rontop == true : apply R-symbol to the top CG3 (index 1)
# rontop == false: apply R-symbol to the bottom-side CG3 
function combine_FTrees(FTreesup::Dict{NTuple{NZ, Int}, Vector{FTree{S, N1, NZ}}},
    FTreesdn::Dict{NTuple{NZ, Int}, Vector{FTree{S, N2, NZ}}},
    rontop::Bool,
    ::Val{M}) where {S<:NonabelianSymm, N1, N2, M, NZ}

    Mn = N1 + N2 - 2*M - 1; NI, NO = N2-M, N1-M
    step1_results = Dict_s1res{S, NI, NO, NZ, Mn}()
    ridx = rontop ? 1 : NO + 1
    for (cspup, ftree_vec_up) in FTreesup
        omup = length(ftree_vec_up)
        for (cspdn, ftree_vec_dn) in FTreesdn
            omdn = length(ftree_vec_dn)
            key = (cspup, cspdn)
            step1_results[key] = Matrix{Step1_result{S, NI, NO, NZ, Mn}}(undef, omup, omdn)
            for i in 1:omup
                for j in 1:omdn
                    res = combine_2FTrees(ftree_vec_up[i], ftree_vec_dn[j], Val(M))
                    apply_Rsymbol!(res, ridx, rsymiofunc_prestep2; verbose=0)
                    step1_results[key][i, j] = res
                end
            end
        end
    end
    return step1_results
end


# TODO: Step2 & 3 can run in parallel, make code easier to parallelize
function do_step2_forall(s1up::Dict_s1res{S, NIU, NOU, NZ},
    s1dn::Dict_s1res{S, NID, NOD, NZ},
    CGT1_oms::CGTom{S, NZ},
    CGT2_oms::CGTom{S, NZ}; 
    verbose=0) where {S<:NonabelianSymm, NIU, NOU, NID, NOD, NZ}

    NKL = NIU + NID + NOU + NOD - 3 # New key length
    NDL = NKL + 1 # The dimension of resulting array
    step2_results = Matrix{Step2_result{S, NIU+NID, NOU+NOD, NZ}}(
        undef, CGT1_oms.totalOM, CGT2_oms.totalOM)

    # This can be parallelized easily
    for I in CartesianIndices(step2_results)
        i1, i2 = I.I
        csp1, up1idx, dn1idx = getominfo(CGT1_oms, i1)
        csp2, up2idx, dn2idx = getominfo(CGT2_oms, i2)

        s1res_up = s1up[(csp2, csp1)][dn2idx, up1idx]
        s1res_dn = s1dn[(csp1, csp2)][dn1idx, up2idx]
        step2_results[i1, i2], _, _ = do_step2(s1res_up, s1res_dn; verbose=verbose)
    end

    return step2_results, NID, NOU
end

function do_step3_forall!(Xsym_arr::Array{Float64, 3},
    step2_res::Matrix{Step2_result{S, NI, NO, NZ}},
    nr_input::Int,
    nr_output::Int,
    nz_in::Int,
    nz_out::Int,
    CGT3_oms::CGTom{S, NZ}; 
    verbose=0) where {S<:NonabelianSymm, NI, NO, NZ}

    @assert size(Xsym_arr, 1) == size(step2_res, 1)
    @assert size(Xsym_arr, 2) == size(step2_res, 2)
    @assert size(Xsym_arr, 3) == CGT3_oms.totalOM
    if verbose > 1 println(nz_in, " ", nz_out) end

    # This can be parallelized easily
    for I in CartesianIndices(step2_res)
        i1, i2 = I.I
        step2res = step2_res[i1, i2]
        s3_res = do_step3(step2res, nr_input, nr_output, nz_in, nz_out; verbose=verbose)
        # Convert a result of step3 into vector form
        Xsym_arr[i1, i2, :] = to_vector(s3_res, CGT3_oms; verbose=verbose)
    end
    return Xsym_arr
end

function get_center_space(sps::NTuple{M, NTuple{NZ, Int}},
    ins::NTuple{NI, NTuple{NZ, Int}},
    outs::NTuple{NO, NTuple{NZ, Int}}) where {M, NI, NO, NZ}

    @assert NI >= 1 && NO >= 1
    if NO == 1 return NI-1, outs[1] 
    elseif NI == 1 return 0, ins[1] 
    else return NI-1, sps[NI-1] end
end

function to_vector(s3_res::Step2_result{S, NI, NO, NZ},
    CGT_oms::CGTom{S, NZ};
    verbose=0) where {S<:NonabelianSymm, NI, NO, NZ}

    totalom = CGT_oms.totalOM
    vec = zeros(Float64, totalom)
    mats_dict = Dict{NTuple{NZ, Int}, Matrix{Float64}}()
    omlist_dict_up = Dict{NTuple{NZ, Int}, OMList{S, NI, NZ}}()
    omlist_dict_dn = Dict{NTuple{NZ, Int}, OMList{S, NO, NZ}}()
    for (sps, arr) in s3_res.coeff
        divfac = s3_res.coeff_nfac[sps]
        csp_idx, csp = get_center_space(sps, s3_res.ins, s3_res.outs)
        if !haskey(omlist_dict_up, csp)
            omlist_dict_up[csp] = getNsave_omlist(S, s3_res.ins, csp)
            @assert !haskey(omlist_dict_dn, csp)
            omlist_dict_dn[csp] = getNsave_omlist(S, s3_res.outs, csp)
            @assert !haskey(mats_dict, csp)
            mats_dict[csp] = zeros(Float64, 
                omlist_dict_up[csp].totalOM, omlist_dict_dn[csp].totalOM)
        end
        up_idx, dn_idx = sps[1:csp_idx-1], sps[csp_idx+1:end]

        omlist_upidx = searchsortedfirst(
            omlist_dict_up[csp].interm_spaces, up_idx)
        omlist_dnidx = searchsortedfirst(
            omlist_dict_dn[csp].interm_spaces, dn_idx)
        up_sidx = omlist_dict_up[csp].cumul[omlist_upidx]
        dn_sidx = omlist_dict_dn[csp].cumul[omlist_dnidx]
        up_size = prod(size(arr)[1:NI-1])
        dn_size = prod(size(arr)[NI:end])

        mats_dict[csp][up_sidx:up_sidx+up_size-1, dn_sidx:dn_sidx+dn_size-1] += 
            reshape(arr, up_size, dn_size) * Float64(divfac)
    end
    # Move matrix elements to vector
    for (i, csp) in enumerate(CGT_oms.spaces)
        if !haskey(mats_dict, csp) continue end
        partial_mat = mats_dict[csp]; l = length(partial_mat)
        si = CGT_oms.om_accumul[i]
        vec[si:si+l-1] = CGT_oms.outgoing_larger ? partial_mat[:] : partial_mat'[:]
    end
    return vec
end

function getominfo(CGT_oms::CGTom{S, NZ}, i::Int) where {S<:NonabelianSymm, NZ}
    @assert i >= 1 && i <= CGT_oms.totalOM
    spidx = searchsortedlast(CGT_oms.om_accumul, i)
    csp = CGT_oms.spaces[spidx]
    start_index = CGT_oms.om_accumul[spidx]
    (om_up, om_dn) = CGT_oms.FTree_oms[spidx]
    rel_idx = i - start_index
    # dnidx moves slower from 0 to om_dn-1, upidx moves faster (0 to om_up-1 om_dn times)
    if CGT_oms.outgoing_larger dnidx, upidx = divrem(rel_idx, om_up)
    else upidx, dnidx = divrem(rel_idx, om_dn) end
    @assert dnidx < om_dn && upidx < om_up
    return csp, upidx+1, dnidx+1
end

is_outgoing_larger(upsp::NTuple{N, NTuple{NZ, Int}},
dnsp::NTuple{N, NTuple{NZ, Int}}) where {N, NZ} = upsp < dnsp

is_outgoing_larger(::NTuple{U, NTuple{NZ, Int}},
::NTuple{D, NTuple{NZ, Int}}) where {U, D, NZ} = U < D

function get_CGTom(::Type{S},
    up1sp::NTuple{U, NTuple{NZ, Int}},
    dn1sp::NTuple{D, NTuple{NZ, Int}},
    is1j::Bool=false) where {S<:NonabelianSymm, U, D, NZ}

    @assert issorted(up1sp) && issorted(dn1sp)
    if is1j
        zerotup = Tuple(0 for _ in 1:NZ)
        return CGTom{S, NZ}(1, [zerotup], [1, 2], [(1, 1)], false)
    end
    # 1. Get candidate for central space
    validout_up = getNsave_validout(S, Tuple(sort(collect(up1sp)))).out_spaces
    validout_dn = getNsave_validout(S, Tuple(sort(collect(dn1sp)))).out_spaces
    central_spaces_lst = sort(intersect(validout_up, validout_dn))

    # 2. For each central space, get outer multiplicity of FTrees (triangles)
    OMcount = 0; om_accumul = Int[1]
    FTree_oms = Tuple{Int, Int}[]

    spaces_dict = Vector{NTuple{NZ, Int}}()
    for sp in central_spaces_lst
        om_up = getNsave_omlist(S, up1sp, sp).totalOM
        om_dn = getNsave_omlist(S, dn1sp, sp).totalOM
        OMcount += om_up * om_dn
        push!(om_accumul, OMcount+1)
        push!(FTree_oms, (om_up, om_dn))
        push!(spaces_dict, sp)
    end
    outgoing_larger = is_outgoing_larger(up1sp, dn1sp)

    return CGTom{S, NZ}(OMcount, spaces_dict, om_accumul, FTree_oms, outgoing_larger)
end

function get_resulting_spaces(up1sp::NTuple{U1, NTuple{NZ, Int}},
    dn1sp::NTuple{D1, NTuple{NZ, Int}},
    up2sp::NTuple{U2, NTuple{NZ, Int}},
    dn2sp::NTuple{D2, NTuple{NZ, Int}},
    legs1::NTuple{M, Int},
    legs2::NTuple{M, Int}) where {U1, D1, U2, D2, NZ, M}

    up3sp = NTuple{NZ, Int}[]
    dn3sp = NTuple{NZ, Int}[]

    for i in 1:U1 if !(i in legs1) push!(up3sp, up1sp[i]) end end
    for i in 1:D1 if !(U1+i in legs1) push!(dn3sp, dn1sp[i]) end end
    for i in 1:U2 if !(i in legs2) push!(up3sp, up2sp[i]) end end
    for i in 1:D2 if !(U2+i in legs2) push!(dn3sp, dn2sp[i]) end end

    return Tuple(sort(up3sp)), Tuple(sort(dn3sp))
end


function Base.show(io::IO, cont::AbstractCG3contract)
    print(io, "$(typeof(cont))\n")
    print(io, "ins=$(cont.ins), outs=$(cont.outs)\n")
    for k in keys(cont.coeff)
        print(io, "coeff[$k]=$(cont.coeff[k])")
        den_nfac = cont.coeff_nfac[k].den
        print(io, den_nfac == 1 ? "\n" : "/$(den_nfac)\n")
    end
end

check_irange(::Type{<:FourCG3s{S}}, i::Int) where {S<:NonabelianSymm} =
    @assert i >= 1 && i <= 4

get_dict_param(::Type{<:FourCG3s{S}}) where {S<:NonabelianSymm} = (4, 4)

create_CG3cont(cont::FourCG3s{S},
    new_ins::NTuple{2, NTuple{NZ, Int}},
    new_outs::NTuple{2, NTuple{NZ, Int}},
    new_coeffs::Dict{NTuple{4, NTuple{NZ, Int}}, Array{BigInt, 4}},
    new_coeff_nfac::Dict{NTuple{4, NTuple{NZ, Int}}, Rational{BigInt}},
    sortcheck::Bool;
    verbose=0) where {S<:NonabelianSymm, NZ} =

    FourCG3s{S, NZ}(new_ins, new_outs, new_coeffs, new_coeff_nfac)




function create_s2res(::Type{S}, 
    ins::NTuple{NI, NTuple{NZ, Int}}, 
    outs::NTuple{NO, NTuple{NZ, Int}}, 
    coeff::Dict{NTuple{M, NTuple{NZ, Int}}, Array{BigInt, D}}, 
    coeff_nfac::Dict{NTuple{M, NTuple{NZ, Int}}, Rational{BigInt}}) where {S<:NonabelianSymm, NI, NO, NZ, M, D}

    @assert M == max(0, D-1)
    for k in keys(coeff)
        if iszero(coeff[k])
            delete!(coeff, k)
            delete!(coeff_nfac, k)
        end
    end
    return Step2_result{S, NI, NO, NZ, M, D}(ins, outs, coeff, coeff_nfac)
end

function create_CG3cont(cont::Step2_result{S, NI, NO, NZ, M, D},
    new_ins::NTuple{NI, NTuple{NZ, Int}},
    new_outs::NTuple{NO, NTuple{NZ, Int}},
    new_coeffs::Dict{NTuple{M, NTuple{NZ, Int}}, Array{BigInt, D}},
    new_coeff_nfac::Dict{NTuple{M, NTuple{NZ, Int}}, Rational{BigInt}},
    sortcheck::Bool;
    verbose=0) where {S<:NonabelianSymm, NI, NO, NZ, M, D}

    @assert M==max(D-1,0)
    create_s2res(S, new_ins, new_outs, new_coeffs, new_coeff_nfac)
end

function rsymiofunc_prestep2(s1res::Step1_result{S, NI, NO, NZ},
    intermsp::NTuple{M1, NTuple{NZ, Int}},
    i::Int) where {S<:NonabelianSymm, NI, NO, NZ, M1}

    @assert NI > 0 && NO > 0
    @assert i == 1 || i == NO+1
    in2 = intermsp[NO]
    if i == 1
        in1 = NO > 1 ? intermsp[1] : s1res.outs[1]
        out = s1res.upper_contracted
    else # i == NO+1
        in1 = NI > 1 ? intermsp[NO+1] : s1res.ins[1]
        out = s1res.lower_contracted
    end
    return in1, in2, out
end

# Contract two CGTs and make it into canonical form
# Key function to get X-symbol
# CGT1 consists of CGT1up and CGT1down, similarly for CGT2
# CGT1legs[i]th leg of CGT1 is contracted with CGT2legs[i]th leg of CGT2
# Index of CGT: [FTree_up.ins..., FTree_down.ins...]
# FTree_up, FTree_down are CGT1up, CGT1down for CGT1, respectively
function contN2canonical(
    CGT1up::FTree{S},
    CGT1down::FTree{S},
    CGT2up::FTree{S},
    CGT2down::FTree{S},
    CGT1legs::NTuple{M, Int},
    CGT2legs::NTuple{M, Int};
    verbose=0) where {S<:NonabelianSymm, M}

    # The number of added zero spaces in input and output
    # This amount of zero spaces are eventually removed
    nz_in, nz_out = 0, 0
    
    n1up, n1down = length(CGT1up.ins), length(CGT1down.ins)
    n2up, n2down = length(CGT2up.ins), length(CGT2down.ins)
    @assert M <= min(n1up+n1down, n2up+n2down)
    # which legs of CGT1up & CGT1down & ... are contracted
    l1up, l1dn, l2up, l2dn = separate_legs(CGT1legs, CGT2legs, n1up, n1down, n2up, n2down)
    # Step 1
    step1_up, nzin1, nzout1 = do_step1(CGT2down, CGT1up, l2dn, l1up; verbose)
    step1_dn, nzin2, nzout2 = do_step1(CGT1down, CGT2up, l1dn, l2up; verbose)
    if verbose > 1
        println("step1_up:")
        println(step1_up)

        println("step1_dn:")
        println(step1_dn)
    end
    nz_in += nzin1 + nzin2; nz_out += nzout1 + nzout2
    if verbose > 1
        println("$(nz_in) dummy 0s in input, $(nz_out) dummy 0s in output")
    end

    # Before step 2, apply R-symbol appropriately to the step1_up and step1_dn
    # Apply R-symbol to the top side CG3 of the step1_up
    apply_Rsymbol!(step1_up, 1, rsymiofunc_prestep2; verbose)
    # Apply R-symbol to the bottom side CG3 of the step1_dn
    apply_Rsymbol!(step1_dn, length(step1_dn.outs)+1, rsymiofunc_prestep2; verbose)

    # Possible plan to implement the step 2:
    # 1. generalize the apply_Rsymbol and apply_Fsymbol functions
    # 2. Implement the function which transforms 4 CG3s contraction to 2 CG3s.
    # 3. Define the struct which will be used in step 3 and 
    # create it from step1_up and step1_dn
    step2_res, nr_input, nr_output = do_step2(step1_up, step1_dn; verbose)

    # Step 3
    # This is just an backward F-symbol application and permutation of triangles.
    # If step 2 is done, implementing it would be easy.
    step3_res = do_step3(step2_res, nr_input, nr_output, nz_in, nz_out; verbose)
    return step3_res
end

# Do the second step of the X-symbol calculation
function do_step2(s1up::Step1_result{S, NIU, NOU, NZ}, 
    s1dn::Step1_result{S, NID, NOD, NZ}; 
    verbose=0) where {S<:NonabelianSymm, NIU, NOU, NID, NOD, NZ}

    @assert s1up.upper_contracted == s1dn.lower_contracted
    @assert s1up.lower_contracted == s1dn.upper_contracted
    int1, int3 = s1dn.upper_contracted, s1dn.lower_contracted

    s1upkeys = sort(collect(keys(s1up.coeff)); by=x->getspaces(s1up, x))
    s1dnkeys = sort(collect(keys(s1dn.coeff)); by=x->getspaces(s1dn, x))
    nkeys1, nkeys2 = length(s1upkeys), length(s1dnkeys)
    NKL = NIU + NOU + NID + NOD - 3 # New key length
    NDL = NKL + 1 # The dimension of resulting array

    new_outs = (s1dn.outs..., s1up.outs...)
    new_ins = (s1up.ins..., s1dn.ins...)
    @assert length(new_outs) == NOU + NOD
    @assert length(new_ins) == NIU + NID

    new_coeff = Dict{NTuple{NKL, NTuple{NZ, Int}}, Array{BigInt, NDL}}()
    new_coeff_nfac = Dict{NTuple{NKL, NTuple{NZ, Int}}, Rational{BigInt}}()

    i1, i2 = 1, 1
    while i1 <= nkeys1
        k1 = s1upkeys[i1]
        out2, in1, int4 = getspaces(s1up, k1)
        deg1 = 0; i1_ = i1 # Degeneracy of input/output/internal spaces
        while i1_ <= nkeys1 && (out2, in1, int4) == getspaces(s1up, s1upkeys[i1_])
            deg1 += 1; i1_ += 1
        end
        
        while i2 <= nkeys2
            k2 = s1dnkeys[i2]
            out1, in2, int2 = getspaces(s1dn, k2)
            deg2 = 0; i2_ = i2
            while i2_ <= nkeys2 && (out1, in2, int2) == getspaces(s1dn, s1dnkeys[i2_])
                deg2 += 1; i2_ += 1
            end

            # Get outer mutliplicity of each CG3, and save them in the tuple sz
            coeff_up, coeff_dn = s1up.coeff[k1], s1dn.coeff[k2]
            sz4, sz1 = size(coeff_up, 1), size(coeff_up, NOU+1)
            sz2, sz3 = size(coeff_dn, 1), size(coeff_dn, NOD+1)
            sz = (sz1, sz2, sz3, sz4)

            FourCG3s_key = (out1, out2, in1, in2, int1, int2, int3, int4)
            # If the data is not cached, calculate it
            FourCG3s_res, nfacs = compute_FourCG3s(S, FourCG3s_key, sz; verbose)

            # Apply it and store the result in Step2_result struct
            for n1 in 0:deg1-1
                for n2 in 0:deg2-1
                    i1_ = i1 + n1; i2_ = i2 + n2
                    key1, key2 = s1upkeys[i1_], s1dnkeys[i2_]
                    arr1, arr2 = s1up.coeff[key1], s1dn.coeff[key2]
                    nfac1::Rational{BigInt} = s1up.coeff_nfac[key1]
                    nfac2::Rational{BigInt} = s1dn.coeff_nfac[key2]

                    for q in keys(FourCG3s_res)
                        fourCG3s_arr::Array{BigInt, 6} = FourCG3s_res[q]
                        fourCG3s_nfac::Rational{BigInt} = nfacs[q]

                        nk = (key2[NOD+1:end]...,
                            key1[NOU+1:end]...,
                            q,
                            key1[1:NOU-1]...,
                            key2[1:NOD-1]...) # new key
                        @assert length(nk) == NKL

                        # Multiply integer arrays
                        arr_multiplied = final_contract_step2(arr1, arr2, fourCG3s_arr,
                        Val(NIU), Val(NOU), Val(NID), Val(NOD))
                        nfac_multiplied = nfac1 * nfac2 * fourCG3s_nfac

                        arr_rp, nfac_rp = arr_relprime(arr_multiplied, nfac_multiplied)
                        if !haskey(new_coeff, nk)
                            new_coeff[nk] = arr_rp
                            new_coeff_nfac[nk] = nfac_rp
                        else
                            new_coeff[nk], new_coeff_nfac[nk] =
                            add_arrs(new_coeff[nk], arr_rp, new_coeff_nfac[nk], nfac_rp)
                        end
                    end
                end
            end
            i2 += deg2
        end
        i1 += deg1; i2 = 1
    end

    step2_result = Step2_result{S, NIU+NID, NOU+NOD, NZ, NKL, NDL}(new_ins, new_outs, new_coeff, new_coeff_nfac)
    return step2_result, NID, NOU
end

# upper step1_result * lower step1_result * fourCG3s result
@generated function final_contract_step2(
    arr1::Array{BigInt, MU}, 
    arr2::Array{BigInt, MD}, 
    fourCG3s_arr::Array{BigInt, 6},
    ::Val{NIU}, ::Val{NOU}, ::Val{NID}, ::Val{NOD}) where {MU, MD, NIU, NOU, NID, NOD}

    @assert MU == NIU + NOU && MD == NID + NOD
    ius = [Symbol(:iu, i) for i in 1:NIU-1]
    ous = [Symbol(:ou, i) for i in 1:NOU-1]
    ids = [Symbol(:id, i) for i in 1:NID-1]
    ods = [Symbol(:od, i) for i in 1:NOD-1]

    result_inds = (:cu, ids..., ius..., :cd, ous..., ods...)
    arr1_inds = (:b4, ous..., :b1, ius...)
    arr2_inds = (:b2, ods..., :b3, ids...)
    fourCG3s_inds = (:b1, :b2, :b3, :b4, :cu, :cd)


    quote
        @tensor result[$(result_inds...)] := arr1[$(arr1_inds...)] * 
                                            arr2[$(arr2_inds...)] * 
                                            fourCG3s_arr[$(fourCG3s_inds...)]
        return result
    end
end

# Get in2, out1, int4 spaces if applied for upper Step1_result
# Get in1, out2, int2 spaces if applied for lower Step1_result
function getspaces(s1res::Step1_result{S, NI, NO, NZ},
    key::Tuple) where {S<:NonabelianSymm, NI, NO, NZ}
    out = NO > 1 ? key[1] : s1res.outs[1]
    in = NI > 1 ? key[NO+1] : s1res.ins[1]
    int = key[NO]
    return out, in, int
end

# Central function to reduce 4 CG3s to 2 CG3s
function compute_FourCG3s(::Type{S},
    spaces::NTuple{8, NTuple{NZ, Int}},
    sz::NTuple{4, Int};
    verbose=0) where {S<:NonabelianSymm, NZ}

    res_dict = Dict{NTuple{NZ, Int}, Array{BigInt, 6}}()
    res_nfac = Dict{NTuple{NZ, Int}, Rational{BigInt}}()

    out1, out2, in1, in2, intsps... = spaces
    if verbose > 1
        println()
        println("input spaces: $in1, $in2")
        println("output spaces: $out1, $out2")
        println("intermediate spaces: $intsps")
        println("outer multiplicities: $sz")
    end
    for I in CartesianIndices(sz)
        arr = zeros(BigInt, sz...)
        arr[I] = 1
        fourCG3s = FourCG3s{S, NZ}((in1, in2), (out1, out2), Dict(intsps => arr), Dict(intsps => 1//1))
        if verbose > 1 println(fourCG3s) end
        res_arr, res_fac = reduce_FourCG3s(fourCG3s)
        if verbose > 1
            println("Modifying $I-th coordinate")
            for k in keys(res_arr)
                println("res_arr[$k]=$(res_arr[k])")
                println("res_fac[$k]=$(res_fac[k])")
            end
        end

        # Save the results in res_dict and res_nfac
        for (k, arr2d) in res_arr
            if !haskey(res_dict, k)
                @assert !haskey(res_nfac, k)
                arr6d_sz = (sz..., size(arr2d)...)
                arr6d = zeros(BigInt, arr6d_sz)
                arr6d[I, :, :] = arr2d
                res_dict[k] = arr6d
                res_nfac[k] = res_fac[k]
            else
                orig_nfac::BigInt = res_nfac[k].den
                new_nfac::BigInt = res_fac[k].den

                final_nfac::BigInt = lcm(orig_nfac, new_nfac)
                res_nfac[k] = 1//final_nfac

                res_dict[k] *= div(final_nfac, orig_nfac)
                res_dict[k][I, :, :] = arr2d * div(final_nfac, new_nfac)
            end
        end
    end
    if verbose > 1 
        println("resulting dictionary: $res_dict")
        println("resulting dividing factors: $res_nfac")
    end
    return res_dict, res_nfac
end

function reduce_FourCG3s(
    fourCG3s::FourCG3s{S, NZ}) where {S<:NonabelianSymm, NZ}

    # First, apply the anomalous F-symbol to the second intermediate space
    fourCG3s = apply_Fsymbol_strange(fourCG3s)
    
    # Next, switch two legs (int4, out2) by applying 2 R-symbols and 1 F-symbol
    apply_Rsymbol!(fourCG3s, 2, (obj, ints, i) -> (ints[1], obj.ins[2], ints[2]))
    ft_switch = (obj, rem, i) -> (obj.ins[2], obj.ins[1], rem[3], rem[1])
    fourCG3s = apply_Fsymbol(fourCG3s, fourCG3s.ins, fourCG3s.outs, 1, true, ft_switch; 
        verbose=0, permute_updown=true)
    apply_Rsymbol!(fourCG3s, 1, (obj, ints, i) -> (obj.ins[1], obj.ins[2], ints[1]))

    # Finally, apply F/R-symbols to the third and fourth CG3s.
    apply_Rsymbol!(fourCG3s, 4, (obj, ints, i) -> (ints[4], obj.outs[2], ints[3]))
    ft_down = (obj, rem, i) -> (obj.outs[1], obj.outs[2], rem[3], rem[2])
    fourCG3s = apply_Fsymbol(fourCG3s, fourCG3s.ins, fourCG3s.outs, 3, true, ft_down; 
        verbose=0, permute_updown=false)

    res_arr, res_nfac = reduce2_twoCG3s(fourCG3s)
    return res_arr, res_nfac
end

function reduce2_twoCG3s(
    fourCG3s::FourCG3s{S, NZ}) where {S<:NonabelianSymm, NZ}

    # This function reduces the 4 CG3s to 2 CG3s
    # The result is stored in res_arr and res_nfac
    # Key is just a unique intermediate space
    res_arr = Dict{NTuple{NZ, Int}, Array{BigInt, 2}}()
    res_nfac = Dict{NTuple{NZ, Int}, Rational{BigInt}}()

    repdim_dict = Dict{NTuple{NZ, Int}, Int}()
    ks = sort(collect(keys(fourCG3s.coeff)))
    for k in ks
        k1, k2, k3, k4 = k
        if k1 != k3 continue end
        d1 = get_dim!(S, repdim_dict, k1)
        d2 = get_dim!(S, repdim_dict, k2)

        arr::Array{BigInt, 4} = fourCG3s.coeff[k]
        nfac::Rational{BigInt} = fourCG3s.coeff_nfac[k]

        _, cg3fac = load_cg3blk(S, BigInt, (k1, k4), [k2])[k2]
        dens = [x.den for x in cg3fac]
        fac_rp, n_rp::Rational{BigInt} = arr_relprime(dens * d2, 1//d1)

        fac_mat = Diagonal(fac_rp)
        @tensor arr_new[o1, o4] := arr[o1, o2, o3, o4] * fac_mat[o2, o3]
        arr_after, fac_after = arr_relprime(arr_new, nfac * n_rp)

        if !haskey(res_arr, k1)
            res_arr[k1] = arr_after
            res_nfac[k1] = fac_after
        else
            # If the key already exists, add the new array
            res_arr[k1], res_nfac[k1] = add_arrs(res_arr[k1], arr_after, res_nfac[k1], fac_after)
        end
    end
    return res_arr, res_nfac
end

function get_dim!(::Type{S}, 
    repdim_dict::Dict{NTuple{NZ, Int}, Int}, 
    key::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ}

    if !haskey(repdim_dict, key) repdim_dict[key] = dimension(S, key) end
    return repdim_dict[key]
end

function apply_Fsymbol_strange(
    fourCG3s::FourCG3s{S, NZ}) where {S<:NonabelianSymm, NZ}

    new_coeff = Dict{NTuple{4, NTuple{NZ, Int}}, Array{BigInt, 4}}()
    new_coeff_nfac = Dict{NTuple{4, NTuple{NZ, Int}}, Rational{BigInt}}()

    in1, in3 = fourCG3s.outs[1], fourCG3s.ins[2]
    # In my algorithm, it has only one entry.
    for (key, coeff) in fourCG3s.coeff
        in2, e, f = key[2], key[1], key[3]
        nfac = fourCG3s.coeff_nfac[key]
        
        # Find the possible list of possible 'q',
        # which can be obtained both from e ⊗ in3 and in1 ⊗ f
        vo1 = getNsave_validout(S, Tuple(sort([e, in3])))
        vo2 = getNsave_validout(S, Tuple(sort([in1, f])))
        qlst = collect(Set(vo1.out_spaces) ∩ Set(vo2.out_spaces))

        for q in qlst
            fsym = getNsave_Fsymbol(S, BigInt, in1, in2, in3, q)
            @assert !isnothing(fsym)
            # arr: Has indices (μ, κ, ν, λ)
            arr, fsym_nfac::Rational{BigInt} = getpartof_step2(fsym, e, f)

            new_key = (key[1], q, key[3], key[4])
            @tensor new_arr[a, ν, λ, b] := coeff[a, μ, κ, b] * arr[μ, κ, ν, λ]
            new_nfac = nfac * fsym_nfac
            new_arr, new_nfac = arr_relprime(new_arr, new_nfac)
            
            @assert !haskey(new_coeff, new_key)
            @assert !haskey(new_coeff_nfac, new_key)
            new_coeff[new_key] = new_arr
            new_coeff_nfac[new_key] = new_nfac
        end
    end
    return FourCG3s{S, NZ}(fourCG3s.ins, fourCG3s.outs, new_coeff, new_coeff_nfac)
end

function getpartof_step2(
    fsym::Fsymbol{S}, 
    e::NTuple{NZ, Int}, 
    f::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ}

    erange, ν, μ = get_range_fsym(fsym.es_list, e)
    frange, λ, κ = get_range_fsym(fsym.fs_list, f)
    fsym_4darr = reshape(fsym.fsym_mat[erange, frange], (ν, μ, λ, κ))

    permuted = permutedims(fsym_4darr, (2, 4, 1, 3)) # To (μ, κ, ν, λ)
    _, νfac = load_cg3blk(S, BigInt, (e, fsym.in3), [fsym.out])[fsym.out]
    _, λfac = load_cg3blk(S, BigInt, (fsym.in1, f), [fsym.out])[fsym.out]
    νfac_mat, λfac_mat = Diagonal(νfac), Diagonal(λfac)

    @tensor result[μ, κ, ν, λ] := permuted[μ, κ, νp, λp] * 
        νfac_mat[νp, ν] * λfac_mat[λp, λ]

    @assert result isa Array{Rational{BigInt}, 4}
    nfac::Integer = lcm(denominator.(result))
    return Array{BigInt, 4}(result * nfac), 1//nfac
end

function get_range_fsym(
    intsp_list::Vector{Tuple{NTuple{NZ, Int}, Tuple{Int, Int}}}, 
    intsp::NTuple{NZ, Int}) where {NZ}

    i = 1
    for (sp, (om1, om2)) in intsp_list
        if sp == intsp return i:i+om1*om2-1, om1, om2 end
        i += om1 * om2
    end
end

function separate_legs(
    CGT1legs::NTuple{M, Int},
    CGT2legs::NTuple{M, Int},
    n1up::Int,
    n1down::Int,
    n2up::Int,
    n2down::Int) where M

    l1up, l1dn, l2up, l2dn = Int[], Int[], Int[], Int[]
    for i in 1:M
        l1, l2 = CGT1legs[i], CGT2legs[i]
        # The index should be in [1, n1up + n1down (=number of legs of CGT1)] 
        @assert 1 <= l1 && l1 <= n1up + n1down
        @assert 1 <= l2 && l2 <= n2up + n2down
        up1, up2 = l1 <= n1up, l2 <= n2up
        # outgoing leg of one CGT is contracted with incoming leg of the other CGT
        @assert up1 != up2
        if up1 push!(l1up, l1) else push!(l1dn, l1 - n1up) end
        if up2 push!(l2up, l2) else push!(l2dn, l2 - n2up) end
    end
    return Tuple(l1up), Tuple(l1dn), Tuple(l2up), Tuple(l2dn)
end



# The first step to get X-symbol
# FTree1 is on the upper position, FTree2 is on the lower position
function do_step1(FTree1::FTree{S, N1}, 
    FTree2::FTree{S, N2}, 
    ci1::NTuple{M, Int}, 
    ci2::NTuple{M, Int};
    verbose=0) where {S<:NonabelianSymm, N1, N2, M}

    @assert M <= N1 && M <= N2
    # Count the number of added zero spaces in input and output
    nz_in, nz_out = 0, 0
    # Get the permutation to move the contraction indices to the left
    perm1 = get_perm_step1(N1, ci1)
    perm2 = get_perm_step1(N2, ci2)

    # Permute the FTree objects
    FTree1_perm = permute!(FTree1, perm1)
    FTree2_perm = permute!(FTree2, perm2)

    # For some exceptional cases, we should introduce dummy 0s
    # If there are no contracted legs
    # This is done by add_singdim function
    if M == 0
        FTree1_perm = add_singdim(FTree1_perm, false; verbose=0)
        FTree2_perm = add_singdim(FTree2_perm, false; verbose=0)
    end

    if N1 == M
        FTree1_perm = add_singdim(FTree1_perm, true; verbose=0)
        nz_out += 1
    end

    if N2 == M
        FTree2_perm = add_singdim(FTree2_perm, true; verbose=0)
        nz_in += 1
    end

    for i=M-1:-1:1
        # If every leg of FTree1_perm is contracted, no need to apply Fsymbol
        # Similarly for FTree2_perm
        ft1 = (a...) -> fsymiofunc_step1(Val(M), a...)
        FTree1_perm = apply_Fsymbol(FTree1_perm, FTree1_perm.ins, 
        FTree1_perm.outs, i, false, ft1)

        ft2 = (a...) -> fsymiofunc_step1(Val(M), a...)
        FTree2_perm = apply_Fsymbol(FTree2_perm, FTree2_perm.ins, 
        FTree2_perm.outs, i, false, ft2) 
    end

    # Get a step1_result object by combining two FTree objects
    s1_result = combine_2FTrees(FTree1_perm, FTree2_perm, Val(max(M, 1)); verbose=verbose)
    if verbose > 1 println("step1 result is ") end
    if verbose > 1 println(s1_result) end
    return s1_result, nz_in, nz_out
end

function combine_2FTrees(FTree1::FTree{S, N1}, 
    FTree2::FTree{S, N2}, 
    ::Val{M}; verbose=0) where {S<:NonabelianSymm, N1, N2, M}
    # Since dummy 0 is added, we can exclude the case N1 == M or N2 == M
    NZ = nzops(S)
    @assert N1 > M && N2 > M
    if verbose > 1
        println("Combining two FTree objects:")
        println("FTree1:")
        println(FTree1)
        println("FTree2:")
        println(FTree2)
        println("M = $M")
    end
    # The number of intermediate legs for the resulting obj.
    # The meaning of Mn is 'M new', similarly for other variables
    Mn = N1 + N2 - 2*M - 1 
    NO, NI = N1 - M, N2 - M

    # Contract the arrays in FTree1 and FTree2. 
    # This is similar to the contractblks_mw function in Fsymbol.jl file
    coeffn = Dict{NTuple{Mn, NTuple{NZ, Int}}, Array{BigInt}}()
    coeff_nfacn = Dict{NTuple{Mn, NTuple{NZ, Int}}, Rational{BigInt}}()

    blks1keys = sort(collect(keys(FTree1.coeff)); by=x->x[1:M-1])
    blks2keys = sort(collect(keys(FTree2.coeff)); by=x->x[1:M-1])
    nkeys1, nkeys2 = length(blks1keys), length(blks2keys)

    i1, i2 = 1, 1
    while i1 <= nkeys1 && i2 <= nkeys2
        k1, k2 = blks1keys[i1], blks2keys[i2]
        if k1[1:M-1] == k2[1:M-1]
            # If the keys match, contract the arrays
            # Find the degeneracy of the keys deg1, deg2
            deg1, deg2 = 0, 0
            i1_, i2_ = i1, i2
            while i1_ <= nkeys1 && blks1keys[i1_][1:M-1] == k1[1:M-1]
                deg1 += 1; i1_ += 1
            end 
            while i2_ <= nkeys2 && blks2keys[i2_][1:M-1] == k2[1:M-1]
                deg2 += 1; i2_ += 1
            end
            
            j1, j2 = i1, i2
            for _=1:deg1
                for _=1:deg2
                    # Contract the j1th array of FTree1 and j2th array of FTree2
                    # and store the result in coeffn, coeff_nfacn
                    key1, key2 = blks1keys[j1], blks2keys[j2]
                    coeff1, coeff2 = FTree1.coeff[key1], FTree2.coeff[key2]
                    dfac1::Rational{BigInt}, dfac2::Rational{BigInt} = 
                        FTree1.coeff_nfac[key1], FTree2.coeff_nfac[key2]

                    # Generate the new key from key1 and key2
                    nkey = (key1[M:end]..., M>1 ? key1[1] : FTree1.ins[end], key2[M:end]...)

                    # Generate the new array by contracting coeff1 and coeff2
                    # with appropriate normalization factor of CG3s
                    coeff_arr = cont_coeffs_s1(S, key1, FTree1.ins, coeff1, coeff2, Val(M); verbose)
                    dfacn::Rational{BigInt} = dfac1 * dfac2

                    coeff_arr, dfacn = arr_relprime(coeff_arr, dfacn)
                    # Save it into coeffn and coeff_nfacn
                    # If the array in coeffn already exists, add the new array
                    # If not, create a new entry
                    if haskey(coeffn, nkey)
                        coeff_fin, dfacfin = 
                        add_arrs(coeffn[nkey], coeff_arr, coeff_nfacn[nkey], dfacn)
                        coeffn[nkey], coeff_nfacn[nkey] = coeff_fin, dfacfin
                    else
                        coeffn[nkey] = coeff_arr
                        coeff_nfacn[nkey] = dfacn
                    end
                    j2 += 1
                end
                j2 = i2; j1 += 1
            end

            i1, i2 = i1_, i2_
        elseif k1[1:M-1] < k2[1:M-1]
            i1 += 1
        else
            i2 += 1
        end
    end
    for (k, arr) in coeffn
        if norm(arr) == 0
            delete!(coeffn, k)
            delete!(coeff_nfacn, k)
        end
    end
    return Step1_result{S, NI, NO, NZ, Mn}(
        FTree2.ins[1:NI], FTree1.ins[1:NO],
        FTree1.outs[1], FTree2.outs[1], 
        coeffn, coeff_nfacn
        )
end

function cont_coeffs_s1(
    ::Type{S},
    key1::Tuple{Vararg{NTuple{NZ, Int}}},
    ins1::Tuple{Vararg{NTuple{NZ, Int}}},
    coeff1::Array{BigInt, D1}, 
    coeff2::Array{BigInt, D2}, 
    ::Val{M}; verbose=0) where {S<:NonabelianSymm, D1, D2, M, NZ}

    @assert D1 >= M && D2 >= M
    @assert D1 == length(ins1) - 1
    # Need to contract with normalization factor of CG3s
    nfac_dens = Vector{Vector{BigInt}}()
    for i in 1:M-1
        in2, out = ins1[end-i+1], key1[i]
        in1 = i==M-1 ? ins1[end-M+1] : key1[i+1]
        _, nfac = load_cg3blk(S, BigInt, (in1, in2), [out])[out]
        push!(nfac_dens, [x.den for x in nfac])
    end
    nfac_mats = [Diagonal(dens) for dens in nfac_dens]
    return cont_coeffs_s1_int(coeff1, coeff2, nfac_mats, Val(M))
end

@generated function cont_coeffs_s1_int(
    coeff1::Array{BigInt, D1}, 
    coeff2::Array{BigInt, D2}, 
    nfac_mats::Vector{Diagonal{BigInt, Vector{BigInt}}},
    ::Val{M}) where {D1, D2, M}

    out_inds1 = [Symbol(:i, i) for i in 1:D1-M+1]
    out_inds2 = [Symbol(:j, j) for j in 1:D2-M+1]
    in1inds = [:i1, (Symbol(:a, i) for i in 2:M)..., (Symbol(:i, i) for i in 2:D1-M+1)...]
    in2inds = [:j1, (Symbol(:b, i) for i in 2:M)..., (Symbol(:j, i) for i in 2:D2-M+1)...]

    conts_facs = [:(nfac_mats[$i][$(in1inds[i+1]), $(in2inds[i+1])])  for i in 1:M-1]
    conts_facs = [:(coeff1[$(in1inds...)]), :(coeff2[$(in2inds...)]), conts_facs...]
    expr = Expr(:call, :*, conts_facs...)
    quote
        @tensor conts_res[$(out_inds1...), $(out_inds2...)] := $(expr)
        return conts_res
    end
end

function do_step3(s2res::Step2_result{S, NI, NO},
    nr_input::Int,
    nr_output::Int,
    nz_in::Int,
    nz_out::Int;
    verbose=0) where {S<:NonabelianSymm, NI, NO}
    @assert NI >= 2 && NO >= 2
    if verbose > 1 println("nr_input=$(nr_input), nr_output=$(nr_output)") end
    if verbose > 1 println(s2res) end

    # First, apply F-symbol in the backward direction to make canonical form
    # In the lower part of Step2_result structure
    in_fsymio = (a...) -> fsymio_in_step3(nr_input, a...)
    out_fsymio = (a...) -> fsymio_out_step3(nr_output, a...)
    for i in 1:nr_input-1
        s2res = apply_Fsymbol(s2res, s2res.ins, s2res.outs, i, true, in_fsymio)
    end

    for i in NI:NI+nr_output-2
        s2res = apply_Fsymbol(s2res, s2res.ins, s2res.outs, i, true, out_fsymio)
    end

    # Remove additional zero spaces
    if nz_in > 0 || nz_out > 0
        s2res = remove_zero_spaces(s2res, nz_in, nz_out; verbose)
    end
    if verbose > 1 
        println("After removing zero spaces:")
        println(s2res) 
    end

    # (Stably) sort the input and output spaces, and permute the object
    final_result = sortio_step3(s2res; verbose)
    if verbose > 1 
        println("The final result after sorting the input and output spaces:")
        println(final_result) 
    end
    return final_result
end

function sortio_step3(s2res::Step2_result{S, NI, NO};
    verbose=0) where {S<:NonabelianSymm, NI, NO}
    @assert NI >= 1 && NO >= 1
    if verbose > 1 
        println("Sorting the input and output spaces") 
        println("Before sorting:")
        println("Input spaces: ", collect(s2res.ins))
        println("Output spaces: ", collect(s2res.outs))
    end
    
    sortperm_ins = sortperm(collect(s2res.ins))
    sortperm_outs = sortperm(collect(s2res.outs))

    if sortperm_ins != collect(1:NI)
        if verbose > 1 println("permute inputs by ", sortperm_ins) end
        s2res = permute!(s2res, Tuple(sortperm_ins), true; verbose)
    end

    if sortperm_outs != collect(1:NO)
        if verbose > 1 println("permute outputs by ", sortperm_outs) end
        s2res = permute!(s2res, Tuple(sortperm_outs), false; verbose)
    end
    return s2res
end

function Base.permute!(s2res::Step2_result{S, NI, NO}, 
    perm::NTuple{N, Int}, 
    is_input::Bool;
    verbose=0) where {S<:NonabelianSymm, NI, NO, N}

    @assert !issorted(perm)
    if is_input @assert length(perm) == NI 
    else @assert length(perm) == NO end
    
    # Decompose the permutation into a list of adjacent-site switch
    switch_list = decompose_perm(perm)
    if verbose > 1 println("Decomposed permutation into adjacent-site switches: ", switch_list) end
    for (i, j) in switch_list
        if verbose > 1 println("Permuting indices $(i), $(j)") end
        s2res = switch_adjacent!(s2res, i, j, is_input; verbose)
    end
    return s2res
end

# If is_input, switch the ith and j(=i+1)th input spaces
# Otherwise, switch the ith and j(=i+1)th output spaces
function switch_adjacent!(s2res::Step2_result{S, NI, NO, NZ, KL}, 
    i1::Int, 
    i2::Int, 
    is_input::Bool; 
    verbose=0) where {S<:NonabelianSymm, NI, NO, NZ, KL}

    @assert i2 == i1 + 1
    if is_input @assert i1 >= 1 && i2 <= NI end
    if !is_input @assert i1 >= 1 && i2 <= NO end

    if is_input
        # Switch the input spaces
        new_ins = (s2res.ins[1:i1-1]...,
            s2res.ins[i2], s2res.ins[i1],
            s2res.ins[i2+1:end]...)
        new_outs = s2res.outs
    else
        # Switch the output spaces
        new_outs = (s2res.outs[1:i1-1]...,
            s2res.outs[i2], s2res.outs[i1],
            s2res.outs[i2+1:end]...)
        new_ins = s2res.ins
    end

    # This is the case when only one R-symbol is needed
    if i1 == 1
        are_same = is_input ? s2res.ins[1] == s2res.ins[2] : 
                              s2res.outs[1] == s2res.outs[2]
        if are_same
            modified_cg3 = is_input ? NI-1 : NI+NO-2
            apply_Rsymbol!(s2res, modified_cg3, rsymiofunc_perm_step3; verbose)
        end
        return create_s2res(S, new_ins, new_outs, 
        s2res.coeff, s2res.coeff_nfac)
    end

    if verbose > 1 println("R_step3") end
    modified_cg3 = is_input ? NI-i1 : NI+NO-1-i1
    apply_Rsymbol!(s2res, modified_cg3, rsymiofunc_perm_step3; verbose)
    if verbose > 1 println("F_step3") end
    s2res_new = apply_Fsymbol(s2res, new_ins, new_outs,
        modified_cg3, true, fsymiofunc_perm_step3; verbose)
    if verbose > 1 println("R_step3") end
    apply_Rsymbol!(s2res_new, modified_cg3+1, rsymiofunc_perm_step3; verbose)
    return s2res_new
end

function rsymiofunc_perm_step3(s2res::Step2_result{S, NI, NO, NZ, KL},
    intermsp::NTuple{KL, NTuple{NZ, Int}},
    i::Int) where {S<:NonabelianSymm, NI, NO, NZ, KL}

    @assert KL == NI + NO - 3; @assert 1 <= i && i <= KL+1
    if i < NI # If ith cg3 is associated with input
        in1 = i < NI-1 ? intermsp[i] : s2res.ins[1]
        in2 = s2res.ins[NI-i+1]
        if i == 1 rout = NO==1 ? s2res.outs[1] : intermsp[NI-1]
        else rout = intermsp[i-1] end

    else # If ith cg3 is associated with output
        in1 = i < NI+NO-2 ? intermsp[i] : s2res.outs[1]
        in2 = s2res.outs[NI+NO-i]
        if i == NI rout = NI==1 ? s2res.ins[1] : intermsp[NI-1]
        else rout = intermsp[i-1] end
    end
    return in1, in2, rout
end

function fsymiofunc_perm_step3(s2res::Step2_result{S, NI, NO, NZ, M},
    rem::NTuple{KL, NTuple{NZ, Int}},
    i::Int) where {S<:NonabelianSymm, NI, NO, NZ, M, KL}

    @assert KL == NI + NO - 4
    if i < NI-1 # If ith cg3 is associated with input
        in2 = i==NI-2 ? s2res.ins[1] : rem[i]
        in1, in3 = s2res.ins[NI-i+1], s2res.ins[NI-i]
        if i == 1 outsp = NO==1 ? s2res.outs[1] : rem[NI-2]
        else outsp = rem[i-1] end

    elseif i > NI-1 # If ith cg3 is associated with output
        in2 = i < NI+NO-3 ? rem[i] : s2res.outs[1]
        in1, in3 = s2res.outs[NI+NO-i], s2res.outs[NI+NO-i-1]
        if i == NI outsp = NI==1 ? s2res.ins[1] : rem[NI-1]
        else outsp = rem[i-1] end
        
    else error("The index i=$(i) is not valid for Step3_result with $(NI) inputs and $(NO) outputs") end
    return in1, in2, in3, outsp
end

check_irange(::Type{<:Step2_result{S, NI, NO}}, i::Int) where {S<:NonabelianSymm, NI, NO} = 
    @assert 1 <= i && i <= NI+NO-3 && i != NI-1

get_dict_param(::Type{<:Step2_result{S, NI, NO}}) where {S<:NonabelianSymm, NI, NO} = 
NI+NO-3, NI+NO-2


function remove_zero_spaces(s2res::Step2_result{S, NI, NO, NZ}, 
    nz_in::Int, 
    nz_out::Int;
    verbose=0) where {S<:NonabelianSymm, NI, NO, NZ}
    @assert NI >= 2 && NO >= 2
    if verbose > 1 println("Removing zero spaces: nz_in=$(nz_in), nz_out=$(nz_out)") end
    if verbose > 1 println("Before removing zero spaces:") end
    if verbose > 1 println(s2res) end

    zero_qlabel = Tuple(0 for _=1:NZ)
    # Find where are zero spaces in inpute and outputs
    idx_in = findall(x->x==zero_qlabel, s2res.ins)
    idx_out = findall(x->x==zero_qlabel, s2res.outs)
    @assert length(idx_in) >= nz_in && length(idx_out) >= nz_out

    removed_in = idx_in[end-nz_in+1:end]
    removed_out = idx_out[end-nz_out+1:end]

    # Remaining indices
    remain_in = sort(collect(setdiff(Set(1:NI), Set(removed_in))))
    remain_out = sort(collect(setdiff(Set(1:NO), Set(removed_out))))

    new_in = Tuple(s2res.ins[i] for i in remain_in)
    new_out = Tuple(s2res.outs[i] for i in remain_out)

    if verbose > 1
        println("Remaining input indices: ", remain_in)
        println("Remaining output indices: ", remain_out)
        println("Remaining input spaces: ", new_in)
        println("Remaining output spaces: ", new_out)
    end

    # If every spaces are eliminated, add a single zero space
    if isempty(new_in) new_in = (zero_qlabel,) end
    if isempty(new_out) new_out = (zero_qlabel,) end
    # New number of inputs, new number of outputs
    NIn, NOn = length(new_in), length(new_out)
    @assert NIn >= 1 && NOn >= 1
    ND = NIn + NOn - 2
    NKL = max(0, ND - 1)

    # Which CG3 survives
    survive_cg3s_in = get_survive_cg3s(NI, remain_in) 
    survive_cg3s_out = get_survive_cg3s(NO, remain_out) .+ (NI - 1)

    survive_intsp_in = survive_cg3s_in[1:end-1] .- 1
    survive_intsp_out = survive_cg3s_out[1:end-1] .- 1

    if verbose > 1
        println("Survived cg3s associated to input: ", survive_cg3s_in)
        println("Survived cg3s associated to output: ", survive_cg3s_out)
        println("Survived intermediate spaces associated to input: ", survive_intsp_in)
        println("Survived intermediate spaces associated to output: ", survive_intsp_out)
    end

    remain_cg3s = sort(vcat(survive_cg3s_out, survive_cg3s_in))
    # If input or output consists of only one space, '
    # there is no center intermediate space
    center_intsp = NIn > 1 && NOn > 1 ? [NI-1] : Int[]
    remain_intsps = sort(vcat(survive_intsp_in, center_intsp, survive_intsp_out))

    if verbose > 1 println("remaining cg3 indices: ", remain_cg3s) end
    if verbose > 1 println("remaining intermediate space indices: ", remain_intsps) end

    new_coeffs = Dict{NTuple{NKL, NTuple{NZ, Int}}, Array{BigInt, ND}}()
    new_coeff_nfac = Dict{NTuple{NKL, NTuple{NZ, Int}}, Rational{BigInt}}()
    if verbose > 1 println("Length of key: $NKL") end
    if verbose > 1 println("Dimension of coeff arrays: $ND") end

    for k in keys(s2res.coeff)
        coeff, nfac = s2res.coeff[k], s2res.coeff_nfac[k]
        if verbose > 1 println("Original key: ", k) end
        if verbose > 1 println("Original coeff: ", coeff) end
        nk = length(remain_intsps)>0 ? k[remain_intsps] : ()
        coeffsz_new = size(coeff)[remain_cg3s]
        if length(remain_cg3s)==0 new_coeff = fill(BigInt(coeff[]))
        else new_coeff = reshape(coeff, coeffsz_new...) end
        new_coeffs[nk] = new_coeff
        new_coeff_nfac[nk] = nfac
    end

    s2res_removed = Step2_result{S, NIn, NOn, NZ, NKL, ND}(
        new_in, new_out, new_coeffs, new_coeff_nfac)
    return s2res_removed
end

function get_survive_cg3s(N::Int,
    remaining::Vector{Int})

    len = length(remaining)
    if len <= 1 return Int[] end
    return [N+1-i for i in remaining[2:end]]
end


function fsymio_out_step3(nr_output::Int,
    s2res::Step2_result{S, NI, NO}, 
    rem::NTuple{M, NTuple{NZ, Int}}, 
    i::Int) where {S<:NonabelianSymm, NI, NO, M, NZ}

    i = i - NI + 1
    @assert M == NI + NO - 4 && 0 < i && i < nr_output
    outsp = rem[NI-2+i]
    in1 = NO - nr_output > 1 ? rem[NI+nr_output-2] : s2res.outs[1]
    in2 = i==nr_output-1 ? s2res.outs[NO-i] : rem[NI-1+i]
    in3 = s2res.outs[NO-i+1]
    return in1, in2, in3, outsp
end

function fsymio_in_step3(nr_input::Int,
    s2res::Step2_result{S, NI, NO}, 
    rem::NTuple{M, NTuple{NZ, Int}}, 
    i::Int) where {S<:NonabelianSymm, NI, NO, M, NZ}

    @assert M == NI + NO - 4 && i < nr_input
    in1 = NI - nr_input > 1 ? rem[nr_input-1] : s2res.ins[1]
    in2 = i==nr_input-1 ? s2res.ins[NI-i] : rem[i]
    in3 = s2res.ins[NI-i+1]
    outsp = i==1 ? rem[NI-2] : rem[i-1]
    return in1, in2, in3, outsp
end

# Add a singleton dimension at the front (or back) to every array in FTree
# Create a new FTree with an additional CG3 q ⊗ 0 -> q
# Front: The position of new 0 in the incoming spaces
function add_singdim(FTree::FTree{S, N, NZ}, front; verbose=0) where {S<:NonabelianSymm, N, NZ}
    @assert N >= 1
    zero_qlabel = Tuple(0 for _=1:NZ)
    new_ins = front ? (zero_qlabel, FTree.ins...) : (FTree.ins..., zero_qlabel)

    new_coeff = Dict{NTuple{N-1, NTuple{NZ, Int}}, Array{BigInt, N}}()
    new_coeff_nfac = Dict{NTuple{N-1, NTuple{NZ, Int}}, Rational{BigInt}}()
    for (key, coeff) in FTree.coeff
        # Get the new key after dummy zero is added
        if N==1 new_key = ()
        else new_key = front ? (key..., FTree.ins[1]) : (FTree.outs[1], key...) end

        # Get the new array after dummy zero is added
        sz = size(coeff)
        if N==1 new_arr = BigInt[1]
        else new_arr = front ? reshape(coeff, (sz..., 1)) : reshape(coeff, (1, sz...)) end
        new_coeff[new_key] = new_arr
        # Lastly, fill in the normalization factor
        new_coeff_nfac[new_key] = FTree.coeff_nfac[key]
    end
    return create_FTree(S, new_ins, FTree.outs, new_coeff, new_coeff_nfac, false; verbose)
end

function add_arrs(
    arr1::Array{BigInt, D}, 
    arr2::Array{BigInt, D}, 
    nfac1::Rational{BigInt}, 
    nfac2::Rational{BigInt}) where {D}

    @assert size(arr1) == size(arr2)
    @assert nfac1.num == 1 && nfac2.num == 1
    den1, den2 = nfac1.den, nfac2.den
    comm_den = lcm(den1, den2)
    new_arr = arr1 * div(comm_den, den1) + arr2 * div(comm_den, den2)
    return arr_relprime(new_arr, 1//comm_den)
end

# NC: Number of contracted legs, N: Number of incoming spaces of FTree
function fsymiofunc_step1(::Val{NC}, 
    FTree::FTree{S, NI}, 
    rem::NTuple{M, NTuple{NZ, Int}}, 
    i) where {S<:NonabelianSymm, NI, NC, M, NZ}
    @assert NI - M == 3
    outsp = i==1 ? FTree.outs[1] : rem[i-1]
    in1 = NC-1>M ? FTree.ins[1] : rem[NC-1]
    in2 = i==NC-1 ? FTree.ins[NI-NC+1] : rem[i]
    in3 = FTree.ins[NI-i+1]
    return in1, in2, in3, outsp
end


function get_perm_step1(n::Int, ci::NTuple{M, Int}) where {M}
    @assert M <= n
    all_nums, excluded = Set(1:n), Set(ci)
    front = sort(collect(setdiff(all_nums, excluded)))
    return (front..., ci...)
end

function getNsave_Xsymbol_1j(::Type{S},
    up1sp::NTuple{U1, NTuple{NZ, Int}},
    dn1sp::NTuple{D1, NTuple{NZ, Int}},
    up2sp::NTuple{U2, NTuple{NZ, Int}},
    dn2sp::NTuple{D2, NTuple{NZ, Int}},
    ctlegs1::NTuple{M, Int},
    ctlegs2::NTuple{M, Int},
    is1j_1::Bool,
    is1j_2::Bool,
    CGTom_res::CGTom{S, NZ};
    verbose) where {S<:NonabelianSymm, NZ, U1, D1, U2, D2, M}

    l1up, l1dn, l2up, l2dn = separate_legs(ctlegs1, ctlegs2, U1, D1, U2, D2)
    @assert length(l1up) == length(l2dn); 
    @assert length(l1dn) == length(l2up); 
    # TODO: Test when both CGTs are 1j-symbol
    if is1j_1 && is1j_2
        if length(up1sp) == 2 
            upsp, dnsp = up1sp, dn2sp
            upct, dnct = ctlegs1, Tuple(x-1 for x in ctlegs2)
        else 
            upsp, dnsp = up2sp, dn1sp 
            upct, dnct = ctlegs2, Tuple(x-1 for x in ctlegs1)
        end

        up1j = getNsave_1jsym(S, BigInt, BigInt, upsp[1])
        mw = qlab2mwz(S, upsp[1]); mwkey1 = (mw, .-mw)
        spdim = dimension(S, upsp[1])
        dn1j = getNsave_1jsym(S, BigInt, BigInt, dnsp[1])

        mwblk_up1j = up1j.blocks[mwkey1]
        @assert size(mwblk_up1j) == (1, 1, 1)
        sign1 = sign(mwblk_up1j[1, 1, 1])

        cw = mwkey1[upct[1]] # contracted weight
        mwkey2 = dnct[1] == 1 ? (cw, .-cw) : (.-cw, cw)
        mwblk_dn1j = dn1j.blocks[mwkey2]
        @assert size(mwblk_dn1j) == (1, 1, 1)
        sign2 = sign(mwblk_dn1j[1, 1, 1])

        elem = sign1 * sign2 * (M == 2 ? 1.0 : 1/sqrt(spdim))
        result = fill(elem, 1, 1, 1)
        return result
    else 
        mfac = 1 # Multiplicative factor when permuting 1j-symbol
        # Only one of the two CGTs is 1j-symbol
        @assert is1j_1 || is1j_2
        if is1j_1 # The first CGT is 1j-symbol
            up1jsp, dn1jsp = up1sp, dn1sp
            upcgtsp, dncgtsp = up2sp, dn2sp
            up1jct, upcgtct = l1up, l2up
            dn1jct, dncgtct = l1dn, l2dn
        else # The second CGT is 1j-symbol
            up1jsp, dn1jsp = up2sp, dn2sp
            upcgtsp, dncgtsp = up1sp, dn1sp
            up1jct, upcgtct = l2up, l1up
            dn1jct, dncgtct = l2dn, l1dn
        end
        U, D = length(upcgtsp), length(dncgtsp)

        # Permute the 1j-symbol if needed
        # After permutation, always the first leg of 1j-symbol is contracted
        incom_1j = !isempty(up1jct)
        if incom_1j 
            @assert length(up1jct) == 1 && isempty(dn1jct)
            if up1jct[1] == 2
                q = up1jsp[2]; mfac *= get_Rsym_sign(S, q)
                up1jct = (1,); up1jsp = (up1jsp[2], up1jsp[1])
            end
        else 
            @assert length(dn1jct) == 1 && isempty(up1jct) 
            if dn1jct[1] == 2
                q = dn1jsp[2]; mfac *= get_Rsym_sign(S, q)
                dn1jct = (1,); dn1jsp = (dn1jsp[2], dn1jsp[1])
            end
        end
        cgtsp_cont = incom_1j ? dncgtsp : upcgtsp
        cgtct = incom_1j ? dncgtct : upcgtct
        @assert length(cgtct) == 1; ctleg_cgt = cgtct[1]
        contracted_sp = cgtsp_cont[ctleg_cgt]

        # bp: before processing
        CGT_oms = get_CGTom(S, upcgtsp, dncgtsp, false)
        CGT_FTrees_up_bp = Dict{NTuple{NZ, Int}, Vector{FTree{S, U, NZ}}}()
        CGT_FTrees_dn_bp = Dict{NTuple{NZ, Int}, Vector{FTree{S, D, NZ}}}()

        # Prepare unit FTrees, csp: central space
        for csp in CGT_oms.spaces
            fill_FTrees!(CGT_FTrees_up_bp, upcgtsp, csp)
            fill_FTrees!(CGT_FTrees_dn_bp, dncgtsp, csp)
        end

        ct_FTrees = incom_1j ? CGT_FTrees_dn_bp : CGT_FTrees_up_bp
        other_FTrees = incom_1j ? CGT_FTrees_up_bp : CGT_FTrees_dn_bp
        combined_res = combine_FTreevecs_1j(S, ct_FTrees, other_FTrees, 
        contracted_sp, CGT_oms, ctleg_cgt, mfac, incom_1j)

        # We got vector of Step2_result, need to permute it to make it canonical form
        permute_io!(combined_res, incom_1j, is1j_1; verbose)

        om = CGT_oms.totalOM
        result = Matrix{Float64}(undef, om, om)
        for i in 1:om
            result[i, :] = to_vector(combined_res[i], CGTom_res; verbose)
        end
        final_size = is1j_1 ? (1, om, om) : (om, 1, om)
        result = reshape(result, final_size)
        return result
    end
end

function combine_FTreevecs_1j(::Type{S}, 
    ct_FTrees::Dict{NTuple{NZ, Int}, Vector{FTree{S, NC, NZ}}},
    other_FTrees::Dict{NTuple{NZ, Int}, Vector{FTree{S, NO, NZ}}},
    contracted_sp::NTuple{NZ, Int},
    CGT_oms::CGTom{S, NZ},
    ctleg_cgt::Int,
    mfac::CT,
    incom_1j::Bool) where {S<:NonabelianSymm, NZ, NC, NO, CT<:Integer}

    NC_ = max(2, NC)
    ct_FTrees_processed = Dict{NTuple{NZ, Int}, Vector{FTree{S, NC_, NZ}}}()
    perm = get_perm_step1(NC, (ctleg_cgt,))
    for (sp, ftree_vec) in ct_FTrees
        ct_FTrees_processed[sp] = Vector{FTree{S, NC_, NZ}}(undef, length(ftree_vec))
        for i in 1:length(ftree_vec)
            # First, permute the leg so that the contrated leg is at the rightmost position
            result = permute!(ftree_vec[i], perm)
            if NC == 1 result = add_singdim(result, true) end
            ct_FTrees_processed[sp][i] = result
            @assert result.ins[end] == contracted_sp
        end
    end

    NU = incom_1j ? NO+1 : NC_-1
    ND = incom_1j ? NC_-1 : NO+1
    om = CGT_oms.totalOM
    results = Vector{Step2_result{S, NU, ND, NZ}}(undef, om)
    for i in 1:om
        csp, upidx, dnidx = getominfo(CGT_oms, i)
        ctidx = incom_1j ? dnidx : upidx
        othidx = incom_1j ? upidx : dnidx
        ct_ftree = ct_FTrees_processed[csp][ctidx]
        oth_ftree = other_FTrees[csp][othidx]
        results[i] = combine_FTrees_1j(S, ct_ftree, oth_ftree, csp, mfac, incom_1j)
    end
    return results
end

function combine_FTrees_1j(::Type{S},
    ct_FTree::FTree{S, NC, NZ},
    oth_FTree::FTree{S, NO, NZ},
    csp::NTuple{NZ, Int},
    mfac::CT,
    incom_1j::Bool) where {S<:NonabelianSymm, NZ, NC, NO, CT<:Integer}

    @assert NC >= 2 && NO >= 1
    NKL = NC + NO - 3
    new_coeff = Dict{NTuple{NKL, NTuple{NZ, Int}}, Array{BigInt, NKL+1}}()
    new_nfac = Dict{NTuple{NKL, NTuple{NZ, Int}}, Rational{BigInt}}()
    for k1 in keys(ct_FTree.coeff)
        arr1, nfac1 = ct_FTree.coeff[k1], ct_FTree.coeff_nfac[k1]
        in1 = NC==2 ? ct_FTree.ins[1] : k1[1]
        in2, out = ct_FTree.ins[end], csp
        CG3f = getNsave_cg3flip(S, BigInt, BigInt, (in1, in2), out)
        flip_mat, flip_nfac::Integer = cg3flip_rightnormalized(CG3f)
        arr1_after = contract_ith(arr1, flip_mat, Val(1))
        nfac1_after = nfac1 // flip_nfac

        for k2 in keys(oth_FTree.coeff)
            arr2, nfac2 = oth_FTree.coeff[k2], oth_FTree.coeff_nfac[k2]
            k1_vec, k2_vec = collect(k1), collect(k2)
            if isempty(k1_vec) k1_vec = NTuple{NZ, Int}[] end
            if isempty(k2_vec) k2_vec = NTuple{NZ, Int}[] end
            nkey = get_newkey_X1j(NC, NO, k1_vec, k2_vec, csp, incom_1j)
            narr, nfac_final = get_newarr_X1j(arr1_after, arr2, nfac1_after, nfac2, Val(incom_1j))
            
            @assert !haskey(new_coeff, nkey)
            @assert !haskey(new_nfac, nkey)
            new_coeff[nkey] = narr*mfac; new_nfac[nkey] = nfac_final
        end
    end
    nin, nout = get_newio_X1j(S, ct_FTree.ins, oth_FTree.ins, incom_1j)
    return create_s2res(S, nin, nout, new_coeff, new_nfac)
end

get_transpose_perm(i::Int) = transpose(reshape(1:i^2, i, i))[:]

function get_conj_perm(cgtom::CGTom)
    perm = Vector{Int}(undef, cgtom.totalOM)
    si = 0
    for (i, j) in cgtom.FTree_oms
        @assert i == j
        perm[si+1:si+i^2] = get_transpose_perm(i) .+ si
        si += i^2
    end
    return perm
end

function get_newkey_X1j(nc::Int,
    no::Int,
    kc::Vector{NTuple{NZ, Int}},
    ko::Vector{NTuple{NZ, Int}},
    csp::NTuple{NZ, Int},
    incom_1j::Bool) where {NZ}

    nkey = Vector{NTuple{NZ, Int}}()
    if incom_1j
        if no >= 2 && no + nc > 3 push!(nkey, csp) end
        append!(nkey, ko); append!(nkey, kc)
    else
        append!(nkey, kc[2:end])
        if !isempty(kc) push!(nkey, kc[1]) end
        if no >= 2 && nc + no > 3 push!(nkey, csp) end
        append!(nkey, ko)
    end
    @assert length(nkey) == nc + no - 3
    return Tuple(nkey)
end

function get_newarr_X1j(arrc::Array{BigInt, NC},
    arro::Array{BigInt, NO},
    nfacc::Rational{BigInt},
    nfaco::Rational{BigInt},
    ::Val{B}) where {NC, NO, B}

    @assert nfacc.num == 1 && nfaco.num == 1
    product = get_product(arrc, arro, Val(B))
    fac_new = nfacc * nfaco
    return arr_relprime(product, fac_new)
end

@generated function get_product(
    arrc::Array{BigInt, NC},
    arro::Array{BigInt, NO},
    ::Val{B}) where {NC, NO, B}

    indc = [Symbol(:ic, i) for i in 1:NC]
    indo = [Symbol(:io, i) for i in 1:NO]
    out_ind = B ? [indc[1], indo..., indc[2:end]...] : 
                  [indc[2:end]..., indc[1], indo...]
    
    quote
        @tensor result[$(out_ind...)] := arrc[$(indc...)] * arro[$(indo...)]
        return result
    end
end


function get_newio_X1j(::Type{S},
    ct_ins::NTuple{NC, NTuple{NZ, Int}},
    oth_ins::NTuple{NO, NTuple{NZ, Int}},
    incom_1j::Bool) where {S<:NonabelianSymm, NZ, NC, NO}

    new_in = incom_1j ? [oth_ins...] : [ct_ins[1:end-1]...]
    new_out = incom_1j ? [ct_ins[1:end-1]...] : [oth_ins...]
    q = ct_ins[end]; dualq = get_dualq(S, q)
    if incom_1j push!(new_in, dualq) else push!(new_out, dualq) end
    return Tuple(new_in), Tuple(new_out)
end

function permute_io!(res_vec::Vector{Step2_result{S, NU, ND, NZ}},
    incom_1j::Bool,
    is1j_1::Bool;
    verbose=0) where {S<:NonabelianSymm, NZ, NU, ND}

    ex = res_vec[1]
    tosort = incom_1j ? collect(ex.ins) : collect(ex.outs)
    len = length(tosort); q = tosort[end]
    @assert issorted(tosort[1:end-1])
    if is1j_1
        i = searchsortedfirst(tosort[1:end-1], q)
        perm = [1:i-1..., len, i:len-1...]
    else
        i = searchsortedlast(tosort[1:end-1], q)
        perm = [1:i..., len, i+1:len-1...]
    end
    if issorted(perm) return end
    for i in 1:length(res_vec)
        res_vec[i] = permute!(res_vec[i], Tuple(perm), incom_1j; verbose)
    end
end

function get_Rsym_sign(::Type{S},
    q::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ}

    dualq = get_dualq(S, q)
    if q != dualq return 1 end
    r = getNsave_Rsymbol(S, BigInt, q, Tuple(0 for _=1:NZ))
    @assert size(r.rsym_mat) == (1, 1)
    @assert abs(r.rsym_mat[1, 1]) == r.nfactor[1].den
    return sign(r.rsym_mat[1, 1])
end

# Tensor operations for contracting CG3s
@generated function contract_newcg3(before_step::SparseArray{FT},
    last_cg3::SparseArray{FT}, ::Val{N}) where {N, FT<:AbstractFloat}
    @assert N >= 2
    # a: the leg to be contracted, o: output leg
    # o1, o2, o3, ...: outer multiplicity legs
    before_inds = [:a, [Symbol(:i, i) for i=3:N]..., :o, [Symbol(:o, i) for i in 1:N-2]...]
    cg3_inds = [:i1, :i2, :a, Symbol(:o, N-1)]
    result_inds = [[Symbol(:i, i) for i=1:N]..., :o, [Symbol(:o, i) for i in 1:N-1]...]

    quote
        @tensor new_cg[$(result_inds...)] := before_step[$(before_inds...)] * last_cg3[$(cg3_inds...)]
        return new_cg
    end
end

@generated function contract_om(contract_res, om_arr, ::Val{N}) where N
    inout_inds = [[Symbol(:i, i) for i in 1:N]..., :o]
    om_inds = [Symbol(:o, i) for i in 1:N-1]
    before_inds = vcat(inout_inds, om_inds)

    quote
        @tensor result[$(inout_inds...)] := contract_res[$(before_inds...)] * om_arr[$(om_inds...)]
        return result
    end
end

# Contract CG3s into a single tensor.
function contract_CG3s(::Type{S},
    ins::Vector,
    out,
    intsps,
    ::Val{N},
    ::Type{FT},
    normalize=false) where {S<:NonabelianSymm, N, FT<:AbstractFloat}
    @assert N >= 2
    @assert length(ins) == N
    @assert length(intsps) == N - 2
    if N == 2
        # Just CG3 case
        return load_cg3_float(S, BigInt, Tuple(ins), out, FT, normalize)[1]::SparseArray{FT}
    end
    ins_before = copy(ins[2:N])
    ins_before[1] = intsps[end]
    before_step = contract_CG3s(S, ins_before, out, intsps[1:end-1], Val(N-1), FT, normalize)
    last_cg3, _, _ = load_cg3_float(S, BigInt, (ins[1], ins[2]), intsps[end], FT, normalize)
    return contract_newcg3(before_step, last_cg3, Val(N))::SparseArray{FT}
end

# Get array from FTree
function FTree2arr(ftree::FTree{S, N},
    ::Type{FT},
    normalize=false) where {S<:NonabelianSymm, N, FT<:AbstractFloat}

    szs = Int[]
    for insp in ftree.ins
        irep = getNsave_irep(S, BigInt, insp)
        push!(szs, dimension(irep))
    end
    out_irep = getNsave_irep(S, BigInt, ftree.outs[1])
    push!(szs, dimension(out_irep))
    if N == 1
        @assert length(szs) == 2 && szs[1] == szs[2]; sz = szs[1]
        elem = haskey(ftree.coeff, ()) ? FT(ftree.coeff[()][] * ftree.coeff_nfac[()]) : 0.0
        return SparseArray(Matrix(I, sz, sz) * elem)::SparseArray{FT}
    end

    arr = SparseArray(zeros(szs...))
    for intsps in keys(ftree.coeff)
        @assert length(intsps) == N - 2
        cont_res = contract_CG3s(S, collect(ftree.ins), ftree.outs[1], intsps, Val(N), FT, normalize)
        om_arr = SparseArray{FT}(ftree.coeff[intsps] * ftree.coeff_nfac[intsps])
        arr += contract_om(cont_res, om_arr, Val(N))
    end
    return arr
end

function contract_arrs(arr1::AbstractArray,
    arr2::AbstractArray,
    c1::NTuple{M, Int},
    c2::NTuple{M, Int}) where M

    # keeped axes
    keep1 = setdiff(1:ndims(arr1), c1)
    keep2 = setdiff(1:ndims(arr2), c2)

    # size check
    @assert all(size(arr1, c1i) == size(arr2, c2i) for (c1i, c2i) in zip(c1, c2)) "size mismatch"

    # reorder the legs of arr1
    perm1 = [collect(keep1)..., c1...]
    arr1p = permutedims(arr1, perm1)

    # reorder the legs of arr2
    perm2 = [c2..., collect(keep2)...]
    arr2p = permutedims(arr2, perm2)

    # reshape arrays
    arr1mat = reshape(arr1p, :, prod(size(arr1p)[length(keep1)+1:end]))
    arr2mat = reshape(arr2p, prod(size(arr2p)[1:length(c2)]), :)

    # perform matrix multiplication
    result_mat = arr1mat * arr2mat

    # Size of the resulting array
    result_shape = [size(arr1p)[1:length(keep1)]..., size(arr2p)[length(c2)+1:end]...]
    return reshape(result_mat, result_shape...)
end

function get_canonical_basis(::Type{S},
    insp::NTuple{NI, NTuple{NZ, Int}},
    outsp::NTuple{NO, NTuple{NZ, Int}},
    om::CGTom{S, NZ};
    verbose=0) where {S<:NonabelianSymm, NI, NO, NZ}

    result = Vector{SparseArray{Float64, NI+NO}}()
    totalOM = om.totalOM
    for i in 1:totalOM
        csp, upidx, dnidx = getominfo(om, i)
        omlist_up = getNsave_omlist(S, insp, csp)
        omlist_dn = getNsave_omlist(S, outsp, csp)
        CGTup = create_unit_FTree(omlist_up, upidx)
        CGTdn = create_unit_FTree(omlist_dn, dnidx)

        CGTup_arr = FTree2arr(CGTup, Float64, true)
        CGTdn_arr = FTree2arr(CGTdn, Float64, true)
        CGT = contract_arrs(CGTup_arr, CGTdn_arr, (length(insp)+1,), (length(outsp)+1,))
        dim_csp = dimension(getNsave_irep(S, BigInt, csp))
        CGT /= sqrt(dim_csp); @assert norm(CGT) ≈ 1
        push!(result, CGT)
    end
    return result
end