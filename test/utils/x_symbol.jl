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
        #println("CGT1: $(cgt1_outsp)->$(cgt1_internalsp)->$(cgt1_insp)")
        #println("CGT2: $(cgt2_outsp)->$(cgt2_internalsp)->$(cgt2_insp)")
        #println("contracted legs: $(cgt1legs) of CGT1, $(cgt2legs) of CGT2")

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
        #println("CGT1: $(cgt1_insp)->$(cgt1_outsp)")
        #println("CGT2: $(cgt2_insp)->$(cgt2_outsp)")
        #println("contracted legs: $(cgt1legs) of CGT1, $(cgt2legs) of CGT2")

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

