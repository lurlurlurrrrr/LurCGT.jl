function test_CGTSVD_identity(::Type{S},
    upsp,
    dnsp;
    verbose=0) where {S<:NonabelianSymm}

    upsp = Tuple(sort(collect(upsp)))
    dnsp = Tuple(sort(collect(dnsp)))
    obj = getNsave_CGTSVD(S, upsp, dnsp, Tuple(1:length(upsp)); save=false, verbose)
    CGTom = get_CGTom(S, upsp, dnsp)
    om = CGTom.totalOM

    @assert size(obj.svd_arr) == (om, om)
    @assert obj.svd_arr ≈ Matrix{Float64}(I, om, om)

    expected = [(CGTom.spaces[i], CGTom.FTree_oms[i][1], CGTom.FTree_oms[i][2]) for i in eachindex(CGTom.spaces)]
    @assert obj.bond_sps == expected
end

function test_CGTSVD_unitary(::Type{S},
    upsp,
    dnsp,
    leftlegs;
    verbose=0) where {S<:NonabelianSymm}

    upsp = Tuple(sort(collect(upsp)))
    dnsp = Tuple(sort(collect(dnsp)))
    obj = getNsave_CGTSVD(S, upsp, dnsp, leftlegs; save=false, verbose)
    om = get_CGTom(S, upsp, dnsp).totalOM
    eye = Matrix{Float64}(I, om, om)

    @assert sum(omL * omR for (_, omL, omR) in obj.bond_sps) == om
    @assert obj.svd_arr * obj.svd_arr' ≈ eye
    @assert obj.svd_arr' * obj.svd_arr ≈ eye
end

function dense_copy(arr)
    reshape(collect(arr), size(arr)...)
end

function normalize_zero_space_axes(::Type{S}, insp, outsp, arr) where {S<:NonabelianSymm}
    zeroq = Tuple(0 for _ in 1:nzops(S))
    insp_vec = collect(insp)
    outsp_vec = collect(outsp)
    arr_dense = dense_copy(arr)

    keep_in = findall(!=(zeroq), insp_vec)
    keep_out = findall(!=(zeroq), outsp_vec)

    if isempty(keep_in)
        if !isempty(insp_vec)
            keep_in = [1]
        end
        new_insp = (zeroq,)
    else
        new_insp = Tuple(insp_vec[keep_in])
    end

    if isempty(keep_out)
        if !isempty(outsp_vec)
            keep_out = [1]
        end
        new_outsp = (zeroq,)
    else
        new_outsp = Tuple(outsp_vec[keep_out])
    end

    nin = length(insp_vec)
    nout = length(outsp_vec)
    drop_axes = vcat(setdiff(collect(1:nin), keep_in), nin .+ setdiff(collect(1:nout), keep_out))
    if !isempty(drop_axes)
        arr_dense = dropdims(arr_dense; dims=Tuple(drop_axes))
    end

    if nin == 0
        arr_dense = reshape(arr_dense, 1, size(arr_dense)...)
    end
    if nout == 0
        arr_dense = reshape(arr_dense, size(arr_dense)..., 1)
    end

    return new_insp, new_outsp, arr_dense
end

function get_final_CGTinfo_stable(::Type{S},
    cgt1_insp::NTuple{NI1, NTuple{NZ, Int}},
    cgt1_outsp::NTuple{NO1, NTuple{NZ, Int}},
    cgt2_insp::NTuple{NI2, NTuple{NZ, Int}},
    cgt2_outsp::NTuple{NO2, NTuple{NZ, Int}},
    cgt1legs::NTuple{M, Int},
    cgt2legs::NTuple{M, Int},
    arr3_cont) where {S<:NonabelianSymm, NI1, NO1, NI2, NO2, M, NZ}

    cgt1_inopen = [cgt1_insp[i] for i in 1:NI1 if !(i in cgt1legs)]
    cgt1_outopen = [cgt1_outsp[i] for i in 1:NO1 if !(i+NI1 in cgt1legs)]

    cgt2_inopen = [cgt2_insp[i] for i in 1:NI2 if !(i in cgt2legs)]
    cgt2_outopen = [cgt2_outsp[i] for i in 1:NO2 if !(i+NI2 in cgt2legs)]

    cgt3_insp = vcat(cgt1_inopen, cgt2_inopen)
    cgt3_outsp = vcat(cgt1_outopen, cgt2_outopen)

    n1ino, n1outo = length(cgt1_inopen), length(cgt1_outopen)
    n2ino, n2outo = length(cgt2_inopen), length(cgt2_outopen)

    perm1 = Tuple(vcat(collect(1:n1ino),
        collect(n1ino+n1outo+1:n1ino+n1outo+n2ino),
        collect(n1ino+1:n1ino+n1outo),
        collect(n1ino+n1outo+n2ino+1:n1ino+n1outo+n2ino+n2outo)))

    #println("perm1: ", perm1)
    arr3_cont = permutedims(arr3_cont, perm1)
    in_sortperm = sortperm(cgt3_insp; alg=MergeSort)
    out_sortperm = sortperm(cgt3_outsp; alg=MergeSort)
    perm2 = Tuple(vcat(in_sortperm, out_sortperm .+ length(cgt3_insp)))
    arr3_cont = permutedims(arr3_cont, perm2)
    #println("perm2: ", perm2)

    cgt3_insp = cgt3_insp[in_sortperm]
    cgt3_outsp = cgt3_outsp[out_sortperm]
    cgt3_insp_t, cgt3_outsp_t, arr3_cont =
        normalize_zero_space_axes(S, Tuple(cgt3_insp), Tuple(cgt3_outsp), arr3_cont)
    CGT3om = get_CGTom(S, cgt3_insp_t, cgt3_outsp_t)

    return cgt3_insp_t, cgt3_outsp_t, CGT3om, arr3_cont
end

function get_random_CGTSVD_input(::Type{S},
    dim_limit=100000,
    qlimit=3;
    verbose=0) where {S<:NonabelianSymm}

    while true
        upsp, dnsp = getRandomCGT(S, dim_limit, qlimit; verbose)
        upsp = Tuple(sort(collect(upsp)))
        dnsp = Tuple(sort(collect(dnsp)))
        upsp, _ = remove_zeros(S, upsp)
        dnsp, _ = remove_zeros(S, dnsp)
        if length(upsp) + length(dnsp) < 2
            continue
        end
        cgtom = get_CGTom(S, upsp, dnsp)
        if cgtom.totalOM > 0
            return upsp, dnsp
        end
    end
end

function get_random_sameq_CGTSVD_input(::Type{S},
    dim_limit=100000,
    qlimit=3;
    verbose=0) where {S<:NonabelianSymm}

    NZ = nzops(S)
    zeroq = Tuple(0 for _ in 1:NZ)
    while true
        qrep = Tuple(rand(0:qlimit) for _ in 1:NZ)
        qrep == zeroq && continue

        rep_dim = LurCGT.dimension(getNsave_irep(S, BigInt, qrep))
        nlegs_max = 0
        total_dim = 1
        while nlegs_max < 8 && total_dim * rep_dim <= dim_limit
            total_dim *= rep_dim
            nlegs_max += 1
        end
        nlegs_max < 3 && continue

        nlegs = rand(3:nlegs_max)
        nup = rand(1:nlegs-1)
        ndn = nlegs - nup
        upsp = ntuple(_ -> qrep, nup)
        dnsp = ntuple(_ -> qrep, ndn)
        cgtom = get_CGTom(S, upsp, dnsp)
        if cgtom.totalOM > 0
            return upsp, dnsp, qrep
        end
    end
end

function get_random_leftlegs(nlegs::Int)
    perm = shuffle(collect(1:nlegs))
    nleft = rand(1:nlegs-1)
    return Tuple(sort(perm[1:nleft]))
end

function get_random_sameq_leftlegs(nup::Int, ndn::Int)
    @assert nup + ndn >= 3
    @assert nup > 1 || ndn > 1

    while true
        leftlegs = get_random_leftlegs(nup + ndn)
        nup_left = count(i -> i <= nup, leftlegs)
        ndn_left = length(leftlegs) - nup_left
        if (nup > 1 && 0 < nup_left < nup) || (ndn > 1 && 0 < ndn_left < ndn)
            return leftlegs
        end
    end
end

function get_CGTSVD_split_basis_direct(::Type{S},
    upsp::NTuple{U, NTuple{NZ, Int}},
    dnsp::NTuple{D, NTuple{NZ, Int}},
    leftlegs::NTuple{L, Int};
    verbose=0) where {S<:NonabelianSymm, U, D, L, NZ}

    obj = getNsave_CGTSVD(S, upsp, dnsp, leftlegs; save=false, verbose)
    zeroq = Tuple(0 for _ in 1:NZ)
    left_up, left_dn, right_up, right_dn = LurCGT.get_split_side_spaces(S, upsp, dnsp, leftlegs)

    cgtom = get_CGTom(S, upsp, dnsp)
    canbasis = get_canonical_basis(S, upsp, dnsp, cgtom; verbose)
    basis_shape = size(canbasis[1])
    splitbasis = Array{Float64, U + D}[]

    for (q, omL, omR) in obj.bond_sps
        dualq = get_dualq(S, q)
        center_up = LurCGT.stable_sort_tuple((q, dualq))
        center_dn = (zeroq,)
        left_dn_q = LurCGT.stable_sort_tuple((left_dn..., q))
        right_dn_dualq = LurCGT.stable_sort_tuple((dualq, right_dn...))
        interm_up = LurCGT.stable_sort_tuple((dualq, left_up...))
        interm_dn = left_dn

        qpos_center = findfirst(==(q), center_up)
        qpos_left = findlast(==(q), left_dn_q)
        dualqpos_interm = findfirst(==(dualq), interm_up)
        dualqpos_right = findfirst(==(dualq), right_dn_dualq)
        @assert !isnothing(qpos_center)
        @assert !isnothing(qpos_left)
        @assert !isnothing(dualqpos_interm)
        @assert !isnothing(dualqpos_right)

        center_ctlegs = (something(qpos_center),)
        left_ctlegs = (length(left_up) + something(qpos_left),)
        interm_ctlegs = (something(dualqpos_interm),)
        right_ctlegs = (length(right_up) + something(dualqpos_right),)

        center_up_s, center_dn_s, center_ctlegs_s =
            LurCGT.standardize_spaces_and_legs(S, center_up, center_dn, center_ctlegs, true)
        left_up_s, left_dn_q_s, left_ctlegs_s =
            LurCGT.standardize_spaces_and_legs(S, left_up, left_dn_q, left_ctlegs, false)
        interm_up_s, interm_dn_s, interm_ctlegs_s =
            LurCGT.standardize_spaces_and_legs(S, interm_up, interm_dn, interm_ctlegs, true)
        right_up_s, right_dn_dualq_s, right_ctlegs_s =
            LurCGT.standardize_spaces_and_legs(S, right_up, right_dn_dualq, right_ctlegs, false)

        center_om = get_CGTom(S, center_up_s, center_dn_s, detect_1j(S, center_up_s, center_dn_s))
        left_om = get_CGTom(S, left_up_s, left_dn_q_s, detect_1j(S, left_up_s, left_dn_q_s))
        right_om = get_CGTom(S, right_up_s, right_dn_dualq_s, detect_1j(S, right_up_s, right_dn_dualq_s))
        @assert center_om.totalOM == 1
        @assert left_om.totalOM == omL
        @assert right_om.totalOM == omR

        center_basis = get_canonical_basis(S, center_up_s, center_dn_s, center_om; verbose)
        #println("center: $(center_up_s)->$(center_dn_s), om=$(center_om.totalOM)")
        left_basis = get_canonical_basis(S, left_up_s, left_dn_q_s, left_om; verbose)
        #println("left: $(left_up_s)->$(left_dn_q_s), om=$(left_om.totalOM)")
        right_basis = get_canonical_basis(S, right_up_s, right_dn_dualq_s, right_om; verbose)
        #println("right: $(right_up_s)->$(right_dn_dualq_s), om=$(right_om.totalOM)")
        dimq = Float64(LurCGT.dimension(getNsave_irep(S, BigInt, q)))

        for b in 1:omR
            for a in 1:omL
                #println("contracted legs: center leg $(center_ctlegs_s[1]), left leg $(left_ctlegs_s[1])")
                arr12_cont = contract_arrs(center_basis[1], left_basis[a], center_ctlegs_s, left_ctlegs_s)
                interm_up_chk, interm_dn_chk, _, arr12_cont = get_final_CGTinfo_stable(
                    S, center_up_s, center_dn_s, left_up_s, left_dn_q_s,
                    center_ctlegs_s, left_ctlegs_s, arr12_cont)
                @assert interm_up_chk == interm_up_s
                @assert interm_dn_chk == interm_dn_s

                #println("contracted legs: interm leg $(interm_ctlegs_s[1]), right leg $(right_ctlegs_s[1])")
                arr123_cont = contract_arrs(arr12_cont, right_basis[b], interm_ctlegs_s, right_ctlegs_s)
                final_up, final_dn, _, arr123_cont = get_final_CGTinfo_stable(
                    S, interm_up_s, interm_dn_s, right_up_s, right_dn_dualq_s,
                    interm_ctlegs_s, right_ctlegs_s, arr123_cont)

                @assert final_up == upsp
                @assert final_dn == dnsp
                final_perm = LurCGT.get_cgtsvd_final_perm(S, upsp, dnsp, leftlegs)
                #println("perm: $(invperm(final_perm))")
                if final_perm != Tuple(1:length(final_perm))
                    arr123_cont = permutedims(arr123_cont, invperm(final_perm))
                end
                push!(splitbasis, reshape(dimq .* dense_copy(arr123_cont), basis_shape...))
            end
        end
    end

    @assert length(splitbasis) == cgtom.totalOM
    return splitbasis
end

function test_CGTSVD_unitary_randinput(::Type{S},
    dim_limit=10000,
    qlimit=3,
    ninput=1000;
    verbose=0) where {S<:NonabelianSymm}

    for i in 1:ninput
        upsp, dnsp = get_random_CGTSVD_input(S, dim_limit, qlimit; verbose)
        leftlegs = get_random_leftlegs(length(upsp) + length(dnsp))
        println("Test input #$i: CGT $(upsp)->$(dnsp), leftlegs=$(leftlegs)")

        obj = getNsave_CGTSVD(S, upsp, dnsp, leftlegs; save=false, verbose)
        om = get_CGTom(S, upsp, dnsp).totalOM
        eye = Matrix{Float64}(I, om, om)

        @assert size(obj.svd_arr) == (om, om)
        @assert sum(omL * omR for (_, omL, omR) in obj.bond_sps) == om
        @assert obj.svd_arr * obj.svd_arr' ≈ eye
        @assert obj.svd_arr' * obj.svd_arr ≈ eye
        println("Test passed")
    end
end

function check_CGTSVD_basischange_case(::Type{S},
    upsp,
    dnsp,
    leftlegs;
    verbose=0) where {S<:NonabelianSymm}

    obj = getNsave_CGTSVD(S, upsp, dnsp, leftlegs; save=false, verbose)
    cgtom = get_CGTom(S, upsp, dnsp)
    om = cgtom.totalOM
    canbasis = get_canonical_basis(S, upsp, dnsp, cgtom; verbose)
    splitbasis = get_CGTSVD_split_basis_direct(S, upsp, dnsp, leftlegs; verbose)

    canbasis_dense = [dense_copy(canbasis[j]) for j in 1:om]
    svd_direct = zeros(Float64, om, om)
    for aidx in 1:om
        for bidx in 1:om
            svd_direct[aidx, bidx] = sum(splitbasis[aidx] .* canbasis_dense[bidx])
        end
    end

    ndiff = norm(svd_direct - obj.svd_arr) / norm(obj.svd_arr)
    @assert ndiff < 1e-10 "CGTSVD basis change mismatch: $(ndiff)"
end

function test_CGTSVD_duplicate_legorder(::Type{SU{2}}; verbose=0)
    upsp = ((1,), (1,), (2,))
    dnsp = ((1,), (1,))
    leftlegs = (2, 4)
    @assert get_CGTom(SU{2}, upsp, dnsp).totalOM > 0
    println("Duplicate-label test input: CGT $(upsp)->$(dnsp), leftlegs=$(leftlegs)")
    check_CGTSVD_basischange_case(SU{2}, upsp, dnsp, leftlegs; verbose)
    println("Test passed")
end

function test_CGTSVD_sameq_randinput(::Type{S},
    dim_limit=10000,
    qlimit=3,
    ninput=1000;
    verbose=0) where {S<:NonabelianSymm}

    for i in 1:ninput
        upsp, dnsp, qrep = get_random_sameq_CGTSVD_input(S, dim_limit, qlimit; verbose)
        leftlegs = get_random_sameq_leftlegs(length(upsp), length(dnsp))
        #upsp = ((2,), (2,), (2,), (2,))
        #dnsp = ((2,), (2,), (2,))
        #leftlegs = (4,)
        
        println("Test input #$i: CGT $(upsp)->$(dnsp) with identical qlabel $(qrep), leftlegs=$(leftlegs)")
        check_CGTSVD_basischange_case(S, upsp, dnsp, leftlegs; verbose)
        println("Test passed")
    end
end

function test_CGTSVD_basischange(::Type{S},
    dim_limit=100000,
    qlimit=3,
    ninput=10;
    verbose=0) where {S<:NonabelianSymm}

    for i in 1:ninput
        upsp, dnsp = get_random_CGTSVD_input(S, dim_limit, qlimit; verbose)
        leftlegs = get_random_leftlegs(length(upsp) + length(dnsp))
        println("Test input #$i: CGT $(upsp)->$(dnsp), leftlegs=$(leftlegs)")
        check_CGTSVD_basischange_case(S, upsp, dnsp, leftlegs; verbose)
        println("Test passed")
    end
end

function getRandomCGT(::Type{S},
    dim_limit=100000,
    qlimit=3;
    verbose=0) where {S<:NonabelianSymm}

    spaces = nothing
    while isnothing(spaces) || spaces[1] == spaces[2]
        spaces = get_spaces(S, dim_limit, qlimit, 200; verbose)
    end
    # TODO: currently, get_spaces is written in old convention. Need to fix it later
    # This function is written in new convention
    cgt_in, cgt_out, _, _, _, _, _, _, _, _ = spaces

    return cgt_in, cgt_out
end

