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

