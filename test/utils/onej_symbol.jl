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

