function test_CGTQR_unitary_case(::Type{S}, upsp, dnsp, leftlegs; verbose=0) where {S<:NonabelianSymm}
    obj = getNsave_CGTQR(S, upsp, dnsp, leftlegs; save=false, verbose)
    om = get_CGTom(S, obj.upsp, obj.dnsp).totalOM
    eye = Matrix{Float64}(I, om, om)

    @assert size(obj.qr_arr) == (om, om)
    @assert sum(omL * omR for (_, omL, omR) in obj.bond_sps; init=0) == om
    @assert obj.qr_arr * transpose(obj.qr_arr) ≈ eye
    @assert transpose(obj.qr_arr) * obj.qr_arr ≈ eye
end

function test_CGTQR_unitary_examples(::Type{S}; verbose=0) where {S<:NonabelianSymm}
    test_CGTQR_unitary_case(S, ((1,),), ((1,),), (1,); verbose)
    test_CGTQR_unitary_case(S, ((0,),), ((1,), (1,)), (1, 2); verbose)
    test_CGTQR_unitary_case(S, ((1,),), ((0,), (1,)), (1, 2); verbose)
end
