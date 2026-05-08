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


