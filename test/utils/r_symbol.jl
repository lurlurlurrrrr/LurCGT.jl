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

function get_random_Rsymbol_real(::Type{S}, ::Type{FT}, qlimit=4) where {S<:NonabelianSymm, FT<:AbstractFloat}
    NZ = nzops(S)
    qlabel = Tuple(rand(0:qlimit) for _=1:NZ)
    vo = getNsave_validout(S, (qlabel, qlabel))
    out = select_out(vo, 100)
    return qlabel, out, Rsymbol_toreal(S, BigInt, qlabel, out)
end

