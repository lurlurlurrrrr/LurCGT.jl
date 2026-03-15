struct Fsymbol{S<:NonabelianSymm, CT, NZ}
    in1::NTuple{NZ, Int}
    in2::NTuple{NZ, Int}   
    in3::NTuple{NZ, Int}
    out::NTuple{NZ, Int}
    es_list::Vector{Tuple{NTuple{NZ, Int}, Tuple{Int, Int}}}
    fs_list::Vector{Tuple{NTuple{NZ, Int}, Tuple{Int, Int}}}
    fsym_mat::Array{CT, 2}
    es_nfactor::Vector{Rational{CT}}
    fs_nfactor::Vector{Rational{CT}}
    # Cached memory footprint in bytes, computed once at construction for LRU cache eviction by size
    size_byte::Int
    
    function Fsymbol{S, CT, NZ}(in1, in2, in3, out, es_list, fs_list, fsym_mat, es_nfactor, fs_nfactor, size_byte::Int=0) where {S, CT, NZ}
        if size_byte == 0
            obj = new{S, CT, NZ}(in1, in2, in3, out, es_list, fs_list, fsym_mat, es_nfactor, fs_nfactor, 0)
            size_byte = Base.summarysize(obj)
        end
        new{S, CT, NZ}(in1, in2, in3, out, es_list, fs_list, fsym_mat, es_nfactor, fs_nfactor, size_byte)
    end
end

# Obtain the F-symbol by contracting (non-normalized) CGTs 
# maximal weight state part only
function getNsave_Fsymbol(::Type{S},
    ::Type{CT},
    in1::NTuple{NZ, Int},
    in2::NTuple{NZ, Int},
    in3::NTuple{NZ, Int},
    out::NTuple{NZ, Int};
    verbose=0) where {S<:NonabelianSymm, CT<:Number, NZ}
    @assert NZ == nzops(S)
    loaded = load_Fsymbol_sqlite(S, CT, in1, in2, in3, out)
    if !isnothing(loaded) return loaded end

    # in1 * in2 -> e, e * in3 -> out
    es_list = find_espaces(S, in1, in2, in3, out, CT; verbose)
    # in2 * in3 -> e, e * in1 -> out
    fs_list = find_espaces(S, in2, in3, in1, out, CT; verbose)

    # Get the outer multiplicity
    mat_size = 0
    for (_, (qin3_om, in12_om)) in es_list
        mat_size += in12_om * qin3_om
    end
    
    # Get the outer multiplicity in other way
    mat_size_f = 0
    for (_, (qin1_om, in23_om)) in fs_list
        mat_size_f += qin1_om * in23_om
    end
    # Outer multiplicity should be the same regardless of the order
    if verbose > 1 println("Outer multiplicity: $mat_size") end

    @assert mat_size == mat_size_f
    # Initialize the F-symbol matrix
    fsym_mat = zeros(CT, mat_size, mat_size)

    # Outer multiplicity label which is closer to the outgoing space comes first
    # ν : interm_space(e) ⊗ in3 -> out
    # μ : in1 ⊗ in2 -> interm_space(e)
    # [in1, in2, in3, ν, μ]
    es_contract_res = Dict{NTuple{3, NTuple{NZ, Int}}, Array{CT, 5}}[]
    es_nfactor = Vector{Rational{CT}}()

    cg3blks_e = load_cg3blk(S, CT, (in1, in2), [e for (e, _) in es_list])

    # Fill in the es_contract_res
    for (e, _) in es_list
        # Load CG3 in1 ⊗ in2 -> e 
        μblk, μfac = cg3blks_e[e]
        # Load CG3 e ⊗ in3 -> out
        νblk, νfac = load_cg3blk(S, CT, (e, in3), [out])[out]

        # Squeeze the 3rd dimension of νblk
        νblk_mw = mwpartof(S, νblk, (e, in3, out), 3)

        blk_contract = contractblks_mw(S, μblk, νblk_mw, true; verbose=verbose)
        push!(es_contract_res, blk_contract)
        νμfac = νfac .* μfac'
        append!(es_nfactor, νμfac[:])
    end

    # Conjugate the blocks of es_contract_res
    conjugate!(S, es_contract_res, in1, in2, in3)


    # λ : interm_space(f) ⊗ in1 -> out
    # κ : in2 ⊗ in3 -> interm_space(f)
    # [in1, in2, in3, λ, κ]
    fs_contract_res = Dict{NTuple{3, NTuple{NZ, Int}}, Array{CT, 5}}[]
    fs_nfactor = Vector{Rational{CT}}()

    cg3blks_f = load_cg3blk(S, CT, (in2, in3), [f for (f, _) in fs_list])

    for (f, _) in fs_list
        # Load CG3 in2 ⊗ in3 -> f
        κblk, κfac = cg3blks_f[f]
        # Load CG3 in1 ⊗ f -> out
        λblk, λfac = load_cg3blk(S, CT, (in1, f), [out])[out]

        # Squeeze the 3rd dimension of λblk
        λblk_mw = mwpartof(S, λblk, (in1, f, out), 3)

        blk_contract = contractblks_mw(S, κblk, λblk_mw, false; verbose=verbose)
        push!(fs_contract_res, blk_contract)
        λκfac =  λfac .* κfac'
        append!(fs_nfactor, λκfac[:])
    end

    # Contract obtained mw states and fill in the F-symbol matrix
    i, j = 1, 1
    for (k1, eblks) in enumerate(es_contract_res)
        _, (m1, m2) = es_list[k1]
        for (k2, fblks) in enumerate(fs_contract_res)
            _, (n1, n2) = fs_list[k2]
            irange = i:(i + m1*m2 - 1)
            jrange = j:(j + n1*n2 - 1)

            for zlabels in keys(eblks)
                if haskey(fblks, zlabels)
                    @tensor fsym[ν, μ, λ, κ] := 
                        eblks[zlabels][i1, i2, i3, ν, μ] * 
                        fblks[zlabels][i1, i2, i3, λ, κ]
                    fsym_mat[irange, jrange] += reshape(fsym, m1*m2, n1*n2)
                end
            end
            j += n1 * n2
        end
        i += m1 * m2
        j = 1
    end
    
    fsym_struct = Fsymbol{S, CT, NZ}(in1, in2, in3, out, es_list, fs_list, 
    fsym_mat, es_nfactor, fs_nfactor)
    save_Fsymbol_sqlite(fsym_struct)
    return fsym_struct
end

function conjugate!(::Type{S}, 
    es_contract_res::Vector{Dict{NTuple{3, NTuple{NZ, Int}}, Array{CT, 5}}}, 
    in1::NTuple{NZ, Int}, 
    in2::NTuple{NZ, Int}, 
    in3::NTuple{NZ, Int}) where {S<:NonabelianSymm, CT<:Number, NZ}
    # Conjugate the blocks of es_contract_res
    irep1 = getNsave_irep(S, BigInt, in1)
    irep2 = getNsave_irep(S, BigInt, in2)
    irep3 = getNsave_irep(S, BigInt, in3)
    for blks in es_contract_res
        for k in keys(blks)
            z1, z2, z3 = k
            m1 = irep1.innerprod[z1]
            m2 = irep2.innerprod[z2]
            m3 = irep3.innerprod[z3]
            @tensor newblk[j1, j2, j3, o1, o2] :=
                blks[k][i1, i2, i3, o1, o2] * 
                m1[i1, j1] * m2[i2, j2] * m3[i3, j3]
            blks[k] = newblk
        end
    end
end

# di: Dropped index, 3 when getting mw part for outgoing leg
# inNout: Tuple (in1, in2, out)
function mwpartof(::Type{S}, 
    blks::Dict{NTuple{3, NTuple{NZ, Int}}, Array{CT}}, 
    inNout::NTuple{3, NTuple{NZ, Int}},
    di::Int = 3) where {S<:NonabelianSymm, CT<:Number, NZ}
    # Extract the maximal weight part of the block
    mw_blk = Dict{NTuple{2, NTuple{NZ, Int}}, Array{CT, 3}}()
    out_mwzval = qlab2mwz(S, inNout[di])
    for k in keys(blks)
        @assert ndims(blks[k]) == 4
        # outgoing weight is equal to the maximal weight
        if k[di] == out_mwzval
            blk = blks[k]
            @assert size(blk, di) == 1
            newkey = Tuple(k[i] for i in 1:3 if i != di)
            mw_blk[newkey] = blk[(i==di ? 1 : Colon() for i=1:4)...]
        end
    end
    return mw_blk
end

# 3rd dimension of blk1 and 1st dimension of blk2_mw are contracted
function contractblks_mw(::Type{S},
    blks1::Dict{NTuple{3, NTuple{NZ, Int}}, Array{CT}}, 
    blks2_mw::Dict{NTuple{2, NTuple{NZ, Int}}, Array{CT, 3}}, 
    es::Bool;
    verbose=0) where {S<:NonabelianSymm, CT<:Number, NZ}

    # 3rd dimension of blk1 and 'ci'th dimension of blk2_mw are contracted
    ci = es ? 1 : 2
    contract_res = Dict{NTuple{3, NTuple{NZ, Int}}, Array{CT, 5}}()
    blks1keys = sort(collect(keys(blks1)); by=x->x[3])
    blks2keys = sort(collect(keys(blks2_mw)); by=x->x[ci])
    nkeys1, nkeys2 = length(blks1keys), length(blks2keys)

    i1, i2 = 1, 1
    while i1 <= nkeys1 && i2 <= nkeys2
        k1, k2 = blks1keys[i1], blks2keys[i2]
        if k1[3] == k2[ci]
            # If the keys match, contract the blocks
            # Find the degeneracy of the blocks deg1, deg2
            deg1, deg2 = 0, 0
            i1_, i2_ = i1, i2
            while i1_ <= nkeys1 && blks1keys[i1_][3] == k1[3]
                deg1 += 1; i1_ += 1
            end
            while i2_ <= nkeys2 && blks2keys[i2_][ci] == k2[ci]
                deg2 += 1; i2_ += 1
            end

            j1, j2 = i1, i2
            for _=1:deg1
                for _=1:deg2
                    # Contract the j1th block of blk1 and j2th block of blk2_mw
                    # and store the result in contract_res
                    key1, key2 = blks1keys[j1], blks2keys[j2]
                    blk1, blk2 = blks1[key1], blks2_mw[key2]
                    @assert ndims(blk1) == 4 && ndims(blk2) == 3

                    # Get new key. Should consider es? fs?
                    new_key = es ? (key1[1], key1[2], key2[2]) : (key2[1], key1[1], key1[2])

                    # Contract the blocks and store the result
                    cont_res = contract_blks(blk1, blk2, es)

                    if !haskey(contract_res, new_key)
                        contract_res[new_key] = cont_res
                    else
                        contract_res[new_key] += cont_res
                    end

                    j2 += 1
                end
                j2 = i2; j1 += 1
            end

            i1, i2 = i1_, i2_
        elseif k1[3] < k2[ci]
            i1 += 1
        else
            i2 += 1
        end
    end

    return contract_res
end

function contract_blks(blk1::Array{CT, 4}, 
    blk2::Array{CT, 3}, 
    es::Bool) where {CT<:Number}
    # Contract the blocks blk1 and blk2
    if es
        @tensor result[i1, i2, i3, ν, μ] := blk1[i1, i2, e, μ] * blk2[e, i3, ν]
    else
        @tensor result[i1, i2, i3, λ, κ] := blk1[i2, i3, f, κ] * blk2[i1, f, λ]
    end
    return result
end

# return type: a dictionary
# Its keys are the qlabels, and its values are tuples of multiplicities
# values are (a, b), where a is the multiplicity in in1 ⊗ in2 -> e
# b is the multiplicity in e ⊗ in3 -> out
function find_espaces(::Type{S},
    in1::NTuple{NZ, Int},
    in2::NTuple{NZ, Int},
    in3::NTuple{NZ, Int},
    out::NTuple{NZ, Int},
    ::Type{CT};
    verbose=0) where {S<:NonabelianSymm, CT<:Number, NZ}
    es = Vector{Tuple{NTuple{NZ, Int}, Tuple{Int, Int}}}()
    in12_out = getNsave_validout(S, minmax(in1, in2))
    for qlabel in in12_out.out_spaces
        qin3_omlist = getNsave_omlist(S, minmax(in3, qlabel), out)
        if qin3_omlist === nothing continue end

        in12_omlist = getNsave_omlist(S, minmax(in1, in2), qlabel)
        in12_om = in12_omlist.totalOM
        qin3_om = qin3_omlist.totalOM
        if verbose > 0
            println("Found espace: $qlabel with multiplicity ($(in12_om), $(qin3_om))")
        end
        push!(es, (qlabel, (qin3_om, in12_om)))
    end
    return es
end

function fsym_leftnormalized(fsym::Fsymbol{S, CT, NZ}) where {S<:NonabelianSymm, CT<:Number, NZ}
    @assert nzops(S) == NZ
    mat_rational = Diagonal(fsym.es_nfactor) * fsym.fsym_mat
    nfac::Integer = lcm(denominator.(mat_rational))
    return Matrix{BigInt}(mat_rational * nfac), nfac
end

function fsym_rightnormalized(fsym::Fsymbol{S, CT, NZ}) where {S<:NonabelianSymm, CT<:Number, NZ}
    @assert nzops(S) == NZ
    mat_rational = fsym.fsym_mat * Diagonal(fsym.fs_nfactor)
    nfac::Integer = lcm(denominator.(mat_rational))
    return Matrix{BigInt}(mat_rational * nfac), nfac
end

Base.size(fsym::Fsymbol{S, CT, NZ}) where {S<:NonabelianSymm, CT, NZ} = size(fsym.fsym_mat)[1]