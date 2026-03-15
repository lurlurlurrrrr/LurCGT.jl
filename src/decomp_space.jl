comm(A, B) = A * B - B * A

# i, j: To assert that weights are valid
function toblk(symm::NTuple{N, Any},
    mat::SparseArray{Int, 2},
    weight_zips::Vector{NTuple{N, Tuple{Vararg{Int}}}},
    omvec::Dict{NTuple{N, Tuple{Vararg{Int}}}, Vector{Int}},
    nthweight::Vector{Int},
    i::Int,
    j::Int) where N
    blks = Dict{NTuple{N, Tuple{Vararg{Int}}}, SparseMatrixCSC{Int}}()
    S = symm[i]
    for (k, v) in mat.data
        a, b = k.I
        wa, na = weight_zips[a], nthweight[a]
        wb, nb = weight_zips[b], nthweight[b]
        for ii in 1:N if ii != i @assert wa[ii] == wb[ii] end end
        if S<:NonabelianSymm @assert wa[i] == wb[i] .- getdz(S, j) end

        ma, mb = omvec[wa], omvec[wb]
        ima, imb = length(ma), length(mb)
        if !haskey(blks, wb) blks[wb] = spzeros(Int, ima, imb) end

        blks[wb][na, nb] = v
    end
    return blks
end

function add_vectors!(ortho_vecs_sparse::SparseMatrixCSC{Float64},
    mult_inds::Vector{Int},
    vectors::Array{BigInt, M},
    stidx::Int,
    inc::Int,
    reps::NTuple{N, Irep},
    w::NTuple{N, Tuple{Vararg{Int}}}) where {N, M}

    @assert M == N + 1
    rep = reps[end]
    a, b = rep.Sz[w[end]]; @assert size(vectors, M) == b - a + 1
    inc *= dimension(rep)
    inc += a - 1
    if N == 1 
        nvec = size(vectors, 2)
        si = stidx + inc
        ortho_vecs_sparse[mult_inds, si:si+nvec-1] = vectors
        return
    end

    for i in 1:size(vectors, M)
        add_vectors!(ortho_vecs_sparse, mult_inds, vectors[(Colon() for _=1:N)..., i], 
        stidx, inc+i-1, reps[1:end-1], w[1:end-1])
    end
end

function add_irops!(irop_3d::SparseArray{Float64, 3},
    irops::Array{SparseMatrixCSC{Int}, N},
    si::Int,
    reps::NTuple{N, Irep},
    w::NTuple{N, Tuple{Vararg{Int}}}) where N

    rep = reps[end]
    a, b = rep.Sz[w[end]]; @assert size(irops, N) == b - a + 1
    si *= dimension(rep)
    si += a - 1
    if N == 1
        nvec = size(irops, 1)
        for i in 1:nvec irop_3d[:, :, si+i] = irops[i] end
        return
    end

    for i in 1:size(irops, N)
        add_irops!(irop_3d, irops[(Colon() for _=1:N-1)..., i], 
        si, reps[1:end-1], w[1:end-1])
    end
end

# Assume that local Hilbert space is orthonormal
# symm: tuple of symmetries 
function decompose_space(symm::NTuple{N, Any},
    weights::NTuple{N, Vector{<:Tuple{Vararg{Int}}}},
    lops::NTuple{N, Vector{<:AbstractMatrix{Int}}}) where N

    # Convert all lowering operators to sparse matrices
    lops_sparse = ntuple(i->SparseArray.(lops[i]), N)
    weights_zip = Vector{NTuple{N, Tuple{Vararg{Int}}}}()
    spdim = length(weights[1])
    for i in 1:spdim
        push!(weights_zip, ntuple(j->weights[j][i], N))
    end

    # Key: weights, value: vector of indices of states with the weights.
    mult_vecs = Dict{NTuple{N, Tuple{Vararg{Int}}}, Vector{Int}}()
    nthweight = Vector{Int}()

    sorted_weights = sort(weights_zip; lt=rev_less_symms, rev=true)

    for i in 1:spdim
        w = weights_zip[i]
        if haskey(mult_vecs, w) push!(mult_vecs[w], i)
        else mult_vecs[w] = [i] end
        push!(nthweight, length(mult_vecs[w]))
    end
    f = (lop, i, j) -> toblk(symm, lop, weights_zip, mult_vecs, nthweight, i, j)
    lops_blk = Tuple([f(lops_sparse[i][j], i, j) for j=1:length(lops_sparse[i])] for i=1:N)

    # Overall structure of this algorithm is similar to the one in clebsch.jl
    # Inner multiplicity count
    remain_inmult = Dict{NTuple{N, Tuple{Vararg{Int}}}, Int}()
    for (k, v) in mult_vecs remain_inmult[k] = length(v) end

    # foundvecs is different to the one in clebsch.jl
    # First key: qlabels of symmetry sectors
    # Second key: weight tuple, value: array of obtained vectors
    # Obtained vectors are not orthogonalized
    # Dim of array: (spdim1, N inner_mults..., outer_mult)
    foundvecs = Vector{Tuple{NTuple{N, Tuple{Vararg{Int}}}, 
        Vector{Tuple{NTuple{N, Tuple{Vararg{Int}}}, Array{BigInt, N+2}}}}}()

    # Key: qlabels tuple, value: vector of norms fo the mwstates
    mwstate_norms = Dict{NTuple{N, Tuple{Vararg{Int}}}, Vector{BigInt}}()
    remain_dim = spdim
    

    # Starting from the highest weight state
    while remain_dim > 0
        for w in sorted_weights
            outer_mult = remain_inmult[w]
            if outer_mult > 0
                qlabels_out = ntuple(i->getqlabel(symm[i], w[i]), N)
                output_ireps = ntuple(i->getNsave_irep(symm[i], BigInt, qlabels_out[i]), N)

                # println("Detected new sector: qlabels=$(qlabels_out), outer_mult=$(outer_mult)")

                # Find a new highest weight state vector from 'foundvecs'
                inmult = length(mult_vecs[w])
                already_found = Matrix{BigInt}(undef, inmult, inmult-remain_inmult[w])
                already_found_norms = BigInt[]
                stidx = 1
                # Fill in already_found and already_found_norms

                # TODO: Modify this part since I no longer use Tmat & Nvec.
#                for (kq, dict) in foundvecs
#                    if haskey(dict, w)
#                        vecs = get_orthogonal_vecs(symm, kq, w, dict[w])
#                        nvecs = prod(size(vecs)[2:end])
#                        @assert nvecs == size(vecs, 2)
#                        
#                        vecs_norms = get_vecs_norms(symm, kq, w, mwstate_norms[kq])
#                        append!(already_found_norms, vecs_norms)
#
#                        already_found[:, stidx:stidx+nvecs-1] = vecs
#                        stidx += nvecs
#                    end
#                end
#                @assert stidx == inmult - remain_inmult[w] + 1
#                for i in 1:size(already_found, 2)
#                    @assert dot(already_found[:, i], already_found[:, i]) == already_found_norms[i]
#                end

                mwstate_norms_sector = Vector{BigInt}(undef, outer_mult)

                mwstates = Matrix{BigInt}(undef, inmult, outer_mult)
                start_index = 1
                for oidx=1:outer_mult
                    mwstates[:, oidx], start_index, mwnorm = 
                        get_mwstate(already_found, already_found_norms, start_index)
                    mwstate_norms_sector[oidx] = mwnorm
                end

                multiplet = get_vectors_sector(output_ireps, lops_blk, mwstates, w)
                vecs_w = Vector{Tuple{NTuple{N, Tuple{Vararg{Int}}}, Array{BigInt, N+2}}}()
                for w in sort(collect(keys(multiplet)); lt=lt_multisymm_weights, rev=true)
                    nvec = prod(size(multiplet[w])[2:end])
                    remain_inmult[w] -= nvec
                    remain_dim -= nvec
                    push!(vecs_w, (w, multiplet[w]))
                end
                push!(foundvecs, (qlabels_out, vecs_w))
                mwstate_norms[qlabels_out] = mwstate_norms_sector
                # println("Found new sector: qlabels=$(qlabels_out), outer_mult=$(outer_mult)")
            end
        end
    end
    ortho_vecs_sparse = sparse(zeros(Float64, spdim, spdim))

    mult_ind = Dict{NTuple{N, Tuple{Vararg{Int}}}, Vector{Tuple{Int, Int}}}()
    oldq = nothing
    stidx = 1
    for (q, vecs_w) in foundvecs
        # println("Multiplet qlabels: ", q)
        reps = Tuple(getNsave_irep(symm[i], BigInt, q[i]) for i=1:N)
        totaldim = prod(dimension(reps[i]) for i=1:N)
        om = size(vecs_w[1][2], N+2)
        for (w, multiplet) in vecs_w
            # println("Weights: ", w)

            w_mult = size(multiplet, 1)
            if !haskey(mult_ind, q) || isempty(mult_ind[q])
                mult_ind[q] = Vector{Tuple{Int, Int}}()
                for i in 1:om
                    inc = (i - 1) * totaldim
                    push!(mult_ind[q], (stidx + inc, stidx + inc + totaldim - 1))
                end
            end

            mult_inds = mult_vecs[w]; @assert length(mult_inds) == w_mult
            for i in 1:om
                inc = (i - 1) * totaldim
                ith_slice = multiplet[(Colon() for _=1:N+1)..., i]
                add_vectors!(ortho_vecs_sparse, mult_inds, ith_slice, stidx+inc, 0, reps, w)
            end
            oldq = q
        end
        stidx += om * totaldim
        # println(stidx)
    end
        
    for (q, vec) in mult_ind
        # println(q, vec)
        outputmats = [cgt_outputmat(getNsave_irep(symm[i], BigInt, q[i]), Float64) for i=1:N]
        #for mat in outputmats display(mat) end
        outputmat = reduce(⊗, outputmats)
        #display(outputmat)
        for (sti, edi) in vec
            ortho_vecs_sparse[:, sti:edi] = ortho_vecs_sparse[:, sti:edi] * outputmat
        end
    end

    @assert ortho_vecs_sparse * ortho_vecs_sparse' ≈ I
    return ortho_vecs_sparse, mult_ind
end

# There is a similar version of this function when getting IROP
function get_vectors_sector(ireps::NTuple{N, Irep},
    lops_blk::NTuple{N, Vector{Dict{NTuple{N, Tuple{Vararg{Int}}}, SparseMatrixCSC{Int}}}},
    mwstates::Matrix{BigInt},
    mw::NTuple{N, Tuple{Vararg{Int}}}) where N

    vecs_sector = Dict{NTuple{N, Tuple{Vararg{Int}}}, Array{BigInt, N+2}}()
    mw_inmult, om = size(mwstates)
    vecs_sector[mw] = reshape(mwstates, mw_inmult, ntuple(_->1, N)..., om)

    get_vectors_sector_!(Val(1), vecs_sector, lops_blk, ireps, mw)
    return vecs_sector
end

# I: start from 1, increased by 1 for every recursion step.
function get_vectors_sector_!(::Val{I},
    vecs_sector::Dict{NTuple{N, Tuple{Vararg{Int}}}, Array{BigInt, M}},
    lops_blk::NTuple{N, Vector{Dict{NTuple{N, Tuple{Vararg{Int}}}, SparseMatrixCSC{Int}}}},
    ireps::NTuple{N, Irep},
    mw::NTuple{N, Tuple{Vararg{Int}}}) where {I, N, M}

    @assert M == N + 2
    lops = lops_blk[I]; irep = ireps[I]
    S, qlabel = symm(irep), irep.qlabel
    @assert mw[I] == qlab2mwz(S, qlabel)

    if S<:NonabelianSymm
        dzs = getdzs(S)
        sorted_weights = sort(collect(keys(irep.Sz)); lt=rev_less, rev=true)
        for (i, w) in enumerate(sorted_weights)
            if i == 1 continue end

            m, n = irep.Sz[w]; w_mult = n - m + 1
            remaining_set = Set(1:w_mult)
            
            for opidx in 1:nlops(S)
                init_weight = Tuple(w .+ dzs[opidx])
                if haskey(irep.Sl[opidx], init_weight)
                    lop_irep = irep.Sl[opidx][init_weight]
                    for j in axes(lop_irep, 2)
                        jth_colm = lop_irep[:, j]
                        if nnz(jth_colm) == 1
                            nzind, nzval = jth_colm.nzind[1], jth_colm.nzval[1]
                            if nzind in remaining_set
                                wt = tup_change_ith(Val(I), mw, w)
                                wtup = tup_change_ith(Val(I), mw, init_weight)
                                lop_blk = lops[opidx][wtup]
                                # If vecs_sector does not have wtup, create a zero array
                                if !haskey(vecs_sector, wt)
                                    sz1 = size(lop_blk, 1)
                                    sz_wtup = size(vecs_sector[wtup])
                                    new_size = [sz1, sz_wtup[2:I]..., w_mult, sz_wtup[I+2:end]...]
                                    for i=I+2:M @assert new_size[i] == 1 end
                                    vecs_sector[wt] = zeros(BigInt, new_size...)
                                end
                                # Fill vecs_sector
                                init_idx = [[Colon() for _=1:I]..., j, [1 for _=I+1:M]...]
                                final_idx = [[Colon() for _=1:I]..., nzind, [1 for _=I+1:M]...]
                                vecs_sector[wt][final_idx...] = div.(reshapeNprod(
                                    vecs_sector[wtup][init_idx...], lop_blk), nzval)
                                delete!(remaining_set, nzind)
                            end
                        end
                    end
                end
            end
            @assert isempty(remaining_set)
        end
    end

    if I < N 
        # For every weights, call get_vectors_sector_! again
        for w in keys(irep.Sz)
            mw_ = tup_change_ith(Val(I), mw, w)
            get_vectors_sector_!(Val(I+1), vecs_sector, lops_blk, ireps, mw_)
        end
    end 
end

function reshapeNprod(arr::Array{T, N}, mat::AbstractMatrix{S}) where {T, S, N}
    @assert N >= 2
    res_mat = mat * reshape(arr, size(arr, 1), :)
    return reshape(res_mat, size(mat, 1), size(arr)[2:end]...)
end

tup_change_ith(::Val{I}, tup::NTuple{N, Tuple{Vararg{Int}}}, e::Tuple{Vararg{Int}}) where {N, I} = 
    tuple(tup[1:I-1]..., e, tup[I+1:end]...)

function get_mwstate(already_found::Matrix{BigInt},
    already_found_norms::Vector{BigInt},
    start_index::Int)

    inmult = size(already_found, 1)
    nvec_found = size(already_found, 2)
    @assert length(already_found_norms) == nvec_found
    mwstate = zeros(BigInt, inmult)
    found = false
    while !found
        innerprod_fac = 1
        mwstate = zeros(BigInt, inmult); mwstate[start_index] = 1
        for i in 1:nvec_found
            overlap = already_found[start_index, i] * innerprod_fac
            found_vec = already_found[:, i]
            multfac, divfac, mwstate = orthogonalize(mwstate, 
                found_vec, overlap, already_found_norms[i])
            innerprod_fac *= multfac
            facs_gcd = gcd(innerprod_fac, divfac)
            mwstate *= div(divfac, facs_gcd)
            innerprod_fac = div(innerprod_fac, facs_gcd)
            if iszero(mwstate) break end
        end
        if iszero(mwstate) start_index += 1; continue
        else found = true end
    end
    _, mwstate = divcfac(mwstate)
    mwstate_norm = dot(mwstate, mwstate)
    first_nzidx = findfirst(!iszero, mwstate)
    if mwstate[first_nzidx] < 0 mwstate = -mwstate end
    return mwstate, start_index + 1, mwstate_norm
end

function get_orthogonal_vecs(symm::NTuple{N, Any},
    qlabels::NTuple{N, Tuple{Vararg{Int}}},
    weight::NTuple{N, Tuple{Vararg{Int}}},
    vecs::Array{BigInt, M}) where {N, M}

    @assert M == N + 2
    for i in 1:N
        irep = getNsave_irep(symm[i], BigInt, qlabels[i])
        Tmat = irep.Tmats[weight[i]]
        vecs = contract_ith(vecs, Tmat, Val(i+1))
    end
    return reshape(vecs, size(vecs, 1), :)
end

function get_vecs_norms(symm::NTuple{N, Any},
    qlabels::NTuple{N, Tuple{Vararg{Int}}},
    weight::NTuple{N, Tuple{Vararg{Int}}},
    mw_norms::Vector{BigInt}) where N

    norms_vec = Vector{Vector{BigInt}}(undef, N+1)
    norms_vec[N+1] = mw_norms
    for i in 1:N
        irep = getNsave_irep(symm[i], BigInt, qlabels[i])
        norms_vec[i] = irep.Nvecs[weight[i]]
    end
    return reduce(⊗, norms_vec)
end


# Test input. N-channel spinful fermionic system
# TODO: Test whether operators from different symmetry commute
function test_input(nchannels::Int=3)
    N = nchannels
    @assert N >= 2 "Need at least 2 channels"
    
    ss = f4down' * f4up
    # Annihilation operations, dim: (channel, spin)
    # FF[1, n]: spin-up annihilation operator for channel n
    # FF[2, n]: spin-down annihilation operator for channel n
    FF = Matrix{SparseMatrixCSC{Int}}(undef, 2, N)
    NN = Matrix{SparseMatrixCSC{Int}}(undef, 2, N)
    # SS[n]: spin lowering operator for channel n
    SS = Vector{SparseMatrixCSC{Int}}(undef, N)
    
    # Total number of fermionic modes: 2 spins × N channels
    total_modes = 2 * N
    
    for i in 1:total_modes
        mats = vcat([z2 for _=1:i-1], [f2], [I2 for _=i+1:total_modes])
        FF[i] = reduce(⊗, mats)
        NN[i] = FF[i]' * FF[i]
    end
    
    for i in 1:N
        mats = vcat([I2 for _=1:2*(i-1)], [ss], [I2 for _=1:2*(N-i)])
        SS[i] = reduce(⊗, mats)
        @assert SS[i] == FF[2, i]' * FF[1, i]
    end

    SZ = Vector{SparseMatrixCSC{Int}}(undef, N)
    for i in 1:N
        mats = vcat([I2 for _=1:2*(i-1)], [sz4], [I2 for _=1:2*(N-i)])
        SZ[i] = reduce(⊗, mats)
    end
    
    # Charge relative to half-filling
    charge = Int.(diag(sum(NN)) .- N)
    charge_t = [(Int(i),) for i in charge]
    spin_z = Int.(diag(sum(SZ)))
    spin_z_t = [(Int(i),) for i in spin_z]
    spin_lowering = sum(SS)

    # Channel lowering operators for SU{N}: N-1 operators
    # chan_l[i] connects channel i to channel i+1
    chan_l = Vector{SparseMatrixCSC{Int}}(undef, N-1)
    for i in 1:(N-1)
        chan_l[i] = FF[1, i+1]'*FF[1, i] + FF[2, i+1]'*FF[2, i]
    end
    
    # Channel z-operators for SU{N}: N-1 operators following Cartan subalgebra
    # Using standard Dynkin basis: H_i = E_{i,i} - E_{i+1,i+1}
    chan_z = Vector{Vector{Int}}(undef, N-1)
    for i in 1:(N-1)
        # H_i: occupation of channel i minus occupation of channel i+1
        chan_z[i] = Int.(diag(sum([NN[1, j] + NN[2, j] for j=1:i]) 
        - i * NN[1, i+1] - i * NN[2, i+1]))
    end
    
    # Verify commutation relations for SU{N} algebra
    @assert comm(spin_lowering, diagm(spin_z)) == 2 * spin_lowering
    for i in 1:(N-1)
        # [H_i, E_i] = 2*E_i (lowering operator i raised by H_i by 2)
        @assert comm(chan_l[i], diagm(chan_z[i])) == (i+1) * chan_l[i]
        
        # [H_i, E_{i-1}] = -E_{i-1} (lowering operator i-1 lowered by H_i by 1)
        if i > 1
            @assert comm(chan_l[i-1], diagm(chan_z[i])) == 0 * chan_l[i-1]
            @assert comm(chan_l[i], diagm(chan_z[i-1])) == -(i-1) * chan_l[i]
        end

        
        # [H_i, E_j] = 0 for |i-j| > 1
        for j in 1:(N-1)
            if abs(i - j) > 1
                @assert comm(chan_l[j], diagm(chan_z[i])) == 0 * chan_l[j]
            end
        end
    end

    symms = (U1, SU{2}, SU{N})
    chan_weights = collect(zip(chan_z...))
    weights = (charge_t, spin_z_t, chan_weights)
    lops = (Matrix{Int}[], [spin_lowering], chan_l)
    # 1. Decompose space and save orthonormal basis into sparse matrix 
    vecs, mult_ind = decompose_space(symms, weights, lops)

    # 2. Get IROP and get matrix elements of operators in the decomposed space
    totalN = diag(sum(NN[i] for i in 1:2*N))
    mwirop_S = sum(SS[i]' for i in 1:N)
    mwirop_F = FF[2, N] # Spin down annihilation operator for the last channel
    mwirop_Z = diagm([i%2==0 ? 1 : -1 for i in totalN])
    mwirop_I = sparse(I, 4^N, 4^N)

    # They are stored in 3D arrays with dim: (local dim, local dim, # of irop ops)
    # println("Getting IROP for spin lowering operator...")
    irop_S, qlabel_S = get_IROP(symms, weights, lops, mwirop_S)
    # println("Getting IROP for fermionic annihilation operator...")
    irop_F, qlabel_F = get_IROP(symms, weights, lops, mwirop_F)
    # println("Getting IROP for fermionic sign operator...")
    irop_Z, qlabel_Z = get_IROP(symms, weights, lops, mwirop_Z)
    # println("Getting IROP for identity operator...")
    irop_I, qlabel_I = get_IROP(symms, weights, lops, mwirop_I)

    # 3. Divide it into CGC and RMT to get final QSpace object
    transf_basis!(irop_S, vecs)
    transf_basis!(irop_F, vecs)
    transf_basis!(irop_Z, vecs)
    transf_basis!(irop_I, vecs)

    # println("Decomposing IROP for spin lowering operator...")
    data_S = decompose_irop(symms, irop_S, mult_ind, qlabel_S)
    # println("Decomposing IROP for fermionic annihilation operator...")
    data_F = decompose_irop(symms, irop_F, mult_ind, qlabel_F)
    # println("Decomposing IROP for fermionic sign operator...")
    data_Z = decompose_irop(symms, irop_Z, mult_ind, qlabel_Z)
    # println("Decomposing IROP for identity operator...")
    data_I = decompose_irop(symms, irop_I, mult_ind, qlabel_I)

    return data_S, data_F, data_Z, data_I
end

# N: The number of symmetries
function decompose_irop(symm::NTuple{N, Any},
    irop::SparseArray{Float64, 3},
    mult_ind::Dict{NTuple{N, Tuple{Vararg{Int}}}, Vector{Tuple{Int, Int}}},
    qlabel::NTuple{N, Tuple{Vararg{Int}}}) where N

    # println(mult_ind)
    CGT_om = 0
    data = Vector{Tuple{NTuple{3, NTuple{N, Tuple{Vararg{Int}}}}, Array{Float64, N+3}}}()
    sorted_qs = sort(collect(keys(mult_ind)); lt=rev_less_symms)
    for k1 in sorted_qs
        om1 = length(mult_ind[k1])
        range1 = mult_ind[k1]
        for k2 in sorted_qs
            RMT, CGT = nothing, nothing
            om2 = length(mult_ind[k2])
            range2 = mult_ind[k2]

            for i1 in 1:om1
                a1, b1 = range1[i1]
                for i2 in 1:om2
                    a2, b2 = range2[i2]
                    block = irop[a1:b1, a2:b2, :]
                    if !iszero(block) 
                        # println("$(i1)th $k1, $(i2)th $k2, $qlabel, norm: $(norm(block))") 
                        if isnothing(RMT) 
                            @assert isnothing(CGT)
                            CGTs = [load_cg3_float(symm[i], BigInt, (k2[i], qlabel[i]), 
                            k1[i], Float64, true)[1] for i=1:N]
                            CGT = permutedims(reduce(⊗, CGTs), (3, 1, 2, 4))
                            # Normalize CGT
                            CGT /= sqrt(size(CGT, 1))
                            CGT_oms = [size(CGTs[i], 4) for i in 1:N]
                            CGT_om = size(CGT, 4); @assert CGT_om == prod(CGT_oms)
                            @assert norm(CGT) ≈ sqrt(CGT_om)
                            RMT = zeros(Float64, om1, om2, 1, CGT_oms...) 
                        end
                        RMT[i1, i2, 1, (Colon() for _=1:N)...] = 
                        reshape(reshape(CGT, :, CGT_om)' * block[:], CGT_oms...)

                        # Check by reconstructing the block from CGT and RMT
                        @assert block ≈ reshape(reshape(CGT, :, CGT_om) * 
                        RMT[i1, i2, 1, (Colon() for _=1:N)...][:], size(block))
                    end
                end
            end
            if !isnothing(RMT) push!(data, ((k1, k2, qlabel), RMT)) end
        end
    end
    return data
end

function transf_basis!(irop::SparseArray{Float64, 3}, 
    vecs::SparseMatrixCSC{Float64})
    for i in 1:size(irop, 3)
        irop[:, :, i] = vecs' * irop[:, :, i] * vecs
    end
    droptol!(irop, 1e-12)
end

function droptol!(sparr::SparseArray{Float64, N}, tol::Float64) where N
    keys_to_drop = [k for (k, v) in nonzero_pairs(sparr) if abs(v) < tol]
    for k in keys_to_drop sparr[k] = 0.0 end
end


function getweight_irop(::Type{S},
    z_ops::Vector{<:Tuple{Vararg{Int}}},
    lowering_ops::Vector{<:AbstractMatrix{Int}},
    mwirop::AbstractMatrix{<:Integer}) where S<:Symmetry

    # Get the qlabel of the maximal weight operator in the multiplet that mwirop belongs to
    NZ = nzops(S); @assert NZ == length(z_ops[1])
    nzidx = findfirst(!=(0), mwirop)
    zops_diags = [[z_ops[i][j] for i=1:length(z_ops)] for j=1:NZ]
    zops = [diagm(zops_diags[j]) for j=1:NZ]
    comm_result = [comm(zops[j], mwirop) for j=1:NZ]
    zvals = Tuple(div(comm_result[i][nzidx], mwirop[nzidx]) for i=1:NZ)
    for i=1:NZ @assert comm_result[i] == zvals[i] * mwirop end
    return zvals
end

function get_irops_sector(ireps::NTuple{N, Irep},
    lowering_ops::NTuple{N, Vector{<:AbstractMatrix{Int}}},
    mwirop::AbstractMatrix{<:Integer},
    mw::NTuple{N, Tuple{Vararg{Int}}}) where N

    irops_sector = Dict{NTuple{N, Tuple{Vararg{Int}}}, Array{SparseMatrixCSC{Int}, N}}()
    irops_sector[mw] = reshape([mwirop], ntuple((_->1, N)...))

    get_irops_sector_!(Val(1), irops_sector, lowering_ops, ireps, mw)
    return irops_sector
end

function get_irops_sector_!(::Val{I},
    irops_sector::Dict{NTuple{N, Tuple{Vararg{Int}}}, Array{SparseMatrixCSC{Int}, N}},
    lowering_ops::NTuple{N, Vector{<:AbstractMatrix{Int}}},
    ireps::NTuple{N, Irep},
    mw::NTuple{N, Tuple{Vararg{Int}}}) where {I, N}

    lops = lowering_ops[I]; irep = ireps[I]
    S, qlabel = symm(irep), irep.qlabel
    @assert mw[I] == qlab2mwz(S, qlabel)

    if S<:NonabelianSymm
        dzs = getdzs(S)
        sorted_weights = sort(collect(keys(irep.Sz)); lt=rev_less, rev=true)
        for (i, w) in enumerate(sorted_weights)
            if i == 1 continue end

            m, n = irep.Sz[w]; w_mult = n - m + 1
            remaining_set = Set(1:w_mult)

            for opidx in 1:nlops(S)
                lop = lops[opidx]
                init_weight = Tuple(w .+ dzs[opidx])
                if haskey(irep.Sl[opidx], init_weight)
                    lop_irep = irep.Sl[opidx][init_weight]
                    for j in axes(lop_irep, 1)
                        jth_colm = lop_irep[:, j]
                        if nnz(jth_colm) == 1
                            nzind, nzval = jth_colm.nzind[1], jth_colm.nzval[1]
                            if nzind in remaining_set
                                wt = tup_change_ith(Val(I), mw, w)
                                wtup = tup_change_ith(Val(I), mw, init_weight)

                                if !haskey(irops_sector, wt)
                                    sz_wtup = size(irops_sector[wtup])
                                    new_size = [sz_wtup[1:I-1]..., w_mult, sz_wtup[I+1:end]...]
                                    for i=I+1:N @assert new_size[i] == 1 end
                                    irops_sector[wt] = Array{SparseMatrixCSC{Int}, N}(undef, new_size...)
                                end
                                slice = [[Colon() for _=1:I-1]..., j, [1 for _=I+1:N]...]
                                for idx in CartesianIndices(irops_sector[wtup][slice...])
                                    idx_init = [idx.I..., j, [1 for _=I+1:N]...]
                                    idx_final = [idx.I..., nzind, [1 for _=I+1:N]...]
                                    comm_res = comm(lop, irops_sector[wtup][idx_init...])
                                    irops_sector[wt][idx_final...] = div.(comm_res, nzval)
                                end
                                delete!(remaining_set, nzind)
                            end
                        end
                    end
                end
            end
            @assert isempty(remaining_set)
        end
    end

    if I < N
        for w in keys(irep.Sz)
            mw_ = tup_change_ith(Val(I), mw, w)
            get_irops_sector_!(Val(I+1), irops_sector, lowering_ops, ireps, mw_)
        end
    end
end

function get_IROP(symm::NTuple{N, Any},
    z_ops::NTuple{N, Vector{<:Tuple{Vararg{Int}}}},
    lowering_ops::NTuple{N, Vector{<:AbstractMatrix{Int}}},
    mwirop::AbstractMatrix{<:Integer}) where N

    spdim = length(z_ops[1])
    # First, get the weight of maximal weight operator
    # Then, apply the lowering operators to get the weights of all operators in the multiplet
    mw = Tuple(getweight_irop(symm[i], z_ops[i], lowering_ops[i], mwirop) for i in 1:N)
    qlabels_out = Tuple(getqlabel(symm[i], mw[i]) for i in 1:N)
    output_ireps = ntuple(i->getNsave_irep(symm[i], BigInt, qlabels_out[i]), N)
    norm2_mwirop = dot(mwirop, mwirop)
    nirop = prod(dimension(output_ireps[i]) for i=1:N)
        
    # println(qlabels_out)
    # println("Norm^2 of the maximal weight operator: ", norm2_mwirop)
    # Get other operators in IROP, don't think outer multiplicity of IROP for now
    irops = get_irops_sector(output_ireps, lowering_ops, mwirop, mw)

    # Orthogonalize and normalize so that every operator has norm equal to mwirop
    # Check orthogonality and proper normalization
    irops_3d = SparseArray(zeros(Float64, spdim, spdim, nirop))

    # First, fill in irops_3d with the obtained irops
    for (w, irop_arr) in irops
        add_irops!(irops_3d, irop_arr, 0, output_ireps, w)
    end

    # Second, orthogonalize and normalize the irops, and fill in irops_3d_ortho
    outputmats = [cgt_outputmat(getNsave_irep(symm[i], BigInt, qlabels_out[i]), Float64) for i=1:N]
    outputmat = SparseArray(reduce(⊗, outputmats))

    multiplied = reshape(irops_3d, spdim^2, nirop) * outputmat
    irops_3d_ortho = reshape(multiplied, spdim, spdim, nirop)

    for i in 1:nirop
        slice_i = irops_3d_ortho[:, :, i]
        for j in i:nirop
            slice_j = irops_3d_ortho[:, :, j]
            @assert dot(slice_i, slice_j) ≈ (i==j ? norm2_mwirop : 0)
        end
    end

    # Return 3d array of dim (local dim, local dim, # of irop ops)
    return irops_3d_ortho, qlabels_out
end
