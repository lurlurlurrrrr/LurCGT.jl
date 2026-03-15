const WeightInfo{NZ} = Dict{Tuple{NTuple{NZ, Int}, NTuple{NZ, Int}}, Tuple{Int, Int}}
const CG3block{NZ, CT} = Dict{NTuple{3, NTuple{NZ, Int}}, Array{CT, 4}}

# S : symmetry type
# RT : numeric type used in lowering operator and norm2s of representations
# CT : numeric type used in cg coefficients
# Just save the trivial and defining irreducible representations,
# and possible alternating product representations of defining rep.
function get_fundamental_ireps(::Type{S}, ::Type{RT}; 
        verbose=0) where {S<:NonabelianSymm, RT<:Number}
	@assert isvalidsymm(S)
    println("Start generating fundamental irreps of symmetry $(S)...")
    # Trivial (1-dimensional) representation
    trivirep = gettrivirep(S, RT)
    save_irep_sqlite(trivirep)
    
	defirep = getdefirep(S, RT)
    save_irep_sqlite(defirep)

    # NZ = number of lowering operators
    NZ = nzops(S)

    # Get alternating product representations
    for n=2:NZ
        Sl_alt, Sz_alt, mw, inner_prod = def_altprod(defirep, n)

        # q-label of the new representation
        # Tuple with length NZ, all elements are zero except the nth element
        irep_alt = get_irep_frommw(S, Sl_alt, Sz_alt, inner_prod, 
        nothing, mw, nothing; verbose)
        save_irep_sqlite(irep_alt)
        if onemore_fundamental_irep(S, n)
            # SO(2N) case when n == NZ
            mw = mw .- Tuple(i==1 ? 4 : 0 for i=1:NZ)
            irep_alt2 = get_irep_frommw(S, Sl_alt, Sz_alt, inner_prod, 
            nothing, mw, nothing; verbose)
            save_irep_sqlite(irep_alt2)
        end
    end
end

onemore_fundamental_irep(::Type{S}, n::Int) where S<:NonabelianSymm = false
onemore_fundamental_irep(::Type{SO{N}}, n::Int) where N = 
(N % 2 == 0 && n == div(N, 2))

# qlabels_in = (q1, q2), qlabels_out = [q3, q4, ...]
# Then, generate CG3s (q1⊗q2->q3), (q1⊗q2->q4), ...
# They are returned as a Dict{qlabel, CG3}
function generate_every_CGT(::Type{S}, ::Type{RT}, ::Type{CT}, 
    qlabels_in::NTuple{2, NTuple{NZ, Int}},
    qlabels_out::Union{Vector{NTuple{NZ, Int}}, Nothing};
    verbose=0,
    assertlev=0,
    save=true) where {S<:NonabelianSymm, RT<:Number, CT<:Number, NZ}

    q1, q2 = qlabels_in
    @assert nzops(S) == NZ
    @assert q1 <= q2
    irep1 = getNsave_irep(S, RT, q1); dim1 = dimension(irep1)
    irep2 = getNsave_irep(S, RT, q2); dim2 = dimension(irep2)
    # remain_inmult : remaining multiplicity of the product representation
    # values will be decreases as CGTs are created
    sl_prod, sz_prod, prod_inmult, prod_winfo = 
    tensor_prod_reps(irep1, irep2; verbose)

    # Basic setting to see which weights should be processed
    vo = getNsave_validout(S, qlabels_in); @assert !isnothing(vo)
    possib_out = vo.out_spaces # This is sorted
    possib_out_sorted = sort(possib_out; by=q->reverse(qlab2mwz(S, q)))

    qlabels_out::Vector{NTuple{NZ, Int}} = 
    sort(isnothing(qlabels_out) ? possib_out : qlabels_out; by=q->reverse(qlab2mwz(S, q)))

    q_smallest = qlabels_out[1]

    pout_idx = findfirst(==(q_smallest), possib_out_sorted)
    pout_set = Set{NTuple{NZ, Int}}(possib_out_sorted[pout_idx:end])
    qout_set = Set{NTuple{NZ, Int}}(qlabels_out)
    # The weights need to be processed
    wout_set = Set{NTuple{NZ, Int}}(qlab2mwz(S, q) for q in pout_set)
    for q in qout_set @assert q ∈ pout_set end

    # Key: q-label of output representation, Value: CG3 object
    result_cg3s = Dict{NTuple{NZ, Int}, CGT{S, CT, NZ, 3}}()

    remain_inmult = copy(prod_inmult)
    # Vector of (qlabel, outer multiplicity)
    result_qlabels = Vector{Tuple{NTuple{NZ, Int}, Int}}()
    # Innerprod matrices of the output representations
    innerprod_dict = Dict{NTuple{NZ, Int}, Dict{NTuple{NZ, Int}, Matrix{RT}}}()
    # Inverse of innerprod matrices of the output representations
    inv_innerprod_dict = Dict{NTuple{NZ, Int}, Dict{NTuple{NZ, Int}, Tuple{Matrix{RT}, Rational{RT}}}}()

    # (Possibly non-orthogonalized) already-found states in the product representation
    foundvecs = Dict{NTuple{NZ, Int}, Matrix{CT}}()
    # Norms of found maximal weight states, key is q-label
    mwstate_norms = Dict{NTuple{NZ, Int}, Vector{BigInt}}()
    # Number of remaining states in the product representation
    remain_dim = dim1 * dim2

    # Initialize the variables defined above
    for (w, wdim) in remain_inmult
        foundvecs[w] = zeros(CT, wdim, wdim) # Inner multiplicity of product space
    end
    # sort weights from highest to the lowest
    sorted_weights = sort(collect(keys(foundvecs)); lt=rev_less, rev=true)

    # Until the product representation is exhausted
    while !isempty(pout_set)
        # Possible maximal weights
        for w in sorted_weights
            # If there is a remaining state in the product representation
            outer_mult = remain_inmult[w]
            if outer_mult > 0
                # w is a maximal weight of some subrepresentation
                # q-label of the new representation
                qlabel_out = getqlabel(S, w)
                push!(result_qlabels, (qlabel_out, 0))
                @assert qlabel_out ∈ pout_set
                getfull = qlabel_out ∈ qout_set

                output_irep = getNsave_irep(S, RT, qlabel_out; verbose)
                # Already founc vectors in the weight space of 'w'
                foundvecs_w = foundvecs[w]

                # Update innerprod_dict, inv_innerprod_dict
                innerprod_dict[qlabel_out] = output_irep.innerprod
                inv_innerprod_dict[qlabel_out] = output_irep.inv_innerprod

                # Orthogonalize the already-found vectors in the weight w
                # key: qlabels
                # value: (innerprod matrix, its inverse, fac) for weight w
                inprod_info = Dict{NTuple{NZ, Int}, Tuple{Matrix{RT}, Matrix{RT}, Rational{RT}}}()
                for (qlabel, _) in result_qlabels
                    if !haskey(innerprod_dict[qlabel], w) continue end
                    innerprod = innerprod_dict[qlabel][w]
                    inv_innerprod, fac = inv_innerprod_dict[qlabel][w]
                    inprod_info[qlabel] = (innerprod, inv_innerprod, fac)
                end

                output_CGT_block = CG3block{NZ, CT}()
                mwstate_norms[qlabel_out] = Vector{BigInt}()
                nfactors = Vector{Rational{BigInt}}(undef, outer_mult)
                # Outer multiplicity index

                # mwstate is obtained by orthogonalizing the unit vector
                # w.r.t the already-found vectors in the product representation
                # start_index means the nonzero index of starting vector
                start_index = 1

                for oidx=1:outer_mult
                    # Get the maximal weight state and next start_index
                    mwstate::Vector{CT}, start_index, mwnorm = 
                    get_mwstate(irep1, irep2, prod_winfo[w], foundvecs_w, remain_inmult[w], 
                    result_qlabels, inprod_info, mwstate_norms, start_index; verbose)

                    push!(mwstate_norms[qlabel_out], mwnorm)
                    nfactors[oidx] = 1//mwnorm
                    # New vectors: Dict{weight, Matrix{CT}}
                    new_vectors = get_multiplet(S, mwstate, prod_inmult, sl_prod, sz_prod, 
                        output_irep, getfull, wout_set, w; verbose, assertlev)

                    # update the variables defined above
                    # wout: possible weight of the irrep whose maximal weight is 'w'
                    for (wout, (m, n)) in output_irep.Sz
                        #if verbose > 1 println("Processing output weight $(wout)") end
                        inmult = n - m + 1
                        st_idx = prod_inmult[wout] - remain_inmult[wout] + 1
                        if haskey(new_vectors, wout)
                            foundvecs[wout][:, st_idx:st_idx+inmult-1] = new_vectors[wout] 
                        end
                        
                        remain_inmult[wout] -= inmult
                        remain_dim -= inmult

                        # Update the CGT blocks if getfull
                        if getfull
                            winfo = prod_winfo[wout]
                            #if verbose > 1 println("For weight $(wout),") end
                            for ((w1, w2), (ist, ied)) in winfo
                                #if verbose > 1 println("w1 : $(w1), w2 : $(w2)") end
                                block_mat = new_vectors[wout][ist:ied, :]
                                m1, n1 = irep1.Sz[w1]; w1_inmult = n1 - m1 + 1
                                m2, n2 = irep2.Sz[w2]; w2_inmult = n2 - m2 + 1
                                if !haskey(output_CGT_block, (w1, w2, wout))
                                    output_CGT_block[(w1, w2, wout)] = 
                                    zeros(CT, w1_inmult, w2_inmult, inmult, outer_mult)
                                end
                                output_CGT_block[(w1, w2, wout)][:, :, :, oidx] =
                                reshape(block_mat, w1_inmult, w2_inmult, inmult)
                            end
                        end
                    end
                    @assert result_qlabels[end][1] == qlabel_out
                    result_qlabels[end] = (qlabel_out, oidx) 
                end
                if getfull
                    cgt = CGT{S, CT, NZ, 3}((q1, q2, qlabel_out), 
                    output_CGT_block, nfactors, ('+', '+', '-'))
                    result_cg3s[qlabel_out] = cgt
                    if save save_cgt_sqlite(S, cgt) end
                end


                delete!(pout_set, qlabel_out)
                delete!(qout_set, qlabel_out)
                delete!(wout_set, w)
            end
            if isempty(pout_set) break end
        end
    end
    @assert isempty(qout_set); @assert isempty(wout_set)

    return result_cg3s
end


# Overall structure is similar to generate_every_CGT function
function omlookup(::Type{S}, ::Type{RT},
    q1::NTuple{NZ, Int},
    q2::NTuple{NZ, Int};
    verbose=0) where {S<:NonabelianSymm, RT<:Number, NZ}

    @assert nzops(S) == NZ; @assert q1 <= q2
    irep1, irep2 = getNsave_irep(S, RT, q1), getNsave_irep(S, RT, q2)
    dim1, dim2 = dimension(irep1), dimension(irep2)

    om_left = Dict{NTuple{NZ, Int}, Int}()
    for z1 in keys(irep1.Sz)
        for z2 in keys(irep2.Sz)
            out_z = z1 .+ z2
            m1, n1 = irep1.Sz[z1]; inmult1 = n1 - m1 + 1
            m2, n2 = irep2.Sz[z2]; inmult2 = n2 - m2 + 1
            if !haskey(om_left, out_z) om_left[out_z] = 0 end
            om_left[out_z] += inmult1 * inmult2
        end
    end

    # Save the outer multiplicities
    # Key: q-label, Value: outer multiplicity
    om_list = Dict{NTuple{NZ, Int}, Int}()
    remain_dim = dim1 * dim2

    sorted_weights = sort(collect(keys(om_left)); lt=rev_less, rev=true)

    while remain_dim > 0
        for w in sorted_weights
            outer_mult = om_left[w]
            if outer_mult > 0
                qlabel_out = getqlabel(S, w)
                om_list[qlabel_out] = outer_mult
                irep_out = getNsave_irep(S, RT, qlabel_out)

                for (wout, (m, n)) in irep_out.Sz
                    inmult = n - m + 1
                    om_left[wout] -= inmult * outer_mult
                    remain_dim -= inmult * outer_mult
                end
            end
        end
    end
    save_omlist_cg3(S, q1, q2, om_list)
end

to_txt(qlabel::NTuple{NZ, Int}) where NZ = join(qlabel, " ")

qlabels_to_txt(qlabels::NTuple{N, NTuple{NZ, Int}}) where {N, NZ} = 
join([to_txt(q) for q in qlabels], ",")

function save_omlist_cg3(::Type{S},
        q1::NTuple{NZ, Int},
        q2::NTuple{NZ, Int},
        om_list::Dict{NTuple{NZ, Int}, Int}) where {S<:NonabelianSymm, NZ}
    @assert q1 <= q2

    for outspace in keys(om_list)
        omlist_obj = create_omlist_cg3(S, q1, q2, outspace, om_list[outspace])
        save_omlist_sqlite(S, omlist_obj)
    end

    # Create ValidOuts object
    out_spaces = sort(collect(keys(om_list)))
    vo_obj = ValidOuts(S, (q1, q2), out_spaces)
    save_validout_sqlite(S, vo_obj)
end

function create_omlist_cg3(::Type{S},
        q1::NTuple{NZ, Int},
        q2::NTuple{NZ, Int},
        outspace::NTuple{NZ, Int},
        cg3_OM::Int) where {S<:NonabelianSymm, NZ}
    incom_spaces = (q1, q2)
    interm_spaces = Vector{NTuple{0, NTuple{NZ, Int}}}([()])
    cg3_oms = [(cg3_OM,)]
    return OMList(S, incom_spaces, outspace, interm_spaces, cg3_oms)
end



function get_mwstate(irep1::Irep{S, NL, NZ, RT}, 
        irep2::Irep{S, NL, NZ, RT}, 
        weight_info::WeightInfo{NZ},
        foundvecs_w::Matrix{CT},
        remain_inmult::Int,
        result_qlabels::Vector{Tuple{NTuple{NZ, Int}, Int}},
        inprod_info::Dict{NTuple{NZ, Int}, Tuple{Matrix{RT}, Matrix{RT}, Rational{RT}}},
        mwstate_norms::Dict{NTuple{NZ, Int}, Vector{BigInt}},
        start_index::Int; 
        verbose) where {S<:NonabelianSymm, NL, NZ, RT<:Number, CT<:Number}

    if verbose > 2 println(already_found_norms) end
    if verbose > 2 println(start_index) end
    if verbose > 2 println(length(already_found_norms)) end
    # Dimensionality of weight space for product space
    inmult = size(foundvecs_w, 1)
    if verbose > 1 println("\nget_mwstate function\n") end 
    if verbose > 1 println("Total multiplicity : $(inmult), remain_inmult : $(remain_inmult)") end
    if verbose > 1 println("Already found : ") end
    if verbose > 1 display(foundvecs_w) end
    # For each already found vector in the weight space of 'w'
    found = false
    mwstate = zeros(CT, inmult)
    added = zeros(CT, inmult)
    while !found
        innerprod_fac = BigInt(1)
        mwstate = zeros(CT, inmult); mwstate[start_index] = 1
        metric_vector, ist, ied = apply_metric(start_index, irep1, irep2, weight_info; verbose=verbose)
        if verbose > 1 println("Initial mwstate : ", mwstate) end

        jst = 1
        gotzero = false
        for (qlabel, ow) in result_qlabels
            @assert length(mwstate_norms[qlabel]) == ow
            if !haskey(inprod_info, qlabel) continue end
            innerprod, inv_innerprod, fac_inprod = inprod_info[qlabel]
            iw = size(innerprod, 1)

            # TODO: document this part. To reduce the intermediate integer factors,
            # I reordered the calculation so it would be hard to understand at once.
            for i=1:ow
                #println("qlabel, i:", qlabel, i)
                mwnorm = mwstate_norms[qlabel][i]
                jed = jst + iw - 1
                overlaps = (metric_vector' * foundvecs_w[ist:ied, jst:jed])' 
                #println(fac_inprod)
                #println(mwnorm)
                coeff, overlap_fac = prodNcfac(inv_innerprod, overlaps, fac_inprod, 1//mwnorm)
                
                fac_ratio = innerprod_fac * overlap_fac

                for j=jst:jed added -= coeff[j - jst + 1] * foundvecs_w[:, j] end

                mwstate *= fac_ratio.den; innerprod_fac *= fac_ratio.den
                mwstate += added * fac_ratio.num
                # Reset 'added' vector
                added .= 0
                if iszero(mwstate) gotzero=true; break end
                divfac = gcd(gcd(mwstate), innerprod_fac)
                mwstate = div.(mwstate, divfac)
                innerprod_fac = div(innerprod_fac, divfac)
                jst += iw
                #println(mwstate)
            end
            if gotzero break end

        end
        if gotzero start_index += 1; continue
        else found = true end
    end

    #println("Final mwstate : ", mwstate)
    _, mwstate = divcfac(mwstate)
    mwstate_norm::BigInt = inprod(BigInt.(mwstate), BigInt.(mwstate), irep1.innerprod, 
        irep2.innerprod, weight_info, nothing; verbose=verbose)
    first_nzidx = findfirst(!iszero, mwstate)
    # Make the first nonzero element positive
    if mwstate[first_nzidx] < 0 mwstate = -mwstate end
    return mwstate, start_index + 1, mwstate_norm
end

function prodNcfac(A::AbstractArray,
    B::AbstractArray,
    afac::Rational,
    bfac::Rational)

    @assert afac.num == 1 && bfac.num == 1
    arr, fac = A * B, afac * bfac
    divfac = gcd(gcd(arr), fac.den)
    nfac = fac * divfac; @assert nfac.num == 1
    return div.(arr, divfac), nfac
end

function apply_metric(start_index::Int, 
        irep1::Irep{S, NL, NZ, RT}, 
        irep2::Irep{S, NL, NZ, RT}, 
        weight_info::WeightInfo{NZ};
        verbose) where {S<:NonabelianSymm, NL, NZ, RT<:Number}
    z1, z2, ist, ied, nzi = get_zinfo(start_index, weight_info)
    inner_prod1 = irep1.innerprod[z1]; l1 = size(inner_prod1)[1]
    inner_prod2 = irep2.innerprod[z2]; l2 = size(inner_prod2)[1]
    cidx1, cidx2 = (nzi - 1) % l1 + 1, div(nzi - 1, l1) + 1
    return view(inner_prod1, :, cidx1) ⊗ view(inner_prod2, :, cidx2), ist, ied
end

function get_zinfo(start_index::Int, 
        weight_info::WeightInfo{NZ}) where {NZ}
    # Get the z-values and length of the vector
    # which is used to apply the metric
    for ((z1, z2), (ist, ied)) in weight_info
        if start_index >= ist && start_index <= ied
            len = ied - ist + 1
            return z1, z2, ist, ied, start_index - ist + 1
        end
    end
    error("Invalid start index $(start_index) for weight info $(weight_info)")
end


function orthogonalize(v::Vector, 
        exv::Vector, 
        overlap::Integer, 
        exvnorm::Integer; 
        verbose) 
    factor = overlap // exvnorm
    v = factor.den * v - factor.num * exv
    cfac, v = divcfac(v)
    return factor.den, cfac, v
end

# If irep corresponding to new multiplet in the product representation already exists,
# use the lowering operator of existing irep to get multiplet vectors.
# This has two advantages
# 1. Faster calculation is possible since we know how to get another vectors
# from lowering operators of existing irep.
# 2. We can prevent the overflow of integer values.

function get_multiplet(::Type{S}, 
        mwstate::Vector{CT},
        prod_inmult::Dict{NTuple{NZ, Int}, Int},
        sl_prod::NTuple{NL, Dict{NTuple{NZ, Int}, SparseMatrixCSC{RT}}},
        sz_prod::Dict{NTuple{NZ, Int}, Tuple{Int, Int}},
        output_irep::Irep{S, NL, NZ, RT},
        getfull::Bool,
        wout_set::Set{NTuple{NZ, Int}},
        mw::NTuple{NZ, Int};
        verbose=0, 
        assertlev=0) where {S<:NonabelianSymm, NL, NZ, RT<:Number, CT<:Number}

    output_vectors = Dict{NTuple{NZ, Int}, Matrix{CT}}()
    output_vectors[mw] = reshape(mwstate, length(mwstate), 1)

    wout_set = getfull ? Set(collect(keys(output_irep.Sz))) : wout_set
    for w in wout_set
        # Get states in the weight space 'w' recursively
        if haskey(output_irep.Sz, w)
            get_wspace!(S, output_vectors, w, prod_inmult, sl_prod, output_irep; verbose)
        end
    end

    if getfull && assertlev > 0
        check_vanishing(output_vectors, sl_prod, output_irep; verbose)
        println("multiplet $(getqlabel(S, mw)) passed")
    end
    return output_vectors
end

function get_wspace!(::Type{S}, 
    output_vectors::Dict{NTuple{NZ, Int}, Matrix{CT}},
    targetw::NTuple{NZ, Int},
    prod_inmult::Dict{NTuple{NZ, Int}, Int},
    sl_prod::NTuple{NL, Dict{NTuple{NZ, Int}, SparseMatrixCSC{RT}}},
    output_irep::Irep{S, NL, NZ, RT};
    verbose) where {S<:NonabelianSymm, NL, NZ, RT<:Number, CT<:Number}

    # If the target weight space is already filled, return
    if haskey(output_vectors, targetw) return end

    # Else (weight space not found yet)
    dzs = getdzs(S)
    size1 = prod_inmult[targetw]
    m, n = output_irep.Sz[targetw]; size2 = n - m + 1
    output_vectors[targetw] = zeros(CT, size1, size2)
    remaining_set = Set(1:size2) # Set of indices of vectors to be filled
    #if verbose > 1 println(remaining_set) end

    # operator index
    for opidx in 1:NL
        init_weight = Tuple(targetw .+ dzs[opidx])
        if haskey(output_irep.Sl[opidx], init_weight)
            get_wspace!(S, output_vectors, init_weight, prod_inmult, sl_prod, output_irep; verbose)
            #if verbose > 1 println("Applying operator $(opidx)") end
            lop_output_irep = output_irep.Sl[opidx][init_weight]
            for j in axes(lop_output_irep, 2)
                jth_colm = lop_output_irep[:, j]
                if nnz(jth_colm) == 1
                    nzind, nzval = jth_colm.nzind[1], jth_colm.nzval[1]
                    if nzind in remaining_set
                        init_vec = output_vectors[init_weight][:, j]
                        output_vectors[targetw][:, nzind] = 
                        div.(sl_prod[opidx][init_weight] * init_vec, nzval)
                        delete!(remaining_set, nzind)
                    end
                end
            end
        end
    end
    @assert isempty(remaining_set) 
end

function check_vanishing(output_vectors::Dict{NTuple{NZ, Int}, Matrix{CT}},
        sl_prod::NTuple{NL, Dict{NTuple{NZ, Int}, SparseMatrixCSC{RT}}},
        output_irep::Irep{S, NL, NZ, RT};
        verbose) where {S<:NonabelianSymm, NL, NZ, RT<:Number, CT<:Number}
    for opi in 1:NZ
        if verbose > 1 println("Testing lowering operator $(opi)") end
        prod_lop = sl_prod[opi]
        irep_lop = output_irep.Sl[opi]
        
        # For each weight in the output representation
        for w in keys(output_irep.Sz)
            if !haskey(irep_lop, w)
                states = output_vectors[w]
                if haskey(prod_lop, w)
                    if verbose > 1
                        println(w)
                        display(states)
                        display(prod_lop[w])
                        println("\n========================\n")
                    end
                    state_after_lop = prod_lop[w] * states
                    @assert iszero(state_after_lop)
                end
            end
        end
    end
    println("Passed the vanishing test for all lowering operators")
end


function ⊗(A::SparseArray{T, N}, B::SparseArray{S, N}) where {T, S, N}
    sizeA = size(A)
    sizeB = size(B)
    
    # Define the size of the resulting sparse array C
    C_size = sizeA .* sizeB
    C = SparseArray{promote_type(T, S), N}(undef, C_size)
    
    for (IA, vA) in nonzero_pairs(A)
        for (IB, vB) in nonzero_pairs(B)
            # Get new indices and new value
            tuple_IC = (Tuple(IB) .- 1) .* sizeA .+ Tuple(IA)
            IC = CartesianIndex(tuple_IC)
            
            C[IC] = vA * vB
        end
    end
    return C
end
⊗(a::AbstractMatrix, b::AbstractMatrix) = kron(b, a)
⊗(a::AbstractVector, b::AbstractVector) = kron(b, a)

# divide vector by common factor (gcd) of nonzero elements
# so that nonzero elements of input vector become relative prime each other
function divcfac(v::SparseVector{<:Integer}) 
	if iszero(v) return 1, v end
	cfac = gcd(v.nzval)
	cfac, div.(v, cfac) 
end

function divcfac(v::Vector{<:Integer})
    if iszero(v) return 1, v end
    cfac = gcd(v)
    return cfac, div.(v, cfac)
end

function is_fundamental(::Type{S}, 
    qlabel::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ}

    if sum(qlabel) <= 1 return true end
    if isSON(S) && sum(qlabel) == 2
        N = defirepdim(S)
        if N % 2 == 1 && qlabel[end] == 2 return true end
        if N % 2 == 0 && (qlabel[end-1] == 2 || qlabel[end] == 2) return true end
    end
    return false
end

function split_qlabel(::Type{S}, 
    qlabel::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ}
    @assert NZ == nzops(S)
    # Split the q-label into two smaller q-labels
    # For example, (a, b, 0, 0) -> (a, b-1, 0, 0) and (0, 1, 0, 0)
    lastidx = findlast(!iszero, qlabel)
    smallq = [0 for _=1:NZ]
    if isSON(S)
        N = defirepdim(S)
        if N % 2 == 1 && lastidx == NZ
            @assert qlabel[NZ] % 2 == 0; smallq[NZ] = 2
        elseif N % 2 == 0 && lastidx >= NZ-1
            @assert (qlabel[NZ-1] + qlabel[NZ]) % 2 == 0
            if lastidx == NZ && qlabel[NZ] % 2 == 1
                smallq[NZ-1] = 1; smallq[NZ] = 1
            else smallq[lastidx] = 2 end
        else smallq[lastidx] = 1 end
    else smallq[lastidx] = 1 end
    smallq = Tuple(smallq)
    bigq = qlabel .- smallq
    return bigq, smallq
end

function getNsave_irep(::Type{S}, 
    ::Type{RT}, 
    qlabel::NTuple{NZ, Int}; 
    verbose=0) where {S<:AbelianSymm, RT<:Number, NZ} 

    @assert nzops(S) == NZ
    Sl = (); Sz = Dict{NTuple{NZ, Int}, Tuple{Int, Int}}()
    Sz[qlabel] = (1, 1)

    innerprod = Dict{NTuple{NZ, Int}, Matrix{RT}}()
    innerprod[qlabel] = ones(RT, 1, 1)

    inv_innerprod = Dict{NTuple{NZ, Int}, Tuple{Matrix{RT}, Rational{RT}}}()
    inv_innerprod[qlabel] = (ones(RT, 1, 1), one(RT))
    return Irep(S, RT, Sl, Sz, innerprod, inv_innerprod, 1)
end

# Recursively get and save ireps. 
function getNsave_irep(::Type{S}, 
    ::Type{RT}, 
    qlabel::NTuple{NZ, Int}; 
    verbose=0) where {S<:NonabelianSymm, RT<:Number, NZ} 

    # If it exists, load the representation and return
    irep = load_irep_sqlite(S, RT, qlabel; verbose) 
    if !isnothing(irep) @assert isa(irep, Irep); return irep end
    if is_fundamental(S, qlabel)
        get_fundamental_ireps(S, RT; verbose=verbose)
        return load_irep_sqlite(S, RT, qlabel; verbose)
    end

    # If not, build the representation from smaller ones
    bigq, smallq = split_qlabel(S, qlabel)
    irep1 = getNsave_irep(S, RT, bigq; verbose) 
    irep2 = getNsave_irep(S, RT, smallq; verbose)
    # Get largest irep from the product rep of irep1 and irep2
    # and return the largest irep in the product rep
    # sl_prod, sz_prod : lowering operators and z-operators of the product representation
    # prod_weight_info: information about weights of product representation
    sl_prod, sz_prod, _, prod_weight_info = tensor_prod_reps(irep1, irep2; verbose)
    irep_new = get_irep_frommw(S, sl_prod, sz_prod, 
        irep1.innerprod, irep2.innerprod, nothing, prod_weight_info; verbose)
    save_irep_sqlite(irep_new; verbose)
    return irep_new
end

# Get the largest irreducible representation from 
# This function is called when 1) getting the alternating product representation,
# or 2) getting the largest irreducible representation from the product representation
# In former case, inner_prod1 is just a dictionary of matrices,
# and inner_prod2, weight_info are 'nothing'.
# In latter case, inner_prod1 and inner_prod2 are inner product matrices of 
# the original irreps, and weight_info is a dictionary of information
function get_irep_frommw(::Type{S}, 
        Sl::NTuple{NL, Dict{NTuple{NZ, Int}, SparseMatrixCSC{RT}}},
        Sz::Dict{NTuple{NZ, Int}, Tuple{Int, Int}},
        inner_prod1,
        inner_prod2,
        mw::Union{NTuple{NZ, Int}, Nothing},
        weight_info::Union{Nothing, Dict{NTuple{NZ, Int}, WeightInfo{NZ}}};
        verbose) where {S<:NonabelianSymm, RT<:Number, NL, NZ}

    # ordered (from highest to lowest) z-values of the original representation
    zvals_ordered = sorted_zvals(Sz)
    mw = isnothing(mw) ? zvals_ordered[1] : mw # Maximal weight of the original representation
    qlabel = getqlabel(S, mw) # q-label of the original representation
    ninnerprod = Dict{NTuple{NZ, Int}, Matrix{RT}}()    

    # Step 1. Get the basis of each weight space of the multiplet
    # Store the basis of the multiplet in this dictionary
    # Keys are just the weights 
    weight_basis_dict = Dict{NTuple{NZ, Int}, Matrix{RT}}()
    # Put maximal weight state in the dictionary
    weight_basis_dict[Tuple(mw)] = reshape(RT[1], (1, 1))
    # The dimension of weight spaces
    weight_cnt = Dict{NTuple{NZ, Int}, Int}()
    tabsNvecs = Dict{Tableau{S}, Vector{RT}}()

    # Assume that the maximal weight has inner multiplicity 1
    mwstate = RT[1]
    mw_tableau = mw_tableau_qlabel(S, qlabel)
    @assert get_weight(mw_tableau) == mw
    tabsNvecs[mw_tableau] = mwstate

    weight_cnt[mw] = 1
    cops = get_crystal_ops(S)
    tabs_queue = [mw_tableau]

    while !isempty(tabs_queue)
        new_queue = Tableau{S}[]
        for tab in tabs_queue
            for i in 1:NZ
                if det_getvector(tab, i, cops, tabsNvecs)
                    new_tab, new_vec = tab, tabsNvecs[tab]
                    while true
                        ow = get_weight(new_tab)
                        new_tab = apply_lowering(new_tab, cops, i)
                        if isnothing(new_tab) break end

                        new_vec = Sl[i][ow] * new_vec
                        _, new_vec = divcfac(new_vec)
                        if haskey(tabsNvecs, new_tab) continue end

                        nw = get_weight(new_tab)
                        if !haskey(weight_cnt, nw) weight_cnt[nw] = 0 end
                        weight_cnt[nw] += 1

                        push!(new_queue, new_tab)
                        tabsNvecs[new_tab] = new_vec
                    end
                end
            end
        end
        tabs_queue = new_queue
    end

    nsaved_vecs = Dict{NTuple{NZ, Int}, Int}()
    for (tab, vec) in tabsNvecs
        w = get_weight(tab); len = length(vec)
        if !haskey(weight_basis_dict, w)
            weight_basis_dict[w] = zeros(RT, len, weight_cnt[w])
        end
        if !haskey(nsaved_vecs, w) nsaved_vecs[w] = 1 end
        weight_basis_dict[w][:, nsaved_vecs[w]] = vec
        nsaved_vecs[w] += 1
    end

    # Step 2. Sort the weights and get the matrix elements of the lowering operators
    # New list of weights
    nSz = Dict{NTuple{NZ, Int}, Tuple{Int, Int}}()
    weights = collect(keys(weight_basis_dict))
    sort!(weights; lt=(a, b) -> less_weight(S, a, b), rev=true)

    # Fill out the weight_inds dictionary
    # nd : the dimension of the new representation
    nd = 0
    for w in weights
        weight_basis = weight_basis_dict[w]
        wdim = size(weight_basis)[2]
        ninnerprod[w] = get_innerprodmat(S, RT, w, weight_basis, wdim, 
            inner_prod1, inner_prod2, weight_info; verbose)
        nSz[w] = (nd+1, nd+wdim)
        nd += wdim
    end
    if verbose > 1 println("List of weights of new irep: $(nSz)") end

    # Then, get the matrix elements of the lowering operators
    # nSl : new Sl operators
    nSl = Tuple(Dict{NTuple{NZ, Int}, SparseMatrixCSC{RT}}() for _=1:NL)

    # For each weight in the new irep
    for w in weights
        # weight_basis.n is the dimension of the weight space
        weight_basis = weight_basis_dict[w]
        inmult_old = size(weight_basis)[2]
        for vind in 1:inmult_old
            # Get the vector of the weight space
            v = weight_basis[:, vind]
            if verbose > 1 println("Weight $(w), vector $(vind) : ", v) end

            for opi = 1:NL
                # If there is no lowering operator for this weight, skip
                if !haskey(Sl[opi], w) continue end
                # Apply the lowering operator to the vector
                v_new = Sl[opi][w] * v

                if verbose > 1 println("Applying Sl[$(opi)] to $(v)") end
                if verbose > 1 println("Result : ", v_new) end

                if iszero(v_new) continue end

                # Express the new vector as the linear combination
                # of the basis vectors of the new weight space
                lowered_w = w .- Tuple(getdzs(S)[opi])
                nweight_basis = weight_basis_dict[lowered_w]
                coeffs = get_linear_coeff(v_new, nweight_basis; verbose=verbose)

                # According to the obtained coefficients
                # fill the matrix element of the lowering operator
                inmult_new = size(nweight_basis)[2]
                wold = Tuple(w)
                if !haskey(nSl[opi], wold)
                    # If the weight is not in the dictionary, create a new entry
                    nSl[opi][wold] = spzeros(RT, inmult_new, inmult_old)
                end
                nSl[opi][wold][:, vind] = coeffs
            end
        end
    end

    inv_innerprod = get_inv_innerprod(S, ninnerprod)
    
    # Return the new irep, which is the largest irreducible representation
    # constructed from the product representation
    # of the original representation and the trivial representation
    dim = dim_from_innerprod(ninnerprod)
    return Irep(S, RT, nSl, nSz, ninnerprod, inv_innerprod, dim) 
end

function get_inv_innerprod(::Type{S},
    ninnerprod::Dict{NTuple{NZ, Int}, Matrix{RT}}; 
    verbose=0) where {S<:NonabelianSymm, RT<:Number, NZ}

    inv_innerprod = Dict{NTuple{NZ, Int}, Tuple{Matrix{RT}, Rational{RT}}}()

    for (zval, inprod_matrix) in ninnerprod
        if all(x->x>=0, getqlabel(S, zval))
            sz = size(inprod_matrix)[1]
            M = matrix(Nemo.QQ, sz, sz, inprod_matrix)
            Minv = Nemo.inv(M)

            fac = lcm(Matrix{BigInt}((denominator.(Minv))))
            B = Matrix{RT}(map_entries(x->Nemo.ZZ(x*fac), Minv))
            inv_innerprod[zval] = (B, 1 // fac)
        end
    end

    return inv_innerprod
end


function get_TNmat(inprod_matrix::Matrix{RT}; 
        verbose=0) where {RT<:Number}
    Nmat = inprod_matrix
    if verbose > 2 println(Nmat) end
    sz = size(inprod_matrix)[1]
    Tmat = Matrix{RT}(I, sz, sz)
    # Orthogonalize jth vector with respect to ith vector
    for j=1:sz
        for i=1:j-1
            transf_mat = sparse(RT, I, sz, sz)
            factor = Nmat[i, j] // Nmat[i, i]
            transf_mat[i, j] = -factor.num
            transf_mat[j, j] = factor.den
            if verbose > 2 display(transf_mat) end
            Tmat = Tmat * transf_mat
            Nmat = transpose(transf_mat) * Nmat * transf_mat
            if verbose > 2 println("Nmat:"); display(Nmat) end
            if verbose > 2 println("Tmat:"); display(Tmat) end
        end
        cfac = gcd(Tmat[:, j])
        if verbose > 2 println("cfac:", cfac) end
        if cfac > 1
            Tmat[:, j] = div.(Tmat[:, j], cfac)
            
            @assert iszero(Nmat[:, j] .% cfac)
            Nmat[:, j] = div.(Nmat[:, j], cfac)
            @assert iszero(Nmat[j, :] .% cfac)
            Nmat[j, :] = div.(Nmat[j, :], cfac)
        end
    end
    @assert isdiag(Nmat)
    return Tmat, diag(Nmat)
end

# Currently this function is not used.
# It was used to calculate TNT and its inverse.
inv_mat(mat::Matrix{RT}; verbose=0) where RT<:Number = 
(s = size(mat, 1); solve_gaussian(s, copy(mat), Matrix{RT}(I, s, s); verbose=verbose))

# Currently this function is not used.
# It was used to calculate TNT and its inverse.
function get_Ninv(Nmatdiag::Vector{RT}) where {RT<:Number}
    comm_mult = lcm(Nmatdiag)
    return [div(comm_mult, x) for x in Nmatdiag], comm_mult
end

# Get 1j-symbol for the given q-label
# nfac in this function is integer in this case
# It is inversed in the very last step
# TODO: If spin representation is implemented, generalize it
function getNsave_1jsym(::Type{S},
        ::Type{RT},
        ::Type{CT},
        q::NTuple{NZ, Int};
        verbose=0) where {S<:NonabelianSymm, RT<:Number, CT<:Number, NZ}
    nfac = 1
    dualq = get_dualq(S, q)
    permute = q > dualq
    qlabels = (q, dualq)

    # If it exists, load and return
    if permute
        new_blocks = Dict{NTuple{2, NTuple{NZ, Int}}, Array{CT}}()
        cgt = getNsave_1jsym(S, RT, CT, dualq; verbose)
        @assert cgt.qlabels == (dualq, q)
        for k in keys(cgt.blocks)
            new_key = (k[2], k[1])
            new_blocks[new_key] = permutedims(cgt.blocks[k], (2, 1, 3))
        end
        return CGT{S, CT, NZ, 2}(qlabels, new_blocks, cgt.nfactor, ('+', '+'), cgt.size_byte) 
    end

    @assert q <= dualq
    cgt = load_cgt_sqlite(S, CT, qlabels)

    # If not, generate it
    irep = getNsave_irep(S, RT, q; verbose=verbose)
    dual_irep = getNsave_irep(S, RT, dualq; verbose=verbose)

    dzs = getdzs(S)
    output_blocks = Dict{NTuple{2, NTuple{NZ, Int}}, Array{CT, 3}}()
    # sort weights from highest to lowest
    # when comparing two weights, last elements are compared
    # if they are the same, compare next-last elements
    sorted_weights = sort(collect(keys(irep.Sz)); lt=rev_less, rev=true)
    mw = sorted_weights[1]
    for (i, w) in enumerate(sorted_weights)
        if verbose > 1 println("\n\nGetting block correspond to ", w) end
        minw = Tuple(-x for x in w) # Minus w
        # Fill first block to [1]
        if i == 1 
            output_blocks[(w, minw)] = reshape(CT[1], 1, 1, 1)
            continue 
        end
        m, n = irep.Sz[w]; mult = n - m + 1

        # Let M be the block we want to fill.
        # It satisfies matrix equation MA = B
        A = zeros(RT, mult, 0)
        B = zeros(RT, mult, 0)
        for opi in 1:NZ # (lowering) operator index
            initw = Tuple(w .+ dzs[opi])
            mininitw = Tuple(-x for x in initw) # Minus initw
            if haskey(irep.Sl[opi], initw)
                pA = transpose(dual_irep.Sl[opi][minw])
                #println(irep.Sl[opi][initw])
                #println(output_blocks[(initw, mininitw, zw)][:, :, 1])
                pB = -irep.Sl[opi][initw] * output_blocks[(initw, mininitw)][:, :, 1]
                A = hcat(A, pA); B = hcat(B, pB)
            end
            @assert size(A) == size(B)
        end
        # M : matrix, n : denominator introduced in the Gaussian elimination
        M, n = solve_gaussian(mult, A, B)
        if verbose > 1
            println("M : ")
            display(M)
            println("n : ", n) 
        end
        for key in keys(output_blocks)
            output_blocks[key] *= n
        end
        nfac *= n^2
        Mvec = reshape(M, mult*mult)
        nfac += partial_inprod(Mvec, Mvec, 
            irep.innerprod[w], dual_irep.innerprod[minw]; verbose) 
        # Fill in the new block
        output_blocks[(w, minw)] = reshape(M, mult, mult, 1) 
    end
    nfac = 1//nfac; den = nfac.den
    cgt = CGT{S, CT, NZ, 2}(qlabels, output_blocks, [nfac], ('+', '+'))
    # If q == dualq, make a R-symbol 
    if q == dualq
        mwblk = output_blocks[(mw, .-mw)]
        @assert size(mwblk) == (1, 1, 1)
        opposite_blk = output_blocks[(.-mw, mw)]
        @assert size(opposite_blk) == (1, 1, 1)
        @assert abs(opposite_blk[1]) == mwblk[1]

        zero_qlab = Tuple(0 for _ in 1:NZ)
        rsym_mat = reshape(CT[sign(opposite_blk[1])*den], 1, 1)
        rsym_struct = Rsymbol{S, CT, NZ}(q, zero_qlab, rsym_mat, [nfac])
        save_Rsymbol_sqlite(rsym_struct)
    end
    zq = Tuple(0 for _ in 1:NZ) # zero qlabel
    save_cgt_sqlite(S, cgt)
    return cgt
end

function load_1jblk(::Type{S},
    ::Type{RT},
    ::Type{CT},
    q::NTuple{NZ, Int}) where {S<:NonabelianSymm, RT<:Number, CT<:Number, NZ}

    sym = getNsave_1jsym(S, RT, CT, q)
    return sym.blocks, sym.nfactor
end

# Solve Ax = b. To use solve_gaussian function, transpose it
function solve_linsys(A::AbstractMatrix{RT},
    b::AbstractVector{RT};
    verbose=0) where {RT<:Number}

    _, n = size(A)
    A_t = transpose(A)
    b_t = transpose(b)
    
    x, xfac = solve_gaussian(n, A_t, b_t; verbose=verbose)
    return reshape(x, n), xfac
end


# Solve the matrix equation MA = B
# M is the matrix we want to fill
# n is the denominator introduced in the Gaussian elimination
function solve_gaussian(mult::Int, 
    A::AbstractMatrix{RT}, 
    B::AbstractMatrix{RT};
    verbose=0) where {RT<:Number}

    if verbose > 1
        println("A : ")
        display(A)
        println("B : ")
        display(B)
    end
    @assert size(A)[1] == mult 
    ncol = size(A)[2]
    n = 1
    for i=1:mult
        if verbose > 1 println("Processing row $(i)") end
        # Find the pivot element
        for j=i:ncol
            if A[i, j] != 0
                # Swap the columns
                if j != i
                    A[:, [i, j]] = A[:, [j, i]]
                    B[:, [i, j]] = B[:, [j, i]]
                end
                break
            end
        end

        pivot = A[i, i]
        for j=1:ncol
            if j == i continue end
            if A[i, j] == 0 continue end
            if verbose > 1 println("Processing column $(j) with pivot $(pivot)") end
            frac = pivot // A[i, j]
            A[:, j] = frac.num * A[:, j] - frac.den * A[:, i]
            B[:, j] = frac.num * B[:, j] - frac.den * B[:, i]
        end
        if verbose > 1
            println("After processing row $(i):")
            println("A : ")
            display(A)
            display(B)
        end
    end
    @assert isdiag(A[:, 1:mult])
    Adiag = diag(A[:, 1:mult])
    M = B[:, 1:mult]
    if verbose > 1 println("Adiag : ", Adiag) end
    if verbose > 1 println("M : ", M) end
    n = lcm(Adiag)
    for i=1:mult M[:, i] *= div(n, Adiag[i]) end
    Mgcd = gcd(M)
    return div.(M, gcd(Mgcd, n)), div(n, gcd(Mgcd, n))
end



function get_innerprodmat(::Type{S}, ::Type{RT}, 
        zval::NTuple{NZ, Int},
        weight_basis::Matrix{RT}, 
        wdim::Int, 
        inner_prod1, 
        inner_prod2, 
        weight_info; verbose) where {S<:NonabelianSymm, RT<:Number, NZ}
    inprod_matrix = zeros(RT, wdim, wdim)
    winfo = isnothing(weight_info) ? nothing : weight_info[zval]
    for i=1:wdim
        for j=i:wdim
            inprod_matrix[i, j] = inprod(weight_basis[:, i], weight_basis[:, j], 
                inner_prod1, inner_prod2, winfo, zval; verbose=verbose)
            if i != j inprod_matrix[j, i] = inprod_matrix[i, j] end
        end
    end
    return inprod_matrix
end

inprod(v1::Vector{RT}, 
    v2::Vector{RT}, 
    basis_innerprod::Dict{NTuple{NZ, Int}, Matrix{RT}},
    ::Nothing,
    ::Nothing, 
    zval::NTuple{NZ, Int};
    verbose=0) where {RT<:Number, NZ} =
dot(v1, basis_innerprod[zval] * v2)

function inprod(V1::Vector, 
    V2::Vector, 
    inner_prod1::Dict{NTuple{NZ, Int}, Matrix{RT}},
    inner_prod2::Dict{NTuple{NZ, Int}, Matrix{RT}},
    winfo::WeightInfo{NZ},
    ::Any;
    verbose=0) where {RT<:Number, NZ} 
    if verbose > 1 println("Calculating inner product of $(V1) and $(V2)") end
    if verbose > 1 println("inner product of 1st space: $(inner_prod1)") end
    if verbose > 1 println("inner product of 2nd space: $(inner_prod2)") end
    if verbose > 1 println("winfo : $(winfo)") end

    result = 0
    for ((z1, z2), (ist, ied)) in winfo
        result += partial_inprod(V1[ist:ied], V2[ist:ied], 
            inner_prod1[z1], inner_prod2[z2]; verbose=verbose)
    end
    return result
end

# Get V1^T * (inner_prod1 ⊗ inner_prod2) * V2
# TODO: Reshape V1 and V2 to matrices and do matrix multiplication
# TODO: Instead of general matrix type, use symmetric matrix type
function partial_inprod(V1::Vector{VT}, 
    V2::Vector{VT}, 
    inner_prod1::Matrix{MT}, 
    inner_prod2::Matrix{MT}; 
    verbose=0) where {VT<:Integer, MT<:Integer}
    s1 = size(inner_prod1)[1]
    s2 = size(inner_prod2)[1]
    @assert length(V1) == length(V2)
    @assert size(inner_prod1)[1] * size(inner_prod2)[1] == length(V1)
    V1c = copy(V1)
    # Apply inner_prod1⊗I to V1c
    for i=1:s2
        range = (i-1)*s1+1:i*s1
        V1c[range] = inner_prod1 * V1c[range]
    end
    # Apply I⊗inner_prod2 to V1c
    for i=1:s1
        range = i:s1:s1*s2
        V1c[range] = inner_prod2 * V1c[range]
    end
    # Then, V1c = (inner_prod1 ⊗ inner_prod2) * V1
    return dot(V2, V1c)
end

# Get the coefficients when the vector v is expressed as a 
# linear combination of the basis vectors of the weight space
function get_linear_coeff(v::Vector{RT}, 
        basis::Matrix{RT}; verbose=0) where RT<:Number
    #println("Getting linear coefficients of the vector ")
    #println(v)
    #println(basis)
    inmult = size(basis)[2]
    if verbose > 1 
        println("\n\nGetting linear coefficients of $(v) in the basis ") 
        for i=1:inmult println(basis[:, i]) end
    end
    
    # Set of indices where the coefficients are not determined yet
    not_determined = Set(1:inmult)
    # nonzero indices of the basis vectors
    nzind_dict = Dict{Int, Set{Int}}()
    for i in not_determined
        nzind_dict[i] = Set(findall(!iszero, basis[:, i]))
    end

    # Coefficients vector. Will be filled with the coefficients
    coeffs = zeros(RT, inmult)

    # Iterate until all coefficients are determined
    while !iszero(v)
        toremove = Set{Int}()
        # Find the basis vector which has unique nonzero index
        for i in not_determined
            uniqind = find_unique_nzind(nzind_dict, i, not_determined)
            if uniqind != 0
                if verbose > 1 println("Found unique nonzero index $(uniqind) in basis vector $(i)") end
                # Get the coefficient of the basis vector
                # TODO: handle the case if the coefficient is not integer
                @assert v[uniqind] % basis[uniqind, i] == 0 
                coeffs[i] = div(v[uniqind], basis[uniqind, i])
                
                v -= coeffs[i] * basis[:, i]
                push!(toremove, i)
            end
        end
        if isempty(toremove) break end
        not_determined = setdiff(not_determined, toremove)
    end
    if !iszero(v)
        nd_vec = sort(collect(not_determined))
        # Solve the remaining linear equation Ax = v
        x, xfac = solve_linsys(basis[:, nd_vec], v; verbose=verbose)
        @assert xfac == 1
        coeffs[nd_vec] = x
    end
    #println(coeffs)
    return coeffs
end

function find_unique_nzind(nzind_dict::Dict{Int, Set{Int}}, i::Int, 
        not_determined::Set{Int}) 
    # If the basis vector has unique nonzero index, return true
    other_inds = Set{Int}()
    for j in not_determined
        if j == i continue end
        union!(other_inds, nzind_dict[j])
    end
    uniqueinds = setdiff(nzind_dict[i], other_inds) 
    return isempty(uniqueinds) ? 0 : first(uniqueinds)
end

# Input : two representation. They should be Irep type
# with same type parameters (S, NL, NZ, RT)
# Return the new representation which is the tensor product of two representations.
function tensor_prod_reps(r1::Irep{S, NL, NZ, RT},
	r2::Irep{S, NL, NZ, RT};
	verbose=0) where {S, NL, NZ, RT<:Number}

    if verbose > 1 println("tensor_prod_reps from $(totxt(r1)), $(totxt(r2))") end
    r1_zvals = sorted_zvals(r1.Sz)
    r2_zvals = sorted_zvals(r2.Sz)
    prod_inmult = Dict{NTuple{NZ, Int}, Int}()
    nSz = Dict{NTuple{NZ, Int}, Tuple{Int, Int}}()
    prod_weight_info = Dict{NTuple{NZ, Int}, WeightInfo{NZ}}()
    # Get the weight space of the product representation
    for z1 in r1_zvals
        for z2 in r2_zvals
            # Get the weight info of the product representation
            # m1, m2: inner multiplicity of weights z1, z2
            a1, b1 = r1.Sz[z1]; m1 = b1 - a1 + 1
            a2, b2 = r2.Sz[z2]; m2 = b2 - a2 + 1
            zsum = z1 .+ z2
            if !haskey(prod_inmult, zsum) prod_inmult[zsum] = 0 end
            if !haskey(prod_weight_info, zsum)
                prod_weight_info[zsum] = WeightInfo{NZ}()
                prod_weight_info[zsum][(z1, z2)] = (1, m1*m2)
            else
                dw = prod_inmult[zsum]
                prod_weight_info[zsum][(z1, z2)] = (dw+1, dw+m1*m2)
            end
            prod_inmult[zsum] += m1 * m2
        end
    end

    # dimension of the product representation
    dp = 0
    for zval in sorted_zvals(prod_weight_info)
        inmult = prod_inmult[zval]
        nSz[zval] = (dp+1, dp+inmult)
        dp += inmult
    end
    if verbose > 1 println("Inner multiplicities are $(prod_inmult)") end 
    if verbose > 1 println("Weight of product rep is $(nSz)") end 
    
    # Get the lowering operators of the product representation
    nSl = Tuple(Dict{NTuple{NZ, Int}, SparseMatrixCSC{RT}}() for _=1:NL)
    dzs = getdzs(S)
    # For each lowering operator
    for opi in 1:NL
        if verbose > 1 println("constructing $(opi)th lowering operator") end
        # For each initial weight in the product representation
        for (zval_init, (a, b)) in nSz
            mult_dom = b - a + 1
            # final weight after applying the lowering operator
            zval_final = Tuple(zval_init .- dzs[opi])
            if !haskey(nSz, zval_final) continue end
            if verbose > 1 println("There is matrix element from $(zval_init) to $(zval_final)") end
            c, d = nSz[zval_final]; mult_final = d - c + 1
            mat = spzeros(RT, mult_final, mult_dom)
            
            # Fill in the matrix element of the lowering operator
            # For each weight info of the product representation
            # sti: start index of the domain
            # similar for edi, stf, edf (initial and final)
            winfo_init = prod_weight_info[zval_init]
            winfo_final = prod_weight_info[zval_final]
            if verbose > 1 println("Info of $(zval_init) (initial weight) is $(winfo_init)") end
            if verbose > 1 println("Info of $(zval_final) (final weight) is $(winfo_final)") end
            for ((z1, z2), (sti, edi)) in winfo_init
                z1c = Tuple(z1 .- dzs[opi])
                if haskey(winfo_final, (z1c, z2))
                    if verbose > 1 println("From ($(z1), $(z2)) to ($(z1c), $(z2))") end
                    if verbose > 1 println("Lowering in the first space") end
                    stf, edf = winfo_final[(z1c, z2)]
                    r1_lop = r1.Sl[opi][z1]
                    if verbose > 2 println("Lowering operator of the first space is $(r1_lop)") end
                    @assert (edf-stf+1) % r1_lop.m == 0
                    idsize = div(edf-stf+1, r1_lop.m)
                    @assert r1_lop.n * idsize == edi-sti+1
                    mat[stf:edf, sti:edi] = r1_lop ⊗ sparse(I, idsize, idsize)
                end

                z2c = Tuple(z2 .- dzs[opi])
                if haskey(winfo_final, (z1, z2c))
                    if verbose > 1 println("From ($(z1), $(z2)) to ($(z1), $(z2c))") end
                    if verbose > 1 println("Lowering in the second space") end
                    stf, edf = winfo_final[(z1, z2c)]
                    r2_lop = r2.Sl[opi][z2]
                    if verbose > 2 println("Lowering operator of the second space is $(r2_lop)") end
                    @assert (edf-stf+1) % r2_lop.m == 0
                    idsize = div(edf-stf+1, r2_lop.m)
                    @assert r2_lop.n * idsize == edi-sti+1
                    mat[stf:edf, sti:edi] = sparse(I, idsize, idsize) ⊗ r2_lop
                end
            end

            nSl[opi][zval_init] = mat
        end
    end

    return nSl, nSz, prod_inmult, prod_weight_info
end

cgt_outputmat(irep::Irep{S}, 
::Type{FT}) where {S<:AbelianSymm, FT<:AbstractFloat} = SparseArray([FT(1);;])

# Tmat: Linear trasnformation encoding the orthogonalization of basis vectors
# Nvec: Norm^2 of the orthogonalized basis vectors
# Used when getting CGC in sparse float tensor form
function cgt_outputmat(irep::Irep,
    ::Type{FT}) where {FT<:AbstractFloat}
    dim = dimension(irep)
    mat = SparseArray(zeros(FT, dim, dim))
    for (z, (a, b)) in irep.Sz
        tmat, nvec = get_TNmat(irep.innerprod[z]) 
        nvec_sqinv = FT(1) ./ sqrt.(Vector{FT}(nvec))
        mat[a:b, a:b] = Matrix{FT}(tmat) * Diagonal(nvec_sqinv)
    end
    return mat
end

function cgt_inputmat(irep::Irep,
    ::Type{FT}) where {FT<:AbstractFloat}
    dim = dimension(irep)
    mat = SparseArray(zeros(FT, dim, dim))
    for (z, (a, b)) in irep.Sz
        tmat, nvec = get_TNmat(irep.innerprod[z]) 
        nvec_sqrt = sqrt.(Vector{FT}(nvec))
        mat[a:b, a:b] = transpose(Diagonal(nvec_sqrt) * inv(Matrix{FT}(tmat)))
    end
    return mat
end

@generated function contract_ith(cgt_float::SparseArray{FT, N}, 
    mat::SparseArray{FT}, 
    ::Val{I}) where {N, I, FT<:AbstractFloat}
    @assert I < N
    inds = [Symbol(:i, i) for i in 1:N]
    out_inds = copy(inds)
    out_inds[I] = :a  # Replace ith index with 'a'
    
    contraction = :(mat[$(inds[I]), a])
    
    quote
        @tensor cgt[$(out_inds...)] := cgt_float[$(inds...)] * $contraction
        return cgt::SparseArray{FT}
    end
end

	
# Convert a CGT object into a float sparse tensor
function to_float(cgt::CGT{S, CT, NZ, N},
    ::Type{FT},
    normalize=true) where {S<:NonabelianSymm, CT<:Number, NZ, N, FT<:AbstractFloat}

    qlabels = cgt.qlabels
    om = length(cgt.nfactor) # outer multiplicity
    ireps = [getNsave_irep(S, BigInt, q) for q in qlabels]
    cgt_size = [dimension(irep) for irep in ireps]
    cgt_float = SparseArray{FT}(zeros(cgt_size..., om))

    Ncolons = [Colon() for _ in 1:N]
    for (k, block) in cgt.blocks
        ranges = UnitRange[]
        for i=1:N
            ithqlabel = k[i]
            stidx, edidx = ireps[i].Sz[ithqlabel]
            push!(ranges, stidx:edidx)
        end
        for j=1:om
            cgt_float[ranges..., j] = Array{FT}(block[Ncolons..., j])
        end
    end

    if normalize
        for i=1:om cgt_float[Ncolons..., i] .*= sqrt(FT(cgt.nfactor[i])) end
    end

    for i=1:N
        dir = cgt.dir[i]
        if dir == '-'
            mat = cgt_outputmat(ireps[i], FT)
        elseif dir == '+'
            mat = cgt_inputmat(ireps[i], FT)
        else
            error("Invalid direction: $dir")
        end
        # Just contract ith index of cgt_float and 1st inex of mat
        cgt_float = contract_ith(cgt_float, mat, Val(i))
    end

    return cgt_float, qlabels, cgt.dir
end

# Load (or generate if needed) a CGT and convert it into a float sparse tensor
function load_cg3_float(::Type{S},
    ::Type{CT},
    qlabels_in::NTuple{2, NTuple{NZ, Int}},
    qlabel_out::NTuple{NZ, Int},
    ::Type{FT},
    normalize=false) where {S<:Symmetry, CT<:Number, NZ, FT<:AbstractFloat}

    @assert NZ == nzops(S)
    dir = ['+', '+', '-']
    qlabels_in_sorted = Tuple(minmax(qlabels_in[1], qlabels_in[2]))
    permute = qlabels_in[1] > qlabels_in[2]

    if S<:NonabelianSymm
        cg3s = getNsave_cg3(S, CT, qlabels_in_sorted, [qlabel_out])
        cgt_float, qlabels, dir = to_float(cg3s[qlabel_out], FT, normalize)
    else
        @assert qlabel_out[1] == add_qn(S, qlabels_in[1][1], qlabels_in[2][1])
        cgt_float = SparseArray(ones(FT, 1, 1, 1, 1))
        qlabels = (qlabels_in_sorted..., qlabel_out)
    end
    if permute
        cgt_float = permutedims(cgt_float, (2, 1, 3, 4))
        qlabels = (qlabels[2], qlabels[1], qlabels[3])
    end
    return cgt_float, qlabels, dir
end
