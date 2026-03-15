# NI: Number of incoming spaces, NO: Number of outgoing spaces
abstract type AbstractCG3contract{S, NI, NO, NZ} end

# Linear combination of fusion trees in TensorKit.jl
# NI is the number of incoming spaces, nonzero positive integer
struct FTree_{S, NI, NO, NZ, M1, M2} <: AbstractCG3contract{S, NI, NO, NZ}
    # Incoming spaces. No need to be sorted
    ins::NTuple{NI, NTuple{NZ, Int}}
    # Outgoing spaces. 
    outs::NTuple{NO, NTuple{NZ, Int}}
    # Dictionary of components. The key is the intermediate spaces
    # and the value is the coefficient of each (non-normalized) fusion tree
    coeff::Dict{NTuple{M2, NTuple{NZ, Int}}, Array{BigInt, M1}}
    coeff_nfac::Dict{NTuple{M2, NTuple{NZ, Int}}, Rational{BigInt}}
end

const FTree{S, N, NZ, M1, M2} = FTree_{S, N, 1, NZ, M1, M2}

Base.copy(ftree::FTree{S}) where S<:NonabelianSymm =
    create_FTree(S, ftree.ins, ftree.outs, copy(ftree.coeff), copy(ftree.coeff_nfac), false)

create_CG3cont(cont::FTree{S, NI},
    new_ins::NTuple{NI, NTuple{NZ, Int}},
    new_outs::NTuple{1, NTuple{NZ, Int}},
    new_coeffs::Dict{NTuple{M2, NTuple{NZ, Int}}, Array{BigInt, M1}},
    new_coeff_nfac::Dict{NTuple{M2, NTuple{NZ, Int}}, Rational{BigInt}},
    sortcheck::Bool;
    verbose=0) where {S<:NonabelianSymm, NI, NZ, M1, M2} =
    create_FTree(S, new_ins, new_outs, new_coeffs, new_coeff_nfac, sortcheck)


function create_unit_FTree(omlist::OMList{S, N, NZ},
    i::Int) where {S<:NonabelianSymm, N, NZ}

    @assert NZ == nzops(S); @assert N >= 1
    @assert i >= 1 && i <= omlist.totalOM

    KL = max(N - 2, 0)
    coeff = Dict{NTuple{KL, NTuple{NZ, Int}}, Array{BigInt, N-1}}()
    coeff_nfac = Dict{NTuple{KL, NTuple{NZ, Int}}, Rational{BigInt}}()
    @assert omlist.incom_spaces == Tuple(sort(collect(omlist.incom_spaces)))

    intsp_idx = searchsortedlast(omlist.cumul, i)
    nz_interm_space = omlist.interm_spaces[intsp_idx]

    for (ii, interm_space) in enumerate(omlist.interm_spaces)
        coeff[interm_space] = zeros(BigInt, omlist.cg3_oms[ii])
        coeff_nfac[interm_space] = 1
        if (interm_space == nz_interm_space)
            coeff[interm_space][i-omlist.cumul[intsp_idx]+1] = BigInt(1)
        end
    end
    return create_FTree(S, omlist.incom_spaces,
        (omlist.out_space,), coeff, coeff_nfac, true)
end

function create_FTree(::Type{S}, 
    ins::NTuple{N, NTuple{NZ, Int}},
    outs::NTuple{1, NTuple{NZ, Int}},
    coeff::Dict{NTuple{M2, NTuple{NZ, Int}}, Array{BigInt, M1}},
    coeff_nfac::Dict{NTuple{M2, NTuple{NZ, Int}}, Rational{BigInt}},
    sortcheck::Bool;
    verbose=0) where {S<:NonabelianSymm, N, NZ, M1, M2}
    @assert NZ == nzops(S); @assert N >= 1
    @assert M1 == max(N - 1, 0); @assert M2 == max(N - 2, 0)

    # If we need to ensure that the incoming spaces are sorted
    sorted = ins == Tuple(sort(collect(ins)))

    # If the incoming spaces are sorted and sortcheck is true, check 
    if sortcheck
        @assert sorted
        omlist = getNsave_omlist(S, ins, outs[1])
        # It should exist
        @assert !isnothing(omlist)
        
        # Set of key of the dictionary coeff
        keys_set = Set(keys(coeff))
        for (i, interm_space) in enumerate(omlist.interm_spaces)
            if verbose > 1 println("Processing intermediate space: ", interm_space) end
            if interm_space in keys_set
                sz = size(coeff[interm_space])
                if verbose > 1 println(sz) end
                if verbose > 1 println(omlist.cg3_oms[i]) end
                @assert sz == omlist.cg3_oms[i]
                delete!(keys_set, interm_space)
            end
        end
        # If we have removed all keys, it means that all intermediate spaces
        # are accounted for
        @assert isempty(keys_set) 
    end
    return FTree{S, N, NZ, M1, M2}(ins, outs, coeff, coeff_nfac)
end

# Generate a random linear combination of fusion trees
# This function is used for testing purposes
function random_FTree(::Type{S}, 
    ins::NTuple{N, NTuple{NZ, Int}},
    outs::NTuple{1, NTuple{NZ, Int}};
    onearr=false) where {S<:NonabelianSymm, N, NZ}

    KL = max(N - 2, 0)
    coeff = Dict{NTuple{KL, NTuple{NZ, Int}}, Array{BigInt, N-1}}()
    coeff_nfac = Dict{NTuple{KL, NTuple{NZ, Int}}, Rational{BigInt}}()
    @assert ins == Tuple(sort(collect(ins)))

    omlist = getNsave_omlist(S, ins, outs[1])
    @assert !isnothing(omlist)
    for (i, interm_space) in enumerate(omlist.interm_spaces)
        coeff[interm_space] = 
        onearr ? ones(Int8, omlist.cg3_oms[i]) : rand(Int8, omlist.cg3_oms[i])
        coeff_nfac[interm_space] = 1
    end

    return create_FTree(S, ins, outs, coeff, coeff_nfac, true)
end

# Decompose a permutation into adjacent transpositions
function decompose_perm(perm::Tuple)
    len = length(perm)
    vec = collect(1:len)
    switch_list = Vector{Tuple{Int, Int}}()

    for i=1:len
        target_elem = perm[i]
        start_idx = findfirst(x -> x == target_elem, vec)
        @assert start_idx !== nothing "Invalid permutation: $perm"
        if start_idx == i continue end  # Already in place
        # Swap the elements
        for k=start_idx-1:-1:i
            vec[k], vec[k+1] = vec[k+1], vec[k]
            push!(switch_list, (k, k+1))
        end
        if Tuple(vec) == Tuple(perm) return switch_list end
    end
end

# FTree in the argument is modified in place
# However, the resulting FTree is different from the original one
function Base.permute!(FTree::FTree{S, N}, perm::NTuple{N, Int};
    verbose=0) where {S<:NonabelianSymm, N}

    # If the permutation is already identity, return the FTree
    if perm == Tuple(1:N)
        if verbose > 1 println("Permutation is identity, returning original FTree") end
        return FTree
    end
    # Decompose the permutation
    switch_list = decompose_perm(perm)
    if verbose > 1 println("Switch list: $(switch_list)") end
    for (i, j) in switch_list
        if verbose > 2 println("Permuting indices $(i), $(j)") end
        @assert i >= 1 && j <= N && j == i + 1
        FTree = permute_adjacent!(FTree, i, j; verbose=verbose)
    end
    return FTree
end

function permute_adjacent!(FTree::FTree{S, N}, i1, i2;
    verbose=0) where {S<:NonabelianSymm, N}
    @assert i2 == i1 + 1
    @assert i1 >= 1 && i2 <= N
    new_ins = (FTree.ins[1:i1-1]..., 
        FTree.ins[i2], FTree.ins[i1], 
        FTree.ins[i2+1:end]...)
    
    if verbose > 1 println("New indices: $(new_ins)") end
    # This is the case when only one R-symbol is needed
    if i1 == 1
        if FTree.ins[1] == FTree.ins[2]
            apply_Rsymbol!(FTree, N-1, rsymiofunc_perm_step1; verbose)
        end
        return create_FTree(S, new_ins, FTree.outs, FTree.coeff, FTree.coeff_nfac, false)
    end

    if verbose > 1 println("R") end
    apply_Rsymbol!(FTree, N - i1, rsymiofunc_perm_step1; verbose)
    if verbose > 1 println("F") end
    FTree_new = apply_Fsymbol(FTree, new_ins, FTree.outs, N - i1, true, fsymiofunc_perm_step1; verbose)
    if verbose > 1 println("R") end
    apply_Rsymbol!(FTree_new, N - i1 + 1, rsymiofunc_perm_step1; verbose)
    return FTree_new
end

function fsymiofunc_perm_step1(FTree::FTree{S, N}, 
    rem::NTuple{M, NTuple{NZ, Int}}, 
    i) where {S<:NonabelianSymm, N, M, NZ}
    outsp = i==1 ? FTree.outs[1] : rem[i-1]
    in2 = i==N-2 ? FTree.ins[1] : rem[i]
    in1, in3 = FTree.ins[N-i+1], FTree.ins[N-i]
    return in1, in2, in3, outsp
end

# rsym input/output function used in permutation step of step 1
function rsymiofunc_perm_step1(FTree::FTree{S, NI},
    intermsp::NTuple{M1, NTuple{NZ, Int}},
    i::Int) where {S<:NonabelianSymm, NI, M1, NZ}

    in1 = i < NI-1 ? intermsp[i] : FTree.ins[1]
    in2 = FTree.ins[NI-i+1]
    rout = i==1 ? FTree.outs[1] : intermsp[i-1]
    return in1, in2, rout
end


# Apply R-symbol to the i-th CG3 in AbstractCG3contract type
function apply_Rsymbol!(CG3cont::AbstractCG3contract{S, NI, NO, NZ},
    i::Int,
    rsymiofunc::Function; 
    verbose=0) where {S<:NonabelianSymm, NI, NO, NZ}

    @assert CG3cont.ins isa NTuple{NI, NTuple{NZ, Int}}
    @assert CG3cont.outs isa NTuple{NO, NTuple{NZ, Int}}
    coeff, coeff_nfac = CG3cont.coeff, CG3cont.coeff_nfac
    if verbose > 1 println("Applying R-symbol to the $(i)th CG3") end
    for intermsp in keys(coeff)
        if verbose > 1 println("Processing intermediate space: $(intermsp)") end
        old_coeff, old_nfac = coeff[intermsp], coeff_nfac[intermsp]
        in1, in2, rout = rsymiofunc(CG3cont, intermsp, i)
        # Nontrivial R-symbol appears only when in1 == in2
        if in1 == in2
            # Prepare the corresponding R-symbol
            if verbose > 1 println("Nontrivial R-symbol: $(in1), $(rout)") end
            rsym = getNsave_Rsymbol(S, BigInt, in1, rout)
            rsym_mat, nfac = rsym_rightnormalized(rsym)

            new_arr = contract_ith(old_coeff, rsym_mat, Val(i))
            new_nfac = old_nfac // nfac

            new_arr, new_nfac = arr_relprime(new_arr, new_nfac)
            coeff[intermsp] = new_arr
            coeff_nfac[intermsp] = new_nfac
        end
    end
end

function arr_relprime(arr, nfac::Rational)
    cfac = gcd(arr)
    n = gcd(cfac, nfac.den)
    return div.(arr, n), nfac * n
end

check_irange(::Type{<:FTree{S, N}}, i::Int) where {S<:NonabelianSymm, N} = 
    @assert i >= 1 && i <= N - 2

# Key length, Array dimension of FTree type
get_dict_param(::Type{<:FTree{S, N}}) where {S<:NonabelianSymm, N} = (N - 2, N - 1) 

# Apply F-symbol to change the ith intermediate space
# For fusion tree, it is applied to the i-th and (i+1)-th CG3s of the fusion tree
function apply_Fsymbol(CG3cont::AbstractCG3contract{S, NI, NO, NZ}, 
    new_ins::NTuple{NI, NTuple{NZ, Int}},
    new_outs::NTuple{NO, NTuple{NZ, Int}},
    i::Int,
    backward::Bool,
    fsymiofunc::Function;
    verbose=0,
    permute_updown=false) where {S<:NonabelianSymm, NZ, NI, NO}

    check_irange(typeof(CG3cont), i)
    KL, AD = get_dict_param(typeof(CG3cont)) # Key length, Array dimension

    new_coeff = Dict{NTuple{KL, NTuple{NZ, Int}}, Array{BigInt, AD}}()
    new_coeff_nfac = Dict{NTuple{KL, NTuple{NZ, Int}}, Rational{BigInt}}()

    permute_tuple = (1:i-1..., i+1, i, i+2:AD...)

    ks_collected = Dict{NTuple{KL-1, NTuple{NZ, Int}}, Set{NTuple{NZ, Int}}}()
    for k in keys(CG3cont.coeff)
        intermsp = (k[1:i-1]..., k[i+1:KL]...)
        if !haskey(ks_collected, intermsp)
            ks_collected[intermsp] = Set()
        end
        push!(ks_collected[intermsp], k[i])
    end

    for (rem, ithset) in ks_collected
        if verbose > 1
            println("Processing intermediate space: $(rem), with set: $(ithset)")
        end
        in1, in2, in3, outsp = fsymiofunc(CG3cont, rem, i)

        if verbose > 1
            println("Corresponding F-symbol: $(in1), $(in2), $(in3) -> $(outsp)")
        end
        Fsymbol = getNsave_Fsymbol(S, BigInt, in1, in2, in3, outsp)
        @assert !isnothing(Fsymbol)
        Fsym_sz = size(Fsymbol)
        if verbose > 1 println("F-symbol size: $(Fsym_sz)") end
        
        firstelem = first(ithset); 
        intsps = (rem[1:i-1]..., firstelem, rem[i:KL-1]...)
        if verbose > 1 println("intsps: $(intsps)") end
        om_tuple = size(CG3cont.coeff[intsps])
        if verbose > 1 println("om_tuple: $(om_tuple)") end
        
        arr_size = (om_tuple[1:i-1]..., Fsym_sz, om_tuple[i+2:AD]...)
        if verbose > 1 println("arr_size: $(arr_size)") end
        arr_before = zeros(BigInt, arr_size)

        intsplst_before = backward ? Fsymbol.fs_list : Fsymbol.es_list
        comm_fac = 0
        for (intsp, _) in intsplst_before
            if intsp in ithset
                intsps = (rem[1:i-1]..., intsp, rem[i:KL-1]...)
                nfac = CG3cont.coeff_nfac[intsps].den
                comm_fac = comm_fac==0 ? nfac : lcm(comm_fac, nfac)
            end
        end
        if verbose > 1 println("comm_fac: $(comm_fac)") end

        # Fill the array before applying the F-symbol
        ii = 1
        for (intsp, (om1, om2)) in intsplst_before
            if intsp in ithset
                if verbose > 1 println("intsp = $(intsp)") end
                intsps = (rem[1:i-1]..., intsp, rem[i:KL-1]...)
                iirange = ii:(ii + om1 * om2 - 1)
                
                idx = ((Colon() for _=1:i-1)..., iirange, (Colon() for _=i+2:AD)...)
                reshape_size = (om_tuple[1:i-1]..., om1*om2, om_tuple[i+2:AD]...)

                @assert comm_fac % CG3cont.coeff_nfac[intsps].den == 0
                mfac = div(comm_fac, CG3cont.coeff_nfac[intsps].den)
                if verbose > 1 println("mfac: $(mfac)\n") end
                coeff_orig = permute_updown ? permutedims(CG3cont.coeff[intsps], permute_tuple) :
                    CG3cont.coeff[intsps]
                arr_before[idx...] = reshape(coeff_orig, reshape_size) * mfac
                delete!(ithset, intsp)
            end
            ii += om1 * om2
        end
        @assert isempty(ithset)
        if verbose > 1 println("Before:"); display(arr_before) end

        # Apply the F-symbol
        fsym_mat, fsym_nfac::Integer = 
            backward ? fsym_leftnormalized(Fsymbol) : fsym_rightnormalized(Fsymbol)
        arr_after = contract_ith(arr_before, backward ? fsym_mat' : fsym_mat, Val(i))
        nfac_after::Integer = comm_fac * fsym_nfac

        if verbose > 1 println("After:"); display(arr_after) end
        if verbose > 1 println("nfac_after: $(nfac_after)") end

        # Save the resulting array in the resulting CG3cont
        ii = 1
        intsplst_after = backward ? Fsymbol.es_list : Fsymbol.fs_list
        for (intsp, (om1, om2)) in intsplst_after
            if verbose > 1 println("intsp = $(intsp)") end
            if verbose > 1 println("om1, om2 = $(om1), $(om2)") end

            iirange = ii:(ii + om1 * om2 - 1)
            idx = ((Colon() for _=1:i-1)..., iirange, (Colon() for _=i+2:AD)...)
            part_arr = arr_after[idx...]
            if !iszero(part_arr)
                intsps = (rem[1:i-1]..., intsp, rem[i:KL-1]...)
                new_arr, new_nfac::Rational = arr_relprime(part_arr, 1//nfac_after)
                reshape_size = (om_tuple[1:i-1]..., om1, om2, om_tuple[i+2:AD]...)
                new_arr_resh = reshape(new_arr, reshape_size)

                new_coeff[intsps] = permute_updown ?
                    permutedims(new_arr_resh, permute_tuple) : new_arr_resh
                new_coeff_nfac[intsps] = new_nfac
            end
            ii += om1 * om2
        end
        if verbose > 1 println("\n\n") end
    end
    return create_CG3cont(CG3cont, new_ins, new_outs, new_coeff, new_coeff_nfac, false; verbose)
end
