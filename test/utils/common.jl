⊗(a::Matrix, b::Matrix) = kron(b, a)
⊗(a::Vector, b::Vector) = kron(b, a)

comm(A, B) = A * B - B * A 

"""
    lowering_matrix(rep::LurCGT.Irep, i::Int, FT=Float64)

Return a `SparseMatrixCSC{FT,Int}` in an orthonormal basis, so its ordinary
adjoint represents the metric adjoint of the stored lowering operator.
For each weight factor the full Gram block as `G = C' * C` using upper
Cholesky factors, and transform each block as `C_dest * F * inv(C_source)`
using a triangular solve. Weight ranges and `ordered_zvalues` remain valid,
but basis vectors within each range change. No input is mutated.
`FT` must be an abstract-float subtype; conversion is approximate. Invalid
indices raise `BoundsError`, and non-positive-definite metrics raise an error.
"""
function lowering_matrix(rep::LurCGT.Irep{S}, i::Int,
        ::Type{FT}=Float64) where {S,FT<:AbstractFloat}
    blocks = rep.Sl[i]
    factors = Dict(w => cholesky(Symmetric(FT.(G))).U for (w,G) in rep.innerprod)
    rows = Int[]; cols = Int[]; values = FT[]
    dz = LurCGT.getdz(S,i)
    for (w,block) in blocks
        dest = w .- dz
        transformed = factors[dest] * Matrix{FT}(block) / factors[w]
        r,c,v = findnz(sparse(transformed))
        append!(rows,r .+ (rep.Sz[dest][1]-1))
        append!(cols,c .+ (rep.Sz[w][1]-1))
        append!(values,v)
    end
    n = LurCGT.dimension(rep)
    return sparse(rows,cols,values,n,n)
end

"""
    ordered_zvalues(rep::LurCGT.Irep) -> Vector{NTuple{NZ,Int}}

Return one stored z-weight per basis state in the global order given by `rep.Sz`,
matching the rows and columns of `lowering_matrix`. Degenerate weights repeat
for every state in their inclusive index range. The returned vector is new;
`rep` is not mutated and weights are not converted to Dynkin labels.
"""
function ordered_zvalues(rep::LurCGT.Irep{S,NL,NZ}) where {S,NL,NZ}
    weights = Vector{NTuple{NZ,Int}}(undef,LurCGT.dimension(rep))
    for (w,(first,last)) in rep.Sz
        for k in first:last
            weights[k] = w
        end
    end
    return weights
end

function generate_incom(::Type{S}, dim_limit, qlimit=4) where S<:NonabelianSymm
    NZ = nzops(S)
    ins = Vector{NTuple{NZ, Int}}()
    dim = 1
    while true
        qlabel = Tuple(rand(0:qlimit) for _=1:NZ)
        irep_dim = LurCGT.dimension(LurCGT.getNsave_irep(S, BigInt, qlabel))
        prod_dim = dim * irep_dim
        if prod_dim > dim_limit break end
        push!(ins, qlabel)
        dim = prod_dim
    end
    return ins
end

function select_out(vo::LurCGT.ValidOuts{S, N, NZ},
    outlimit=100) where {S<:NonabelianSymm, N, NZ}
    outlist = vo.out_spaces
    outdims = [dimension(S, out) for out in outlist]
    if minimum(outdims) > outlimit
        i = argmin(outdims)
        return outlist[i]
    end
    while true
        out = rand(vo.out_spaces)
        if dimension(S, out) <= outlimit return out end
    end
end

function select_out(::Type{S},
    vec::Vector{NTuple{NZ, Int}},
    outlimit=100) where {S<:NonabelianSymm, NZ}

    outdims = [dimension(S, out) for out in vec]
    if minimum(outdims) > outlimit
        i = argmin(outdims)
        return vec[i]
    end
    while true
        out = rand(vec)
        if dimension(S, out) <= outlimit return out end
    end
end

@generated function contract_om_step3(arr_up, 
    arr_dn, 
    om_arr, 
    ::Val{NI}, 
    ::Val{NO}) where {NI, NO}

    result_inds = [[Symbol(:in, i) for i=1:NI]..., [Symbol(:out, i) for i=1:NO]...]
    up_inds = [[Symbol(:in, i) for i=1:NI]..., :c, [Symbol(:oup, i) for i=1:NI-1]...]
    dn_inds = [[Symbol(:out, i) for i=1:NO]..., :c, [Symbol(:odn, i) for i=1:NO-1]...]
    om_inds = [[Symbol(:oup, i) for i=1:NI-1]..., [Symbol(:odn, i) for i=1:NO-1]...]

    quote
        @tensor result[$(result_inds...)] :=
            arr_up[$(up_inds...)] * arr_dn[$(dn_inds...)] * om_arr[$(om_inds...)]
        return result
    end
end

function get_spaces(::Type{S}, 
    dim_limit::Int, 
    qlimit::Int, 
    outlimit::Int; 
    verbose) where S<:NonabelianSymm

    NZ = nzops(S)
    # Get input and output of the first CGT
    incom = sort(generate_incom(S, dim_limit, qlimit))
    len = length(incom)
    if len < 3 return nothing end
    # The number of outgoing legs of CGT1
    cgt1_up = rand(1:len-1)
    # The number of incoming legs of CGT1
    cgt1_dn = len - cgt1_up

    shuffle_len = shuffle(1:len)
    
    cgt1_outsp = incom[sort(shuffle_len[1:cgt1_up])]
    cgt1_insp = incom[sort(shuffle_len[cgt1_up+1:end])]

    # If there is no common factor of product of input and product of output, skip
    if cgt1_up > 1
        vo_out = getNsave_validout(S, Tuple(cgt1_outsp))
        possible_out = Set(vo_out.out_spaces)
    else possible_out = Set(cgt1_outsp[1]) end

    if cgt1_dn > 1
        vo_in = getNsave_validout(S, Tuple(cgt1_insp))
        possible_in = Set(vo_in.out_spaces)
    else possible_in = Set(cgt1_insp[1]) end

    intersect_outsps = intersect(possible_out, possible_in)
    if isempty(intersect_outsps) return nothing end
    cgt1_internalsp = rand(collect(intersect_outsps))
    
    # Number of contracted legs
    nctlegs = rand(1:len)
    # The number of contracted output spaces from CGT1
    nctlegs_out1 = rand(max(nctlegs-cgt1_dn,0):min(nctlegs, cgt1_up))
    # The number of contracted input spaces from CGT1
    nctlegs_in1 = nctlegs - nctlegs_out1

    # Pick 'nctlegs_out1' legs from the output spaces of CGT1
    contracted_1out2in = cgt1_outsp[sort(shuffle(1:cgt1_up)[1:nctlegs_out1])]
    # Pick 'nctlegs_in1' legs from the input spaces of CGT1
    contracted_1in2out = cgt1_insp[sort(shuffle(1:cgt1_dn)[1:nctlegs_in1])]

    # List of contracted spaces
    contracted_spaces = vcat(contracted_1out2in, contracted_1in2out)
    # Product of dimensions of contracted spaces
    contracted_dim = prod([LurCGT.dimension(LurCGT.getNsave_irep(S, BigInt, sp))
                            for sp in contracted_spaces])


    # Get the input and output of the second CGT
    # The number of contracted input spaces from CGT2 (=nctlegs_out1)
    nctlegs_in2 = nctlegs_out1
    # The number of contracted output spaces from CGT2 (=nctlegs_in1)
    nctlegs_out2 = nctlegs_in1
    remain_dim = 10 * div(dim_limit, contracted_dim)
    found = false
    find_cnt = 0

    cgt2_insp, cgt2_outsp = [], []
    cgt2_internalsp = nothing

    # Possible minimum length of incomnew2
    len2_min = max(3-nctlegs_in2-nctlegs_out2,
                    nctlegs_in2==0 || nctlegs_out2==0 ? 1 : 0)

    while true
        incomnew2 = sort(generate_incom(S, remain_dim, qlimit))
        len2 = length(incomnew2)
        if len2 < len2_min continue end

        nopenlegs_in2_min = nctlegs_in2>0 ? 0 : 1
        nopenlegs_in2_max = nctlegs_out2>0 ? len2 : len2-1
        # The number of open input legs of CGT2
        nopenlegs_in2 = rand(nopenlegs_in2_min:nopenlegs_in2_max)
        # The number of open output legs of CGT2
        nopenlegs_out2 = len2 - nopenlegs_in2
        @assert nctlegs_in2 + nopenlegs_in2 >= 1 && nctlegs_out2 + nopenlegs_out2 >= 1
        @assert nctlegs_in2 + nopenlegs_in2 + nctlegs_out2 + nopenlegs_out2 >= 3

        shuffle_len2 = shuffle(1:len2)
        # Open input spaces of CGT2
        opensp_in2 = incomnew2[sort(shuffle_len2[1:nopenlegs_in2])]
        # Open output spaces of CGT2
        opensp_out2 = incomnew2[sort(shuffle_len2[nopenlegs_in2+1:end])]

        # Open input, output spaces of CGT2, respectively
        cgt2_insp = sort(vcat(opensp_in2, contracted_1out2in))
        cgt2_outsp = sort(vcat(opensp_out2, contracted_1in2out))

        # Check whether the output spaces are valid
        if length(cgt2_insp) > 1
            vo_in2 = getNsave_validout(S, Tuple(cgt2_insp))
            possible_in = Set(vo_in2.out_spaces)
        else possible_in = Set(cgt2_insp[1]) end

        if length(cgt2_outsp) > 1
            vo_out2 = getNsave_validout(S, Tuple(cgt2_outsp))
            possible_out = Set(vo_out2.out_spaces)
        else possible_out = Set(cgt2_outsp[1]) end

        intersect_outsps = intersect(possible_out, possible_in)
        if !isempty(intersect_outsps) 
            found = true
            cgt2_internalsp = rand(collect(intersect_outsps))
            break 
        end

        if find_cnt > 100 break end
        find_cnt += 1
    end
    if !found return nothing end

    # cgt1_insp, cgt1_outsp, cgt2_insp, cgt2_outsp are now ready
    # also contracted_1in2out and contracted_1out2in
    # also cgt1_internalsp and cgt2_internal
    @assert length(cgt1_insp) >= 1 && length(cgt1_outsp) >= 1
    @assert length(cgt2_insp) >= 1 && length(cgt2_outsp) >= 1
    @assert length(cgt1_insp) + length(cgt1_outsp) >= 3
    @assert length(cgt2_insp) + length(cgt2_outsp) >= 3

    cgt1legs = Int[]
    cgt2legs = Int[]
    contracted_1in2out_dict = Dict{NTuple{NZ, Int}, Int}()
    contracted_1out2in_dict = Dict{NTuple{NZ, Int}, Int}()

    for sp in contracted_1in2out
        if verbose > 1 println(sp) end
        if haskey(contracted_1in2out_dict, sp)
            contracted_1in2out_dict[sp] += 1
        else
            contracted_1in2out_dict[sp] = 1
        end
    end

    for sp in contracted_1out2in
        if verbose > 1 println(sp) end
        if haskey(contracted_1out2in_dict, sp)
            contracted_1out2in_dict[sp] += 1
        else
            contracted_1out2in_dict[sp] = 1
        end
    end

    if verbose > 1
        println("contracted_1in2out: $(contracted_1in2out)")
        println("contracted_1out2in: $(contracted_1out2in)")
        println("contracted_1in2out_dict: $(contracted_1in2out_dict)")
        println("contracted_1out2in_dict: $(contracted_1out2in_dict)")
    end

    for (sp, cnt) in contracted_1in2out_dict
        idx1_in = sample(findall(==(sp), cgt1_insp), cnt; replace=false)
        idx2_out = sample(findall(==(sp), cgt2_outsp), cnt; replace=false)
        append!(cgt1legs, idx1_in.+length(cgt1_outsp))
        append!(cgt2legs, idx2_out)
    end
    for (sp, cnt) in contracted_1out2in_dict
        idx1_out = sample(findall(==(sp), cgt1_outsp), cnt; replace=false)
        idx2_in = sample(findall(==(sp), cgt2_insp), cnt; replace=false)
        append!(cgt1legs, idx1_out)
        append!(cgt2legs, idx2_in.+length(cgt2_outsp))
    end

    if verbose > 1
        println(cgt2_insp)
        println(cgt2_outsp)
    end
    return cgt1_outsp, cgt1_insp, cgt2_outsp, cgt2_insp, cgt1_internalsp, cgt2_internalsp, 
    cgt1legs, cgt2legs, contracted_1out2in, contracted_1in2out
end

"""
    check_generator_residual(residual, scale, context; atol, rtol)

Require a finite Frobenius/Euclidean residual norm below `atol + rtol*scale`.
Throw an error with relation context and the measured bound on failure. This
uses an absolute bound for zero-valued relations, including nested commutators.
"""
function check_generator_residual(residual, scale, context; atol, rtol)
    measured = norm(residual)
    bound = atol + rtol*scale
    isfinite(measured) && isfinite(bound) && measured <= bound ||
        error("$context: residual $measured exceeds bound $bound")
    return nothing
end

"""
    validate_crystal_generators(rep, FT=Float64; reference=nothing, atol=100eps(FT), rtol=1000eps(FT))

Validate exact lowering weight shifts and metric-aware sparse self/mixed
commutators and lowering Serre relations. Fit each self-commutator to stored
z-coordinates and anchor its action to `getqlabel`. A full-rank weight matrix
is required unless `reference` supplies coefficients from a faithful irrep of
the same symmetry; this permits checking SO(4) spin fundamentals. When supplied,
compare the fitted coefficients with that reference. Return a matrix with one
coefficient row per lowering index for comparison across irreps. No inputs are mutated. Failures identify the irrep,
operator indices, relation, and residual. Use BigFloat for higher-precision
checks; an irrep with deficient weight rank cannot establish its own reference.
"""
function validate_crystal_generators(rep::LurCGT.Irep{S,NL,NZ},
        ::Type{FT}=Float64; reference=nothing, atol=100eps(FT), rtol=1000eps(FT)) where {S,NL,NZ,FT<:AbstractFloat}
    weights = ordered_zvalues(rep)
    n = length(weights)
    context = "$S $(rep.qlabel)"
    Z = FT[weights[k][a] for k in 1:n, a in 1:NZ]
    full_rank = n >= NZ && all(abs.(diag(qr(Z).R)) .> rtol*norm(Z))
    reference === nothing || size(reference) == (NL,NZ) ||
        throw(ArgumentError("$context: wrong reference coefficient dimensions"))
    full_rank || reference !== nothing ||
        error("$context: insufficient weight rank for generator coefficients")
    matrices = [lowering_matrix(rep,i,FT) for i in 1:NL]
    coefficients = Matrix{FT}(undef,NL,NZ)
    for i in 1:NL
        for (w,block) in rep.Sl[i]
            dest = w .- LurCGT.getdz(S,i)
            haskey(rep.Sz,dest) || error("$context node $i: lowering leaves stored weights at $w")
            rows,cols,values = findnz(block)
            for (r,c,v) in zip(rows,cols,values)
                iszero(v) && continue
                weights[rep.Sz[dest][1]+r-1] == weights[rep.Sz[w][1]+c-1] .- LurCGT.getdz(S,i) ||
                    error("$context node $i: incorrect weight shift at $w entry ($r,$c)")
            end
        end
        m = matrices[i]; H = m'*m-m*m'; h = diag(H)
        scale = 2norm(m)^2
        check_generator_residual(H-spdiagm(0=>h),scale,"$context node $i diagonal self-commutator";atol,rtol)
        c = full_rank ? Z \ h : FT.(reference[i,:])
        if reference !== nothing
            expected_c = FT.(reference[i,:])
            check_generator_residual(c-expected_c,max(norm(c),norm(expected_c)),"$context node $i reference coefficients";atol,rtol)
        end
        check_generator_residual(Z*c-h,max(norm(Z)*norm(c),norm(h)),"$context node $i z-coefficient fit";atol,rtol)
        expected = FT[LurCGT.getqlabel(S,w)[i] for w in weights]
        check_generator_residual(h-expected,max(scale,norm(expected)),"$context node $i coroot normalization";atol,rtol)
        coefficients[i,:] = c
    end
    for i in 1:NL, j in 1:NL
        i == j && continue
        x,y = matrices[i],matrices[j]
        check_generator_residual(x'*y-y*x',2norm(x)*norm(y),"$context nodes ($i,$j) mixed commutator";atol,rtol)
        aij = LurCGT.getqlabel(S,Tuple(LurCGT.getdz(S,j)))[i]
        R = y; scale = norm(y)
        for _ in 1:1-aij
            R = x*R-R*x
            scale *= 2norm(x)
        end
        check_generator_residual(R,scale,"$context nodes ($i,$j) Serre relation";atol,rtol)
    end
    return coefficients
end
