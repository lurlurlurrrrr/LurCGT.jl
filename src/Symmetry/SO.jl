# For SO(N), the symmetries are greatly different for even and odd N
# Start from odd N first
isvalidsymm(::Type{SO{N}}) where N = 4 <= N

isSON(::Type{<:SO}) = true
isSON(::Any) = false

# dimension of defining irep
defirepdim(::Type{SO{N}}) where N = N

# Text expression of the symmetry.
totxt(::Type{SO{N}}) where N = "SO$(N)"

# number of raising operators (for SO(N), it is floor(N / 2))
nlops(::Type{SO{N}}) where N = div(N, 2)

# number of z-operators (for SO(N), it is floor(N / 2))
nzops(::Type{SO{N}}) where N = div(N, 2)

maxalt(::Type{SO{N}}) where N = div(N, 2)

getsl_triv(::Type{SO{N}}, ::Type{RT}) where {N, RT<:Number} =
Tuple(Dict{NTuple{div(N, 2), Int}, SparseMatrixCSC{RT}}() for _=1:div(N, 2))

getsz_triv(::Type{SO{N}}) where N = Dict(Tuple(0 for _=1:div(N, 2))=>(1, 1))

function getsl_def(::Type{SO{N}}, ::Type{RT}) where {N, RT<:Number} 
    NZ = div(N, 2); Nodd = (N % 2 == 1)
    def_sl = Tuple(Dict{NTuple{div(N, 2), Int}, SparseMatrixCSC{RT}}() for _=1:div(N, 2))
    def_sz = getsz_def(SO{N})
    def_zvals = sorted_zvals(def_sz)
    for i=1:NZ
        dz = getdz(SO{N}, i)
        val = spzeros(RT, 1, 1); val[1, 1] = RT(1)
        sw = def_zvals[!Nodd&&i==NZ ? i-1 : i]
        def_sl[i][sw] = val
        def_sl[i][.-(sw.-dz)] = val
    end
    return def_sl
end

# We need to consider spin representation, so weights are multiplied by 2
function getsz_def_vec_(::Type{SO{N}}) where N 
    szs = Vector{NTuple{div(N, 2), Int}}()
    for i in 1:N
        sz = zeros(Int, div(N, 2))
        if N % 2 == 1 && i == N push!(szs, Tuple(sz)); continue; end
        dv, rem = divrem(i+1, 2)
        sz[dv] = rem == 0 ? 2 : -2
        push!(szs, Tuple(sz))
    end
    sort!(szs; lt=rev_less, rev=true)
    return szs
end

getsz_def_vec(::Type{SO{N}}) where N = [collect(t) for t in getsz_def_vec_(SO{N})]

getsz_def(::Type{SO{N}}) where N = Dict(w => (i, i) for (i, w) in enumerate(getsz_def_vec_(SO{N})))


function getdz(::Type{SO{N}}, lop::Int) where N
    NZ = div(N, 2); Nodd = (N % 2 == 1)
    sz_lst = getsz_def_vec_(SO{N})
    if Nodd || lop < NZ return sz_lst[lop] .- sz_lst[lop+1] end
    return sz_lst[lop-1] .- sz_lst[lop+1]
end

function crystal_chars_map(::Type{SO{N}}) where N
    charmap = Dict{Int, Int}(); 
    NZ = div(N, 2); Nodd = N % 2
    for i in 1:NZ charmap[i] = i end
    for i in NZ+1:NZ*2 charmap[i-2*NZ-1] = i + Nodd end
    if Nodd == 1 charmap[0] = NZ+1 end
    return charmap
end

# This function should be modified
function mw_column(::Type{<:SO{N}}, l::Int, aux::Bool) where N
    col = collect(l:-1:1)
    if aux && l == div(N, 2) col[1] = -col[1] end
    return col
end


getdzs(::Type{SO{N}}) where N = [collect(getdz(SO{N}, i)) for i=1:nlops(SO{N})]

function charlist(::Type{SO{N}}) where N
    lst = vcat(collect(1:div(N, 2)), collect(-div(N, 2):-1))
    if N % 2 == 1 push!(lst, 0) end
    return lst
end

# f-operations defined for crystal of Tableau
function get_fops_std(::Type{SO{N}}) where N
    fops = Vector{Dict{Int, Int}}()
    NZ = div(N, 2); Nodd = (N % 2 == 1)
    for i in 1:NZ-1
        fop = Dict{Int, Int}()
        fop[i] = i + 1; fop[-i-1] = -i
        push!(fops, fop)
    end
    fop_last = Dict{Int, Int}()
    if Nodd fop_last[NZ] = 0; fop_last[0] = -NZ
    else fop_last[NZ-1] = -NZ; fop_last[NZ] = -NZ+1 end
    push!(fops, fop_last)
    return fops
end

function qlab2mwz(::Type{SO{N}}, qlabel::NTuple{NZ, Int}) where {N, NZ}
    @assert div(N, 2) == NZ
    if N % 2 == 1
        z = fill(qlabel[NZ], NZ)
        for i in 1:NZ-1 z[NZ+1-i] += 2*sum(qlabel[i:NZ-1]) end
    else
        z = fill(sum(qlabel[NZ-1:NZ]), NZ)
        z[1] -= 2 * qlabel[NZ-1]
        for i in 1:NZ-2 z[NZ+1-i] += 2*sum(qlabel[i:NZ-2]) end
    end
    return Tuple(z)
end

function getqlabel(::Type{SO{N}}, z::NTuple{NZ, Int}) where {N, NZ}
    Nodd = (N % 2 == 1)
    @assert NZ == nzops(SO{N})
    qlabel = zeros(Int, NZ)
    for i in 1:NZ-2
        wdiff = z[NZ+1-i] - z[NZ-i]
        @assert wdiff % 2 == 0
        qlabel[i] = div(wdiff, 2)
    end
    qlabel[NZ-1:NZ] = determine_last(SO{N}, z[1], z[2])
    return Tuple(qlabel)
end

function determine_last(::Type{SO{N}}, z1::Int, z2::Int) where N
    Nodd = (N % 2 == 1)
    if Nodd
        a2 = z2 - z1; b = z1
        @assert a2 % 2 == 0; a = div(a2, 2)
    else
        a2, b2 = z1 + z2, z2 - z1
        @assert a2 % 2 == 0 && b2 % 2 == 0
        a = div(a2, 2); b = div(b2, 2)
    end
    return [a, b]
end

function preprocess(::Type{SO{N}}, qlabel::Vector{Int}) where N
    Nodd = (N % 2 == 1); NZ = div(N, 2)
    @assert length(qlabel) == NZ
    res = copy(qlabel)
    if Nodd @assert qlabel[NZ] % 2 == 0; res[NZ] = div(qlabel[NZ], 2) 
    else 
        @assert (qlabel[NZ-1] + qlabel[NZ]) % 2 == 0
        mi, ma = minmax(qlabel[NZ-1], qlabel[NZ])
        res[NZ-1] = mi; res[NZ] = div(ma - mi, 2)
    end
    return res
end

get_auxarg(::Type{SO{N}}, qlabel::NTuple{NZ, Int}) where {N, NZ} = 
    N % 2 == 0 && qlabel[NZ-1] < qlabel[NZ]

# Weight comparision, inputs are a form of shape of the Young tableau
function less_weight(::Type{SO{N}}, w1::NTuple{NZ, Int}, w2::NTuple{NZ, Int}) where {N, NZ}
	@assert NZ == div(N, 2)
	return reverse(w1) < reverse(w2)
end