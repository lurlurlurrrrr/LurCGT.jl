isvalidsymm(::Type{SU{N}}) where N = N >= 2

isSU2(::Type{SU{2}}) = true
isSU2(::Any) = false

isSUN(::Type{<:SU}) = true
isSUN(::Any) = false

# text expression of the symmetry. Needed to construct file/folder name
totxt(::Type{SU{N}}) where N = "SU$(N)"

# dimension of defining irep
defirepdim(::Type{SU{N}}) where N = N

# number of raising operators (for SU(N), it is N - 1)
nlops(::Type{SU{N}}) where N = N - 1

# number of z-operators (for SU(N), it is N - 1)
nzops(::Type{SU{N}}) where N = N - 1

getsl_triv(::Type{SU{N}}, ::Type{RT}) where {N, RT<:Number} =
Tuple(Dict{NTuple{N-1, Int}, SparseMatrixCSC{RT}}() for _=1:N-1)

getsz_triv(::Type{SU{N}}) where N = Dict(Tuple(0 for _=1:N-1)=>(1, 1))

function getsl_def(::Type{SU{N}}, ::Type{RT}) where {N, RT<:Number} 
	def_sl = Tuple(Dict{NTuple{N-1, Int}, SparseMatrixCSC{RT}}() for _=1:N-1)
	def_sz = getsz_def(SU{N})
	def_zvals = sorted_zvals(def_sz)
	for i=1:N-1
		sz = def_zvals[i]
		val = spzeros(RT, 1, 1)
		val[1, 1] = RT(1) 
		def_sl[i][sz] = val
	end
	return def_sl
end

function getsz_def(::Type{SU{N}}, i) where N 
	sz = zeros(Int, N - 1)
	for j in 1:N-1
		if j >= i
			sz[j] = 1
		elseif j == i - 1
			sz[j] = -i + 1
		else
			sz[j] = 0
		end
	end
	return Tuple(sz)
end

getsz_def(::Type{SU{N}}) where N = Dict(Tuple(getsz_def(SU{N}, i)) => (i, i) for i=1:N)
getsz_def_vec(::Type{SU{N}}) where N = [collect(getsz_def(SU{N}, i)) for i=1:N]

# Change of z-values when lowering operator is applied.
function getdz(::Type{SU{N}}, lop::Int) where N
	sz = getsz_def(SU{N})
	zvals = sorted_zvals(sz)
	return zvals[lop] .- zvals[lop+1]
end

crystal_chars_map(::Type{SU{N}}) where N = Dict(i => i for i=1:N)

# This function should be modified
mw_column(::Type{<:SU}, l::Int) = collect(l:-1:1)

getdzs(::Type{SU{N}}) where N = [collect(getdz(SU{N}, i)) for i=1:nlops(SU{N})]

charlist(::Type{SU{N}}) where N = collect(1:N)

function get_fops_std(::Type{SU{N}}) where N
	fops = Vector{Dict{Int, Int}}()
	for i in 1:N-1
		fop = Dict{Int, Int}(); fop[i] = i + 1;
		push!(fops, fop)
	end
	return fops
end

function qlab2mwz(::Type{SU{N}}, qlabel::NTuple{NZ, Int}) where {N, NZ}
	@assert N - 1 == NZ
	sz_defs = [collect(getsz_def(SU{N}, i)) for i in 1:NZ]
	z = zeros(Int, NZ)
	z_added = zeros(Int, NZ)
	for i in 1:NZ
		z_added .+= sz_defs[i]
		z .+= z_added .* qlabel[i]
	end
	return Tuple(z)
end

function getqlabel(::Type{SU{N}}, z::NTuple{NZ, Int}) where {N, NZ}
	@assert N - 1 == NZ
	w = zeros(Int, NZ)
	for i=NZ:-1:2
        @assert (z[i] - z[i-1]) % i == 0
		w[i] = div(z[i] - z[i-1], i)
	end
	w[1] = z[1]
	return Tuple(w)
end

maxalt(::Type{SU{N}}) where N = N - 1

ytnrows(::Type{SU{N}}) where N = N

# Weight comparision, inputs are a form of shape of the Young tableau
function less_weight(::Type{SU{N}}, w1::NTuple{NZ, Int}, w2::NTuple{NZ, Int}) where {N, NZ}
	@assert NZ == N - 1
	return reverse(w1) < reverse(w2)
end

function fundamental_qlabels(::Type{SU{N}}) where N
	funda_lst = Vector{NTuple{N-1, Int}}()
	for i in 1:N-1
		q = Tuple(i == j ? 1 : 0 for j=1:N-1)
		push!(funda_lst, q)
	end
	return funda_lst
end