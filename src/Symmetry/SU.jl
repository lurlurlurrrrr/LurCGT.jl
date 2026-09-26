"""
    isvalidsymm(::Type{SU{N}}) -> Bool

Return whether `SU{N}` is supported by the generic symmetry machinery.

`N` is the defining dimension. Supported SU groups have `N >= 2`, rank
`N - 1`, and defining representation dimension `N`.
"""
isvalidsymm(::Type{SU{N}}) where N = N >= 2

isSU2(::Type{SU{2}}) = true
isSU2(::Any) = false

isSUN(::Type{<:SU}) = true
isSUN(::Any) = false

# text expression of the symmetry. Needed to construct file/folder name
totxt(::Type{SU{N}}) where N = "SU$(N)"

"""
    defirepdim(::Type{SU{N}}) -> Int

Return the dimension of the defining vector representation.
"""
# dimension of defining irep
defirepdim(::Type{SU{N}}) where N = N

"""
    nlops(::Type{SU{N}}) -> Int
    nzops(::Type{SU{N}}) -> Int

Return the number of simple lowering operators and z-weight coordinates.

Both counts equal the SU(N) rank `N - 1`. These values determine qlabel tuple
lengths, tableau weight lengths, and the number of sparse lowering-operator
dictionaries stored in an `Irep`.
"""
# number of raising operators (for SU(N), it is N - 1)
nlops(::Type{SU{N}}) where N = N - 1

# number of z-operators (for SU(N), it is N - 1)
nzops(::Type{SU{N}}) where N = N - 1

getsl_triv(::Type{SU{N}}, ::Type{RT}) where {N, RT<:Number} =
Tuple(Dict{NTuple{N-1, Int}, SparseMatrixCSC{RT}}() for _=1:N-1)

"""
    getsz_triv(::Type{SU{N}}) -> Dict

Return the weight-sector dictionary for the trivial SU(N) representation.

The only sector is the all-zero z-weight, and its basis index range is `(1, 1)`.
"""
getsz_triv(::Type{SU{N}}) where N = Dict(Tuple(0 for _=1:N-1)=>(1, 1))

"""
    getsl_def(::Type{SU{N}}, ::Type{RT}) -> NTuple{N-1, Dict}

Build sparse lowering-operator blocks for the defining SU(N) representation.

`RT` is the scalar type of the sparse matrices. The returned tuple has one
dictionary per simple lowering operator; each dictionary maps a source z-weight
sector to the 1x1 matrix block that applies that lowering operation.
"""
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

"""
    getsz_def(::Type{SU{N}}, i) -> NTuple{N-1, Int}

Return the z-weight of basis vector `i` in the defining representation.

`i` is one-based and should lie in `1:N`. The coordinate convention is the one
used by crystal/tableau generation and by `qlab2mwz`/`getqlabel`.
"""
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

"""
    getdz(::Type{SU{N}}, lop::Int) -> NTuple{N-1, Int}

Return the z-weight change produced by simple lowering operator `lop`.

`lop` is in `1:nlops(SU{N})`. The value is the difference between adjacent
defining-representation weights and is used to move between irrep sectors.
"""
# Change of z-values when lowering operator is applied.
function getdz(::Type{SU{N}}, lop::Int) where N
	sz = getsz_def(SU{N})
	zvals = sorted_zvals(sz)
	return zvals[lop] .- zvals[lop+1]
end

crystal_chars_map(::Type{SU{N}}) where N = Dict(i => i for i=1:N)

"""
    mw_column(::Type{<:SU}, l::Int) -> Vector{Int}

Return the canonical maximal-weight tableau column of height `l`.

The returned characters descend from `l` to `1` and are used when translating
qlabels into tableau shapes.
"""
# This function should be modified
mw_column(::Type{<:SU}, l::Int) = collect(l:-1:1)

getdzs(::Type{SU{N}}) where N = [collect(getdz(SU{N}, i)) for i=1:nlops(SU{N})]

"""
    charlist(::Type{SU{N}}) -> Vector{Int}

Return the tableau character alphabet for SU(N), ordered as `1:N`.
"""
charlist(::Type{SU{N}}) where N = collect(1:N)

"""
    get_fops_std(::Type{SU{N}}) -> Vector{Dict{Int, Int}}

Return the standard crystal lowering maps on SU(N) tableau characters.

The vector has one dictionary per simple lowering operator. Operator `i` maps
character `i` to `i + 1`; absent dictionary entries mean the operator cannot
act on that character.
"""
function get_fops_std(::Type{SU{N}}) where N
	fops = Vector{Dict{Int, Int}}()
	for i in 1:N-1
		fop = Dict{Int, Int}(); fop[i] = i + 1;
		push!(fops, fop)
	end
	return fops
end

"""
    qlab2mwz(::Type{SU{N}}, qlabel::NTuple{N-1, Int}) -> NTuple{N-1, Int}

Convert an SU(N) highest-weight qlabel to maximal-weight z-coordinates.

`qlabel[i]` is the coefficient of the `i`th fundamental weight in Dynkin-label
coordinates. The returned z-weight is used as the highest sector key in `Irep`
storage.
"""
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

"""
    getqlabel(::Type{SU{N}}, z::NTuple{N-1, Int}) -> NTuple{N-1, Int}

Convert maximal-weight z-coordinates back to an SU(N) qlabel.

`z` must lie on the SU(N) qlabel lattice; assertions check tuple length and the
required divisibility relations.
"""
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

"""
    less_weight(::Type{SU{N}}, w1::NTuple{N-1, Int}, w2::NTuple{N-1, Int}) -> Bool

Order two SU(N) weight tuples for deterministic sector traversal.

`w1` and `w2` are compared lexicographically after reversing coordinate order,
matching the tableau-weight convention used elsewhere in the package.
"""
# Weight comparision, inputs are a form of shape of the Young tableau
function less_weight(::Type{SU{N}}, w1::NTuple{NZ, Int}, w2::NTuple{NZ, Int}) where {N, NZ}
	@assert NZ == N - 1
	return reverse(w1) < reverse(w2)
end

"""
    fundamental_qlabels(::Type{SU{N}}) -> Vector{NTuple{N-1, Int}}

Return the qlabels of the fundamental SU(N) representations.

Each qlabel is a unit vector in Dynkin-label coordinates. These are used as
building blocks by tensor-product and irrep-generation routines.
"""
function fundamental_qlabels(::Type{SU{N}}) where N
	funda_lst = Vector{NTuple{N-1, Int}}()
	for i in 1:N-1
		q = Tuple(i == j ? 1 : 0 for j=1:N-1)
		push!(funda_lst, q)
	end
	return funda_lst
end
