"""
    isvalidsymm(::Type{Sp{N}}) -> Bool

Return whether `Sp{N}` is supported by this implementation.

`N` is the defining dimension and must be positive even. The rank, qlabel
length, and number of simple lowering operators are all `div(N, 2)`.
"""
isvalidsymm(::Type{Sp{N}}) where N = 2 <= N && N % 2 == 0

isSpN(::Type{<:Sp}) = true
isSpN(::Any) = false

# dimension of defining irep
defirepdim(::Type{Sp{N}}) where N = N

# Text expression of the symmetry.
totxt(::Type{Sp{N}}) where N = "Sp$(N)"

# number of raising operators (for Sp(N), it is N / 2)
nlops(::Type{Sp{N}}) where N = div(N, 2)

# number of z-operators (for Sp(N), it is N / 2)
nzops(::Type{Sp{N}}) where N = div(N, 2)

maxalt(::Type{Sp{N}}) where N = div(N, 2)

getsl_triv(::Type{Sp{N}}, ::Type{RT}) where {N, RT<:Number} =
Tuple(Dict{NTuple{div(N, 2), Int}, SparseMatrixCSC{RT}}() for _=1:div(N, 2))

"""
    getsz_triv(::Type{Sp{N}}) -> Dict

Return the trivial Sp(N) weight-sector dictionary.

The only sector is the all-zero z-weight of length `div(N, 2)` with one basis
state at range `(1, 1)`.
"""
getsz_triv(::Type{Sp{N}}) where N = Dict(Tuple(0 for _=1:div(N, 2))=>(1, 1))

"""
    getsl_def(::Type{Sp{N}}, ::Type{RT}) -> NTuple{div(N,2), Dict}

Build sparse lowering-operator blocks for the defining Sp(N) representation.

`RT` is the scalar type of the sparse matrices. The returned tuple has one
dictionary per simple lowering operator and records the source z-weight sectors
where the defining-representation lowering block is nonzero.
"""
function getsl_def(::Type{Sp{N}}, ::Type{RT}) where {N, RT<:Number} 
    def_sl = Tuple(Dict{NTuple{div(N, 2), Int}, SparseMatrixCSC{RT}}() for _=1:div(N, 2))
    def_sz = getsz_def(Sp{N})
    def_zvals = sorted_zvals(def_sz)
    for i=1:div(N, 2)
        sz = def_zvals[i]
        val = spzeros(RT, 1, 1)
        val[1, 1] = RT(1) 
        def_sl[i][sz] = val
        if i < div(N, 2)
            sz = def_zvals[N-i]
            val = spzeros(RT, 1, 1)
            val[1, 1] = RT(1) 
            def_sl[i][sz] = val
        end
    end
    return def_sl
end

"""
    getsz_def(::Type{Sp{N}}, i) -> NTuple{div(N,2), Int}

Return the z-weight of basis vector `i` in the defining Sp(N) representation.

`i` is one-based in `1:N`. The first half of the defining basis carries positive
weights and the second half carries the corresponding negated weights.
"""
function getsz_def(::Type{Sp{N}}, i) where N 
    @assert 1 <= i <= N
    if i > div(N, 2) ii = N - i + 1 else ii = i end

    sz = zeros(Int, div(N, 2))
    for j in 1:div(N, 2)
        if j >= ii
            sz[j] = 1
        elseif j == ii - 1
            sz[j] = -j
        else
            sz[j] = 0
        end
    end
    if i > div(N, 2) sz *= -1 end
    return Tuple(sz)
end
getsz_def(::Type{Sp{N}}) where N = Dict(Tuple(getsz_def(Sp{N}, i)) => (i, i) for i=1:N)
getsz_def_vec(::Type{Sp{N}}) where N = [collect(getsz_def(Sp{N}, i)) for i=1:N]

"""
    getdz(::Type{Sp{N}}, lop::Int) -> NTuple{div(N,2), Int}

Return the z-weight change produced by Sp(N) lowering operator `lop`.

`lop` must be in `1:nlops(Sp{N})`. The result is computed from adjacent sorted
defining-representation weights and drives irrep-sector traversal.
"""
function getdz(::Type{Sp{N}}, lop::Int) where N
    @assert 1 <= lop <= div(N, 2)
    sz = getsz_def(Sp{N})
    zvals = sorted_zvals(sz)
    return zvals[lop] .- zvals[lop+1]
end

"""
    crystal_chars_map(::Type{Sp{N}}) -> Dict{Int, Int}

Map Sp(N) tableau characters to defining-representation basis positions.

Positive characters identify the first half of the defining basis and negative
characters identify the mirrored second half.
"""
function crystal_chars_map(::Type{Sp{N}}) where N
    charmap = Dict{Int, Int}()
    for i in 1:div(N, 2) charmap[i] = i end
    for i in div(N, 2)+1:N charmap[i-N-1] = i end
    return charmap
end

"""
    mw_column(::Type{<:Sp}, l::Int) -> Vector{Int}

Return the canonical maximal-weight tableau column of height `l`.

The column characters descend from `l` to `1`, matching the symplectic tableau
convention used to seed crystal generation.
"""
# This function should be modified
mw_column(::Type{<:Sp}, l::Int) = collect(l:-1:1)

getdzs(::Type{Sp{N}}) where N = [collect(getdz(Sp{N}, i)) for i=1:nlops(Sp{N})]

"""
    charlist(::Type{Sp{N}}) -> Vector{Int}

Return the ordered Sp(N) tableau character alphabet.

The list contains positive characters followed by negative characters in
descending absolute value.
"""
charlist(::Type{Sp{N}}) where N = vcat(collect(1:div(N, 2)), collect(-div(N, 2):-1))

"""
    get_fops_std(::Type{Sp{N}}) -> Vector{Dict{Int, Int}}

Return standard crystal lowering maps for Sp(N) tableau characters.

For non-terminal simple roots, the map lowers both positive and mirrored
negative characters. The final simple root maps the terminal positive character
to its negative partner.
"""
# f-operations defined for crystal of tableau
function get_fops_std(::Type{Sp{N}}) where N
    fops = Vector{Dict{Int, Int}}()
    for i in 1:div(N, 2)
        fop = Dict{Int, Int}()
        if i < div(N, 2) fop[i] = i + 1; fop[-i-1] = -i
        else fop[i] = -i end
        push!(fops, fop)
    end
    return fops
end 

"""
    qlab2mwz(::Type{Sp{N}}, qlabel::NTuple{div(N,2), Int}) -> NTuple{div(N,2), Int}

Convert an Sp(N) highest-weight qlabel to maximal-weight z-coordinates.

`qlabel[i]` is the coefficient of the `i`th fundamental weight. The returned
z-weight is the key used for highest sectors in generated irreps.
"""
function qlab2mwz(::Type{Sp{N}}, qlabel::NTuple{NZ, Int}) where {N, NZ}
    @assert div(N, 2) == NZ
    sz_defs = [collect(getsz_def(Sp{N}, i)) for i in 1:NZ]
    z = zeros(Int, NZ)
    z_added = zeros(Int, NZ)
    for i in 1:NZ
        z_added .+= sz_defs[i]
        z .+= z_added .* qlabel[i]
    end
    return Tuple(z)
end

"""
    getqlabel(::Type{Sp{N}}, z::NTuple{div(N,2), Int}) -> NTuple{div(N,2), Int}

Convert maximal-weight z-coordinates back to an Sp(N) qlabel.

`z` must satisfy the divisibility relations of the package's symplectic weight
lattice; assertions enforce those conditions.
"""
function getqlabel(::Type{Sp{N}}, z::NTuple{NZ, Int}) where {N, NZ}
    @assert div(N, 2) == NZ
    w = zeros(Int, NZ)
    for i=NZ:-1:2
        @assert (z[i] - z[i-1]) % i == 0
        w[i] = div(z[i] - z[i-1], i)
    end
    @assert z[1] % 1 == 0
    w[1] = z[1]
    return Tuple(w)
end 

"""
    less_weight(::Type{Sp{N}}, w1::NTuple{div(N,2), Int}, w2::NTuple{div(N,2), Int}) -> Bool

Order Sp(N) weight tuples for deterministic sector traversal.

The comparison reverses coordinate order before lexicographic comparison,
matching the tableau-shape convention.
"""
# Weight comparision, inputs are a form of shape of the Young tableau
function less_weight(::Type{Sp{N}}, w1::NTuple{NZ, Int}, w2::NTuple{NZ, Int}) where {N, NZ}
	@assert NZ == div(N, 2)
	return reverse(w1) < reverse(w2)
end

"""
    fundamental_qlabels(::Type{Sp{N}}) -> Vector{NTuple{div(N,2), Int}}

Return the fundamental Sp(N) qlabels.

Each qlabel is a unit vector in fundamental-weight coordinates. These qlabels
seed tensor-product decomposition and irrep generation.
"""
function fundamental_qlabels(::Type{Sp{N}}) where N
	funda_lst = Vector{NTuple{div(N, 2), Int}}()
	for i in 1:div(N, 2)
		q = Tuple(i == j ? 1 : 0 for j=1:div(N, 2))
		push!(funda_lst, q)
	end
	return funda_lst
end
