isvalidsymm(::Type{G2}) = true

totxt(::Type{G2}) = "G2"

defirepdim(::Type{G2}) = 7

nlops(::Type{G2}) = 2
nzops(::Type{G2}) = 2

maxalt(::Type{G2}) = 2

const G2_DEF_CHARS = (1, 2, 3, 0, -3, -2, -1)
const G2_DEF_WEIGHTS = (
    (1, 1),
    (-1, 1),
    (2, 0),
    (0, 0),
    (-2, 0),
    (1, -1),
    (-1, -1),
)
const G2_FOPS = (
    Dict(1 => 2, 3 => 0, 0 => -3, -2 => -1),
    Dict(2 => 3, -3 => -2),
)

getsl_triv(::Type{G2}, ::Type{RT}) where {RT<:Number} =
    Tuple(Dict{NTuple{2, Int}, SparseMatrixCSC{RT}}() for _=1:2)

getsz_triv(::Type{G2}) = Dict((0, 0) => (1, 1))

getsz_def(::Type{G2}, i::Int) = G2_DEF_WEIGHTS[i]

getsz_def(::Type{G2}) = Dict(getsz_def(G2, i) => (i, i) for i=1:defirepdim(G2))

getsz_def_vec(::Type{G2}) = [collect(getsz_def(G2, i)) for i=1:defirepdim(G2)]

function getsl_def(::Type{G2}, ::Type{RT}) where {RT<:Number}
    def_sl = Tuple(Dict{NTuple{2, Int}, SparseMatrixCSC{RT}}() for _=1:2)
    char_weight = Dict(char => weight for (char, weight) in zip(G2_DEF_CHARS, G2_DEF_WEIGHTS))
    for (lop, fop) in enumerate(G2_FOPS)
        for (src, _) in fop
            val = spzeros(RT, 1, 1)
            val[1, 1] = RT(1)
            def_sl[lop][char_weight[src]] = val
        end
    end
    def_sl[1][(0, 0)] = RT[2;;]
    return def_sl
end

function getdz(::Type{G2}, lop::Int)
    @assert 1 <= lop <= 2
    return lop == 1 ? (2, 0) : (-3, 1)
end

crystal_chars_map(::Type{G2}) = Dict(char => i for (i, char) in enumerate(G2_DEF_CHARS))

function mw_column(::Type{G2}, l::Int)
    @assert 1 <= l <= 2
    return l == 1 ? [1] : [2, 1]
end

getdzs(::Type{G2}) = [collect(getdz(G2, i)) for i=1:nlops(G2)]

charlist(::Type{G2}) = collect(G2_DEF_CHARS)

get_fops_std(::Type{G2}) = [copy(G2_FOPS[1]), copy(G2_FOPS[2])]

function qlab2mwz(::Type{G2}, qlabel::NTuple{NZ, Int}) where {NZ}
    @assert NZ == 2
    a, b = qlabel
    return (a, a + 2 * b)
end

function getqlabel(::Type{G2}, z::NTuple{NZ, Int}) where {NZ}
    @assert NZ == 2
    z1, z2 = z
    @assert iseven(z2 - z1)
    return (z1, div(z2 - z1, 2))
end

function less_weight(::Type{G2}, w1::NTuple{NZ, Int}, w2::NTuple{NZ, Int}) where {NZ}
    @assert NZ == 2
    return reverse(w1) < reverse(w2)
end

fundamental_qlabels(::Type{G2}) = [(1, 0), (0, 1)]

function getdefirep(::Type{G2}, ::Type{RT}) where {RT<:Number}
    def_sl = getsl_def(G2, RT)
    def_sz = getsz_def(G2)
    innerprod = get_identity_innerprod(G2, RT, def_sz)
    innerprod[(0, 0)] = RT[2;;]
    inv_innerprod = get_identity_invinprod(G2, RT, def_sz)
    inv_innerprod[(0, 0)] = (RT[1;;], 1//2)
    dim = defirepdim(G2)
    return Irep{G2, 2, 2, RT}(def_sl, def_sz, innerprod, inv_innerprod, (1, 0), dim)
end
