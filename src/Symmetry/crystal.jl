# Tableau object consist of length(shape) rows, with shape[i] boxes in the i-th row.
# colread is a vector representing the column reading of the tableau.
# Each column is read from bottom to top, and columns are read from right to left.
struct Tableau{S<:NonabelianSymm}
    shape::Vector{Int}
    colread::Vector{Int}
end

struct Crystal_ops{S<:NonabelianSymm}
    f::Vector{Dict{Int, Int}}  # lowering operators
    e::Vector{Dict{Int, Int}}  # raising operators
    ϕ::Dict{Tuple{Int, Int}, Int} # maximal apply number for lowering operators
    ϵ::Dict{Tuple{Int, Int}, Int} # maximal apply number for raising operators
end

is_valid_tableau(::Tableau{S}) where {S<:NonabelianSymm} = true

function mw_tableau_qlabel(::Type{S}, qlabel::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ}
    @assert NZ == nlops(S)
    shape = Tableau_shape_from_qlabel(S, collect(qlabel))
    aux = get_auxarg(S, qlabel)
    return mw_tableau_shape(S, shape, aux...)
end

get_auxarg(::Type{S}, qlabel::NTuple{NZ, Int}) where {S<:NonabelianSymm, NZ} = ()

function mw_tableau_shape(::Type{S}, shape::Vector{Int}, aux...) where S <: NonabelianSymm
    colread = Int[]
    before_val = 0
    for i in length(shape):-1:1
        ncol = shape[i] - before_val
        if ncol == 0 continue end
        col = mw_column(S, i, aux...)
        for j in 1:ncol append!(colread, col) end
        before_val = shape[i]
    end
    return Tableau{S}(shape, colread)
end

# This function can be modified for different symmetries
function Tableau_shape_from_qlabel(::Type{S}, qlabel::Vector{Int}) where S <: NonabelianSymm
    @assert length(qlabel) == nlops(S)
    preprocessed = preprocess(S, qlabel)
    shape = cumsum(reverse(preprocessed))
    return append!(reverse(shape))
end

preprocess(::Type{S}, qlabel::Vector{Int}) where S <: NonabelianSymm = qlabel

# Pretty print for Tableau
function Base.show(io::IO, tab::Tableau)
    println(io, "Shape: ", tab.shape)
    # Currently, it is represented in column reading form. 
    # Print in row form
    println(io, "Table:")
    nrows = length(tab.shape)
    rows = Vector{Vector{Int}}(undef, nrows)

    for i in 1:nrows rows[i] = Vector{Int}() end
    a = mw_tableau_shape(SU{1000}, tab.shape).colread

    for ii in 1:length(tab.colread)
        col = a[ii]; push!(rows[col], tab.colread[ii])
    end

    for i in 1:nrows
        row = rows[i]
        for elem in row
            print(io, lpad(elem, 4))
        end
        println(io)
    end
end

function get_crystal_ops(::Type{S}) where S <: NonabelianSymm
    char_list = charlist(S)
    f = get_fops_std(S)
    e = get_eops_std(f)
    ϕ = maximal_apply_num(char_list, f)
    ϵ = maximal_apply_num(char_list, e)
    
    return Crystal_ops{S}(f, e, ϕ, ϵ)
end

function maximal_apply_num(char_list::Vector{Int}, ops::Vector{Dict{Int, Int}})
    # Nonzero values are stored in a dictionary
    # max_apply[(i, char)] = n means ϕ_i(char) = n
    max_apply = Dict{Tuple{Int, Int}, Int}()
    for (i, op) in enumerate(ops)
        apply_num = Dict{Int, Int}()
        for char in char_list
            count = 0
            current = char
            while haskey(op, current)
                current = op[current]
                count += 1
            end
            if count > 0 max_apply[(i, char)] = count end
        end
    end
    return max_apply
end

# Generate the raising operators from the lowering operators
function get_eops_std(fops_std::Vector{Dict{Int, Int}})
    eops_std = Vector{Dict{Int, Int}}()
    for fop in fops_std
        eop = Dict{Int, Int}()
        for (k, v) in fop
            eop[v] = k
        end
        push!(eops_std, eop)
    end
    return eops_std
end

# In-place application of 'opi'th lowering operator to tableau
function apply_lowering!(tab::Tableau{S}, 
    cops::Crystal_ops{S},
    opi::Int) where S <: NonabelianSymm
    # Apply lowering operator lop to tableau tab

    f, e = cops.f, cops.e
    ϕ, ϵ = cops.ϕ, cops.ϵ
    colread = tab.colread
    max_ind, max_val = 0, -Inf
    current_val = 0
    for i in 1:length(colread)
        current_val += get(ϕ, (opi, colread[i]), 0)
        if i > 1 current_val -= get(ϵ, (opi, colread[i-1]), 0) end
        if current_val > max_val
            max_val = current_val
            max_ind = i
        end
    end
    if !haskey(f[opi], colread[max_ind]) return false end
    oldval = colread[max_ind]
    colread[max_ind] = f[opi][oldval]
    if !is_valid_tableau(tab)
        colread[max_ind] = oldval
        return false
    end
    return true
end

# Non in-place application of 'opi'th lowering operator to tableau
function apply_lowering(tab::Tableau{S}, 
    cops::Crystal_ops{S},
    opi::Int) where S <: NonabelianSymm
    new_tab = Tableau{S}(copy(tab.shape), copy(tab.colread))
    success = apply_lowering!(new_tab, cops, opi)
    if success return new_tab else return nothing end
end

function get_weight(tab::Tableau{S}) where S <: NonabelianSymm
    cnt::Dict{Int, Int} = get_count(tab)
    NZ = nzops(S)
    szs = getsz_def_vec(S)
    chardict = crystal_chars_map(S)
    @assert length(szs[1]) == NZ
    weight = zeros(Int, NZ)
    for char in keys(cnt)
        weight += cnt[char] * szs[chardict[char]]
    end
    return Tuple(weight)
end

function det_getvector(tab::Tableau{S}, 
    i::Int,
    cops::Crystal_ops{S},
    tabdict::Dict{Tableau{S}, Vector{RT}}) where {S<:NonabelianSymm, RT}

    while true
        tab = apply_lowering(tab, cops, i)
        if tab === nothing break end
        if !haskey(tabdict, tab) return true end
    end
    return false
end

function get_count(tab::Tableau{S}) where S <: NonabelianSymm
    cnt::Dict{Int, Int} = Dict{Int, Int}()
    for val in tab.colread
        cnt[val] = get(cnt, val, 0) + 1
    end
    return cnt
end

# Needed to define hash for Tableau to be used as Dict key
function Base.hash(tab::Tableau{S}, h::UInt) where S<:NonabelianSymm
    h′ = h
    for e in tab.colread h′ = hash(e, h′) end
    return h′
end

function Base.:(==)(t1::Tableau{S}, t2::Tableau{S}) where S<:NonabelianSymm
    return t1.shape == t2.shape && t1.colread == t2.colread
end
