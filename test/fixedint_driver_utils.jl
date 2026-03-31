using LurCGT

function parse_fixedint_type(name::AbstractString)
    name == "Int64" && return Int64
    name == "Int128" && return Int128
    throw(ArgumentError("numtype must be Int64 or Int128"))
end

function parse_nonabelian_symmetry(symmetry_name::AbstractString)
    symmetry_name == "G2" && return LurCGT.G2
    m = match(r"^(SU|SO|Sp)(\d+)$", symmetry_name)
    isnothing(m) && throw(ArgumentError("symmetry must look like SU3, SO5, Sp4, or G2"))
    family, rank_text = m.captures
    rank = parse(Int, rank_text)
    if family == "SU"
        return Core.apply_type(LurCGT.SU, rank)
    elseif family == "SO"
        return Core.apply_type(LurCGT.SO, rank)
    end
    return Core.apply_type(LurCGT.Sp, rank)
end

fixedint_data_root_from_env() = get(ENV, "FIXEDINT_DATA_ROOT", joinpath(@__DIR__, "fixedint_data"))
