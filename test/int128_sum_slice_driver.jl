using LurCGT

function parse_nonabelian_symmetry(symmetry_name::AbstractString)
    m = match(r"^(SU|SO|Sp)(\d+)$", symmetry_name)
    isnothing(m) && throw(ArgumentError("symmetry must look like SU3, SO5, or Sp4"))
    family, rank_text = m.captures
    rank = parse(Int, rank_text)
    if family == "SU"
        return Core.apply_type(LurCGT.SU, rank)
    elseif family == "SO"
        return Core.apply_type(LurCGT.SO, rank)
    end
    return Core.apply_type(LurCGT.Sp, rank)
end

function run_int128_cgt_sum_slice_cli(symmetry_name::AbstractString,
    n1::Int,
    n2::Int;
    verbose=1)

    symmetry = parse_nonabelian_symmetry(symmetry_name)
    return LurCGT.run_int128_cgt_sum_slice(symmetry, n1, n2; verbose)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 3
        println("Usage: julia --project=. test/int128_sum_slice_driver.jl <SYMMETRY> <N1> <N2>")
        println("Example: julia --project=. test/int128_sum_slice_driver.jl SU3 5 2")
        exit(1)
    end

    symmetry_name = ARGS[1]
    n1 = parse(Int, ARGS[2])
    n2 = parse(Int, ARGS[3])
    summary = run_int128_cgt_sum_slice_cli(symmetry_name, n1, n2; verbose=1)
    exit(summary.failed_pairs == 0 ? 0 : 1)
end
