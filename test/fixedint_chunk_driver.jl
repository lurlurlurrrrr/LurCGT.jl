using LurCGT

include("fixedint_driver_utils.jl")

function run_fixedint_chunk_cli(numtype_name::AbstractString,
    symmetry_name::AbstractString,
    d1min::Int,
    d1max::Int,
    d2min::Int,
    d2max::Int,
    m1::Int,
    m2::Int,
    chunk1::Int,
    chunk2::Int;
    update_catalog=false,
    merge_local_ireps=false,
    verbose=1)

    RT = parse_fixedint_type(numtype_name)
    symmetry = parse_nonabelian_symmetry(symmetry_name)
    return LurCGT.run_fixedint_cgt_chunk(
        symmetry,
        RT,
        d1min,
        d1max,
        d2min,
        d2max,
        m1,
        m2,
        chunk1,
        chunk2;
        base_dir=fixedint_data_root_from_env(),
        save=true,
        update_catalog=update_catalog,
        merge_local_ireps=merge_local_ireps,
        verbose=verbose,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 10
        println("Usage: julia --project=. test/fixedint_chunk_driver.jl <NUMTYPE> <SYMMETRY> <D1MIN> <D1MAX> <D2MIN> <D2MAX> <M1> <M2> <CHUNK1> <CHUNK2>")
        println("Example: julia --project=. test/fixedint_chunk_driver.jl Int128 SU3 1 20 1 20 4 4 1 2")
        println("Run test/fixedint_irep_catalog_driver.jl first to create the catalog.")
        exit(1)
    end

    numtype_name = ARGS[1]
    symmetry_name = ARGS[2]
    values = parse.(Int, ARGS[3:end])
    summary = run_fixedint_chunk_cli(numtype_name, symmetry_name, values...; verbose=1)
    println("Chunk $(summary.dim_range1) vs $(summary.dim_range2): $(summary.passed_pairs) passed, $(summary.failed_pairs) failed")
    exit(summary.failed_pairs == 0 ? 0 : 1)
end
