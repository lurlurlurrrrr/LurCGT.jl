using LurCGT

include("fixedint_driver_utils.jl")

function run_fixedint_catalog_cli(numtype_name::AbstractString,
    symmetry_name::AbstractString,
    maxdim::Int;
    verbose=1)

    RT = parse_fixedint_type(numtype_name)
    symmetry = parse_nonabelian_symmetry(symmetry_name)
    return LurCGT.update_fixedint_irrep_catalog(
        symmetry,
        RT;
        maxdim=maxdim,
        base_dir=fixedint_data_root_from_env(),
        save=true,
        verbose=verbose,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 3
        println("Usage: julia --project=. test/fixedint_irep_catalog_driver.jl <NUMTYPE> <SYMMETRY> <MAXDIM>")
        println("Example: julia --project=. test/fixedint_irep_catalog_driver.jl Int64 SU3 20")
        exit(1)
    end

    numtype_name = ARGS[1]
    symmetry_name = ARGS[2]
    maxdim = parse(Int, ARGS[3])
    catalog = run_fixedint_catalog_cli(numtype_name, symmetry_name, maxdim; verbose=1)
    println("Saved $(length(catalog.accepted)) accepted irreps to $(LurCGT.fixedint_catalog_path(parse_nonabelian_symmetry(symmetry_name), parse_fixedint_type(numtype_name); base_dir=fixedint_data_root_from_env()))")
end
