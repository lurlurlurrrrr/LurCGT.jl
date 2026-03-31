using LurCGT

include("fixedint_driver_utils.jl")

function parse_fixedint_catalog_mode(name::AbstractString)
    name == "summary" && return :summary
    name == "accepted" && return :accepted
    name == "scanned" && return :scanned
    name == "all" && return :all
    throw(ArgumentError("mode must be summary, accepted, scanned, or all"))
end

function print_fixedint_catalog(io::IO, catalog; mode::Symbol=:accepted)
    println(io, "Catalog: $(catalog.numtype) $(catalog.symmetry)")
    println(io, "Accepted: $(length(catalog.accepted))")
    println(io, "Scanned: $(length(catalog.scanned))")

    if mode == :accepted || mode == :all
        for entry in catalog.accepted
            println(io, "accepted dim=$(entry.dim) q=$(entry.qlabel)")
        end
    end

    if mode == :scanned || mode == :all
        qlabels = sort!(collect(keys(catalog.scanned)); by=q -> (catalog.scanned[q].dim, q))
        for q in qlabels
            status = catalog.scanned[q]
            println(io, "$(status.status) dim=$(status.dim) q=$(q) reason=$(status.reason)")
        end
    end

    return nothing
end

function run_fixedint_catalog_read_cli(numtype_name::AbstractString,
    symmetry_name::AbstractString;
    mode="accepted",
    io::IO=stdout)

    RT = parse_fixedint_type(numtype_name)
    symmetry = parse_nonabelian_symmetry(symmetry_name)
    base_dir = fixedint_data_root_from_env()
    path = LurCGT.fixedint_catalog_path(symmetry, RT; base_dir=base_dir)
    isfile(path) || throw(ArgumentError("No catalog found at $path. Run test/fixedint_irep_catalog_driver.jl first."))

    catalog = LurCGT.load_fixedint_catalog(symmetry, RT; base_dir=base_dir)
    parsed_mode = parse_fixedint_catalog_mode(mode)
    println(io, "Path: $path")
    print_fixedint_catalog(io, catalog; mode=parsed_mode)
    return catalog
end

if abspath(PROGRAM_FILE) == @__FILE__
    if !(length(ARGS) == 2 || length(ARGS) == 3)
        println("Usage: julia --project=. test/fixedint_catalog_read_driver.jl <NUMTYPE> <SYMMETRY> [MODE]")
        println("Example: julia --project=. test/fixedint_catalog_read_driver.jl Int128 SU3 accepted")
        println("Modes: summary, accepted, scanned, all")
        exit(1)
    end

    mode = length(ARGS) == 3 ? ARGS[3] : "accepted"
    run_fixedint_catalog_read_cli(ARGS[1], ARGS[2]; mode=mode)
end
