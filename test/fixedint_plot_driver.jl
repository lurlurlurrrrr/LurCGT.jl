using LurCGT
using Plots

include("fixedint_driver_utils.jl")

function run_fixedint_plot_cli(numtype_name::AbstractString,
    symmetry_name::AbstractString;
    output_path=nothing)

    RT = parse_fixedint_type(numtype_name)
    symmetry = parse_nonabelian_symmetry(symmetry_name)
    base_dir = fixedint_data_root_from_env()
    cells = LurCGT.collect_fixedint_plot_cells(symmetry, RT; base_dir=base_dir)

    plot_dir = joinpath(base_dir, numtype_name, LurCGT.totxt(symmetry), "plots")
    mkpath(plot_dir)
    final_path = isnothing(output_path) ? joinpath(plot_dir, "dim_success.png") : output_path

    xs = Int[cell.dim1 for cell in cells]
    ys = Int[cell.dim2 for cell in cells]
    colors = [cell.status == :passed ? :green : :red for cell in cells]
    plt = scatter(
        xs,
        ys;
        markercolor=colors,
        markersize=6,
        legend=false,
        xlabel="dim(q1)",
        ylabel="dim(q2)",
        title="$(numtype_name) $(LurCGT.totxt(symmetry)) CG chunk results",
    )
    savefig(plt, final_path)
    return final_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    if !(length(ARGS) == 2 || length(ARGS) == 3)
        println("Usage: julia --project=. test/fixedint_plot_driver.jl <NUMTYPE> <SYMMETRY> [OUTPUT_PATH]")
        println("Example: julia --project=. test/fixedint_plot_driver.jl Int128 SU3")
        exit(1)
    end

    output_path = length(ARGS) == 3 ? ARGS[3] : nothing
    final_path = run_fixedint_plot_cli(ARGS[1], ARGS[2]; output_path=output_path)
    println("Saved plot to $final_path")
end
