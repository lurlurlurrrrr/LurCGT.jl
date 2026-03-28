using LurCGT

include("int128_frontier_utils.jl")

function run_int128_cgt_frontier(::Type{S};
    max_tested::Int=typemax(Int),
    verbose=1) where {S<:NonabelianSymm}

    frontier = find_int128_cgt_frontier(S; max_tested, verbose)
    println("Int128 CGT frontier for $(S):")
    for (q1, q2) in frontier
        println("  ($(q1), $(q2))")
    end
    return frontier
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("Example:")
    println("  include(\"test/int128_frontier_driver.jl\")")
    println("  run_int128_cgt_frontier(SU{3}; max_tested=5000, verbose=1)")
end
