using Test
using LurCGT

include("int128_frontier_utils.jl")

@testset "Int128 frontier qlabel steps" begin
    @test sort(int128_frontier_cover_steps(SU{3})) == [(0, 1), (1, 0)]
    @test sort(int128_frontier_cover_steps(Sp{4})) == [(0, 1), (1, 0)]
    @test sort(int128_frontier_cover_steps(SO{5})) == [(0, 2), (1, 0)]
    @test sort(int128_frontier_cover_steps(SO{6})) == [(0, 0, 2), (0, 1, 1), (0, 2, 0), (1, 0, 0)]
end

@testset "Int128 frontier search" begin
    cap = ((0, 1), (2, 1))
    frontier = find_int128_cgt_frontier(
        SU{3};
        pass_predicate = (q1, q2) -> all(q1 .<= cap[1]) && all(q2 .<= cap[2]),
        max_tested = 200,
        verbose = 0,
    )
    @test frontier == [cap]
end

@testset "Int128 frontier irep comparison" begin
    rep = getNsave_irep(SU{2}, BigInt, (1,))
    rep_same = LurCGT.Irep{SU{2}, 1, 1, BigInt}(
        rep.Sl,
        rep.Sz,
        rep.innerprod,
        rep.inv_innerprod,
        rep.qlabel,
        rep.dimension,
        rep.size_byte + 1,
    )
    rep_diff = LurCGT.Irep{SU{2}, 1, 1, BigInt}(
        rep.Sl,
        rep.Sz,
        rep.innerprod,
        rep.inv_innerprod,
        rep.qlabel,
        rep.dimension + 1,
        rep.size_byte,
    )
    @test int128_frontier_irep_equal_except_size_byte(rep, rep_same)
    @test !int128_frontier_irep_equal_except_size_byte(rep, rep_diff)
    @test int128_frontier_irep_matches(SU{2}, (1,); verbose=0)
end

@testset "Int128 frontier real generator smoke" begin
    @test int128_frontier_passes(SU{2}, (1,), (1,); verbose=0)
end
