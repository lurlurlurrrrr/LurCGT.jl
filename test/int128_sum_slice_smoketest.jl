using Test
using LurCGT

@testset "Int128 sum-slice qlabels" begin
    @test LurCGT.int128_sum_slice_qlabels(SU{3}, 2) == [(0, 2), (1, 1), (2, 0)]
    @test LurCGT.int128_sum_slice_qlabels(SO{5}, 2) == [(0, 2), (2, 0)]
    @test LurCGT.int128_sum_slice_qlabels(SU{2}, 0) == [(0,)]
end

@testset "Int128 sum-slice pairs" begin
    @test LurCGT.int128_sum_slice_pairs(SU{2}, 2, 1) == [((1,), (2,))]
    @test LurCGT.int128_sum_slice_pairs(SU{3}, 1, 1) == [
        ((0, 1), (0, 1)),
        ((0, 1), (1, 0)),
        ((1, 0), (1, 0)),
    ]
end

@testset "Int128 sum-slice real generator smoke" begin
    summary = LurCGT.run_int128_cgt_sum_slice(SU{2}, 1, 1; verbose=0)
    @test summary.n1 == 1
    @test summary.n2 == 1
    @test summary.total_pairs == 1
    @test summary.passed_pairs == 1
    @test summary.failed_pairs == 0
    @test isempty(summary.failures)
end

@testset "Int128 sum-slice quiet output" begin
    pipe = Pipe()
    redirect_stdout(pipe) do
        LurCGT.run_int128_cgt_sum_slice(SU{2}, 1, 1; verbose=0)
    end
    close(pipe.in)
    output = read(pipe, String)
    @test !occursin("Dict{Tuple", output)
end
