using LurCGT
using SparseArrayKit
using TensorOperations
using LinearAlgebra
using SparseArrays
using Random
using StatsBase
using Test
using TOML

include("test_utils.jl")

function imported_packages(src_root)
    imported = Set{String}()
    for path in readdir(src_root; join=true)
        if isdir(path)
            union!(imported, imported_packages(path))
            continue
        end
        endswith(path, ".jl") || continue
        for line in eachline(path)
            match(r"^\s*(using|import)\s+([A-Za-z0-9_]+)", line) === nothing && continue
            pkg = match(r"^\s*(using|import)\s+([A-Za-z0-9_]+)", line).captures[2]
            startswith(pkg, ".") && continue
            push!(imported, pkg)
        end
    end
    return imported
end

@testset "package smoke" begin
    @test isdefined(LurCGT, :SU)
    @test isdefined(LurCGT, :getNsave_irep)
    @test isabelian(U1)
end

@testset "CGTSVD tests" begin
    test_CGTSVD_sameq_randinput(SU{2}, 5000, 5, 100; verbose=0)
    test_CGTSVD_sameq_randinput(SU{3}, 100000, 2, 50; verbose=0)
    test_CGTSVD_basischange(SU{2}, 5000, 8, 30; verbose=0)
    test_CGTSVD_basischange(SU{3}, 20000, 3, 10; verbose=0)

    test_CGTSVD_unitary_randinput(SU{2}, 5000, 8, 30; verbose=0)
    test_CGTSVD_unitary_randinput(SU{3}, 20000, 3, 10; verbose=0)
end


@testset "declared deps cover imports" begin
    project = TOML.parsefile(joinpath(dirname(@__DIR__), "Project.toml"))
    declared = Set(keys(get(project, "deps", Dict())))
    stdlibs = Set(["Base64", "Dates", "DelimitedFiles", "Downloads", "FileWatching",
                   "Libdl", "LinearAlgebra", "Logging", "Markdown", "Mmap",
                   "Pkg", "Printf", "Random", "REPL", "Serialization",
                   "SHA", "Sockets", "SparseArrays", "Statistics", "TOML",
                   "Tar", "Test", "UUIDs"])

    imports = imported_packages(joinpath(dirname(@__DIR__), "src"))
    undeclared = sort!(collect(setdiff(imports, union(declared, stdlibs))))
    @test isempty(undeclared)
end

function test_Fsym_unitarity(::Type{S}, test_inputs=100, qlimit=4) where S<:NonabelianSymm
    for i in 1:test_inputs
        println("Test #$i started:")
        ins, out, fsym_mat = get_random_Fsymbol_real(S, Float64, qlimit)
        sz = size(fsym_mat, 1)
        @assert fsym_mat * fsym_mat' ≈ Matrix(I, sz, sz) atol=1e-10 rtol=1e-10 "F-symbol is not unitary for input $ins and output $out"
        println("Test #$i passed. Input labels: $ins, output label: $out")
    end
end

function test_Rsym_unitarity(::Type{S}, test_inputs=30, qlimit=4) where S<:NonabelianSymm
    for i in 1:test_inputs
        in, out, rsym_mat = get_random_Rsymbol_real(S, Float64, qlimit)
        sz = size(rsym_mat, 1)
        @assert rsym_mat * rsym_mat' ≈ Matrix(I, sz, sz) atol=1e-10 rtol=1e-10 "R-symbol is not unitary for input $ins and output $out"
        println("Test #$i passed. Input label: $in, output label: $out")
    end
end

@testset "CGTperm test" begin
    test_CGTperm(SU{2}, 10000, 8, 100; verbose=0)
    test_CGTperm(SU{3}, 500000, 3, 20; verbose=0)
end

@testset "Conjugation of CGT test (same input and output)" begin
    test_CGT_conj_sameq(SU{3}, 1000, 3, 100, 20; verbose=0)
    test_CGT_conj_sameq(SU{2}, 100, 8, 20, 100; verbose=0)
end

@testset "Conjugation of CGT test" begin
    test_CGT_conj(SU{3}, 500000, 3, 100, 20; verbose=0)
    test_CGT_conj(SU{2}, 5000, 8, 20, 100; verbose=0)
end


# Test pentagon equation for F-symbols
@testset "F-symbol pentagon equation tests" begin
    test_pentagon_randinput(SU{2}, BigFloat, 6, 3000; verbose=0)
    test_pentagon_randinput(SU{3}, BigFloat, 2, 100; verbose=0)
end

@testset "X-symbol with 1j tests" begin
    test_Xsym_1j(SU{2}, 1000, 10, 20, 1000; verbose=0)
    test_Xsym_1j(SU{3}, 300000, 3, 100, 20; verbose=0)
end

@testset "R-symbol from 1j tests" begin
    test_rsym_from_1j(SU{2}, 15, 8)
    test_rsym_from_1j(SU{3}, 5, 5)
    test_rsym_from_1j(Sp{4}, 3, 6)
end

sp6_dim = Dict(
    (0,0,0) => 1, (1,0,0) => 6, (0,0,1) => 14, (0,1,0) => 14,
    (2,0,0) => 21, (3,0,0) => 56, (1,1,0) => 64, (1,0,1) => 70,
    (0,0,2) => 84, (0,2,0) => 90, (0,1,1) => 126, (4,0,0) => 126,
    (2,1,0) => 189, (2,0,1) => 216, (0,0,3) => 330, (1,2,0) => 350,
    (1,0,2) => 378, (0,3,0) => 385, (3,1,0) => 448, (1,1,1) => 512,
    (3,0,1) => 525, (0,1,2) => 594, (0,2,1) => 616, (2,2,0) => 924,
    (2,0,2) => 1078, (1,3,0) => 1344, (2,1,1) => 1386, (1,2,1) => 2205,
    (1,1,2) => 2240
)

# Test whether the computed dimensions of Sp(6) irreps match known values
@testset "Sp(6) irrep dimensionality tests" begin
    for (qlabel, dim) in sp6_dim
        rep = getNsave_irep(Sp{6}, BigInt, qlabel)
        computed_dim = LurCGT.dimension(rep)
        @assert computed_dim == dim 
        println("Sp(6) irrep $qlabel: dimension $computed_dim verified.")
    end
end

@testset "F-symbol unitarity tests" begin
    test_Fsym_unitarity(SU{3}, 20, 4)
    test_Fsym_unitarity(SU{2}, 100, 10)
end

@testset "R-symbol unitarity tests" begin
    test_Rsym_unitarity(SU{2}, 30, 10)
    test_Rsym_unitarity(SU{3}, 20, 4)
end

@testset "permutation of FTree test" begin
    test_permutation_randinput(SU{2}, 10000, 20, 100; verbose=0)
    test_permutation_randinput(SU{3}, 100000, 100, 20; verbose=0)
end

@testset "X-symbol (normalized) test" begin
    test_Xsym(SU{3}, 10000, 3, 20; verbose=0)
    test_Xsym(SU{2}, 1000, 8, 100; verbose=0)
    LurCGT.merge_all_to_global(SU{2})
    LurCGT.merge_all_to_global(SU{3})
end


