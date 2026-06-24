using LurCGT
using SparseArrayKit
using TensorOperations
using LinearAlgebra
using SparseArrays
using Random
using StatsBase
using Test
using TOML
using SQLite
using DBInterface

include("test_utils.jl")
include("sqlite_cleanup.jl")

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
    @test :merge_all_to_global in names(LurCGT)
    @test :merge_table_to_global in names(LurCGT)
    @test :sqlite_stats in names(LurCGT)
    @test :delete_closed_local_sqlite_dbs in names(LurCGT)
    @test :finalize_sqlite! in names(LurCGT)
    @test :finalize_all_sqlite! in names(LurCGT)
    @test hash(Z{3}) == hash((0, 3))
    @test hash(U1) == hash((1,))
    @test hash(SU{3}) == hash((2, 3))
    @test hash(Sp{4}) == hash((3, 4))
    @test hash(SO{5}) == hash((4, 5))
    @test hash(G2) == hash((5,))
    @test LurCGT.irep_cache_key(SU{2}, Int, (1,)) == (SU{2}, Int, (1,))
    @test LurCGT.cgt_cache_key(SU{2}, Int, ((1,), (1,), (0,))) == (SU{2}, Int, ((1,), (1,), (0,)))
    @test LurCGT.fsymbol_cache_key(SU{2}, Int, (1,), (1,), (1,), (1,)) == (SU{2}, Int, (1,), (1,), (1,), (1,))
    @test LurCGT.rsymbol_cache_key(SU{2}, Int, (1,), (1,)) == (SU{2}, Int, (1,), (1,))
    @test LurCGT.xsymbol_cache_key(SU{2}, ((1,),), ((0,),), ((1,),), ((0,),), (1,), (2,)) == (SU{2}, ((1,),), ((0,),), ((1,),), ((0,),), (1,), (2,))
    @test LurCGT.omlist_cache_key(SU{2}, ((1,), (0,)), (1,)) == (SU{2}, ((1,), (0,)), (1,))
    @test LurCGT.validout_cache_key(SU{2}, ((1,), (0,))) == LurCGT.validout_cache_key(SU{2}, ((0,), (1,)))
    @test LurCGT.cgtperm_cache_key(SU{2}, ((1,),), ((0,),), (1, 2)) == (SU{2}, ((1,),), ((0,),), (1, 2))
    @test LurCGT.conjperm_cache_key(SU{2}, ((1,),)) == (SU{2}, ((1,),))
    @test LurCGT.cgtsvd_cache_key(SU{2}, ((1,),), ((0,),), (1,)) == (SU{2}, ((1,),), ((0,),), (1,))
    @test LurCGT.cg3flip_cache_key(SU{2}, Int, ((1,), (1,)), (0,)) == (SU{2}, Int, ((1,), (1,)), (0,))
    cache_hit_key = (:cache_hit_no_string_key,)
    cache_hit_obj = (size_byte=1, value=:hit)
    lock(LurCGT.IREP_CACHE_LOCK) do
        LurCGT.IREP_CACHE[cache_hit_key] = cache_hit_obj
    end
    @test LurCGT._cached_load(
        () -> error("SQLite string key should not be generated on cache hit"),
        SU{2}, "irreps", cache_hit_key, LurCGT.IREP_CACHE, LurCGT.IREP_CACHE_LOCK) === cache_hit_obj
    conjperm = getNsave_Conjperm(SU{2}, ((1,), (1,)))
    @test conjperm isa Conjperm
    @test conjperm === getNsave_Conjperm(SU{2}, ((1,), (1,)))
    @test conjperm.perm == get_conj_perm(get_CGTom(SU{2}, ((1,), (1,)), ((1,), (1,))))
end

@testset "real-valued lowering operators" begin
    symm = (SU{2},)
    weights = ([(2,), (0,), (-2,)],)
    lowering_ops = ([Float64[0 0 0; sqrt(2) 0 0; 0 sqrt(2) 0]],)

    ortho_vecs, space_list = decompose_space(symm, weights, lowering_ops)
    @test size(ortho_vecs) == (3, 3)
    @test haskey(space_list, ((2,),))

    mwirop = Float64[0 sqrt(2) 0; 0 0 sqrt(2); 0 0 0]
    irop, qlabel = get_IROP(symm, weights, lowering_ops, mwirop)
    @test qlabel == ((2,),)
    @test size(irop) == (3, 3, 3)
end

@testset "SQLite environment overrides" begin
    mktempdir() do tmp
        try
            local_dir = joinpath(tmp, "local-db")
            global_dir = joinpath(tmp, "global-db")
            node_local_dir = joinpath(tmp, "node-local-db")
            node_global_dir = joinpath(tmp, "node-global-db")
            empty_local_dir = joinpath(tmp, "empty-local-db")
            empty_global_dir = joinpath(tmp, "empty-global-db")
            fresh_node_local_dir = joinpath(tmp, "fresh-node-local-db")
            fresh_node_global_dir = joinpath(tmp, "fresh-node-global-db")
            process_id = LurCGT.process_local_id()
            symm_name = LurCGT.totxt(SU{2})
            local_copy_key = "sqlite_server_local_copy"
            global_copy_key = "sqlite_server_global_copy"
            payload = UInt8[0x01, 0x02, 0x03]

            withenv("LURCGT_LOCALDB_DIR" => local_dir,
                    "LURCGT_GLOBALDB_DIR" => global_dir) do
                try
                    LurCGT.close_all_sqlite_dbs()
                    local_db = LurCGT.get_sqlite_db(SU{2}, :local)
                    global_db = LurCGT.get_sqlite_db(SU{2}, :global)
                    DBInterface.execute(local_db,
                        "INSERT OR REPLACE INTO irreps (key, data) VALUES (?, ?)",
                        [local_copy_key, payload])
                    DBInterface.execute(global_db,
                        "INSERT OR REPLACE INTO irreps (key, data) VALUES (?, ?)",
                        [global_copy_key, payload])
                finally
                    LurCGT.close_all_sqlite_dbs()
                end
            end

            withenv("LURCGT_RUN_MODE" => "server",
                    "LURCGT_LOCALDB_DIR" => local_dir,
                    "LURCGT_GLOBALDB_DIR" => global_dir,
                    "LURCGT_LOCALDB_DIR_NODE" => node_local_dir,
                    "LURCGT_GLOBALDB_DIR_NODE" => node_global_dir) do
                try
                    LurCGT.close_all_sqlite_dbs()
                    local_node_path = LurCGT.sqlite_db_path(SU{2}, :local; process_id=process_id)
                    global_node_path = LurCGT.sqlite_db_path(SU{2}, :global; process_id=process_id)
                    @test !isfile(local_node_path)
                    @test !isfile(global_node_path)

                    local_db = LurCGT.get_sqlite_db(SU{2}, :local)
                    global_db = LurCGT.get_sqlite_db(SU{2}, :global)

                    @test isfile(local_db.file)
                    @test isfile(global_db.file)
                    @test LurCGT.sqlite_lock_dir() == joinpath(normpath(abspath(node_global_dir)), "locks")
                    @test first(DBInterface.execute(local_db,
                        "SELECT COUNT(*) AS n FROM irreps WHERE key = ?",
                        [local_copy_key])).n == 0
                    @test first(DBInterface.execute(global_db,
                        "SELECT COUNT(*) AS n FROM irreps WHERE key = ?",
                        [global_copy_key])).n == 1
                finally
                    LurCGT.close_all_sqlite_dbs()
                end
            end

            withenv("LURCGT_RUN_MODE" => "server",
                    "LURCGT_LOCALDB_DIR" => empty_local_dir,
                    "LURCGT_GLOBALDB_DIR" => empty_global_dir,
                    "LURCGT_LOCALDB_DIR_NODE" => fresh_node_local_dir,
                    "LURCGT_GLOBALDB_DIR_NODE" => fresh_node_global_dir) do
                try
                    LurCGT.close_all_sqlite_dbs()
                    local_node_path = LurCGT.sqlite_db_path(SU{2}, :local; process_id=process_id)
                    global_node_path = LurCGT.sqlite_db_path(SU{2}, :global; process_id=process_id)
                    @test !isfile(local_node_path)
                    @test !isfile(global_node_path)

                    local_db = LurCGT.get_sqlite_db(SU{2}, :local)
                    global_db = LurCGT.get_sqlite_db(SU{2}, :global)

                    @test isfile(local_db.file)
                    @test isfile(global_db.file)
                    @test dirname(local_db.file) == joinpath(normpath(abspath(fresh_node_local_dir)), symm_name)
                    @test dirname(global_db.file) == normpath(abspath(fresh_node_global_dir))
                    @test first(DBInterface.execute(local_db, "SELECT COUNT(*) AS n FROM irreps")).n == 0
                    @test first(DBInterface.execute(global_db, "SELECT COUNT(*) AS n FROM irreps")).n == 0
                finally
                    LurCGT.close_all_sqlite_dbs()
                end
            end

            withenv("LURCGT_RUN_MODE" => "invalid") do
                @test_throws ArgumentError LurCGT.sqlite_run_mode()
            end

            withenv("LURCGT_RUN_MODE" => "server",
                    "LURCGT_LOCALDB_DIR" => joinpath(tmp, "unused-source-local"),
                    "LURCGT_GLOBALDB_DIR" => joinpath(tmp, "merge-source-global"),
                    "LURCGT_LOCALDB_DIR_NODE" => joinpath(tmp, "merge-node-local"),
                    "LURCGT_GLOBALDB_DIR_NODE" => joinpath(tmp, "merge-node-global")) do
                try
                    LurCGT.close_all_sqlite_dbs()
                    local_db = LurCGT.get_sqlite_db(SU{2}, :local)
                    source_global_path = LurCGT.sqlite_db_source_path(SU{2}, :global)
                    node_global_path = LurCGT.sqlite_db_path(SU{2}, :global)
                    merge_key = "sqlite_server_source_merge"
                    @test !isfile(source_global_path)
                    @test !isfile(node_global_path)

                    DBInterface.execute(local_db,
                        "INSERT OR REPLACE INTO irreps (key, data) VALUES (?, ?)",
                        [merge_key, payload])
                    result = LurCGT.merge_all_to_global(SU{2}; tables=("irreps",), verbose=0)

                    @test result.merged == 1
                    @test isfile(source_global_path)
                    @test !isfile(node_global_path)
                    LurCGT.close_all_sqlite_dbs()
                    source_db = SQLite.DB(source_global_path)
                    try
                        @test first(DBInterface.execute(source_db,
                            "SELECT COUNT(*) AS n FROM irreps WHERE key = ?",
                            [merge_key])).n == 1
                    finally
                        SQLite.close(source_db)
                    end
                finally
                    LurCGT.close_all_sqlite_dbs()
                end
            end

            withenv("LURCGT_LOCALDB_DIR" => joinpath(tmp, "finalize-local"),
                    "LURCGT_GLOBALDB_DIR" => joinpath(tmp, "finalize-global")) do
                try
                    LurCGT.close_all_sqlite_dbs()
                    local_db = LurCGT.get_sqlite_db(SU{2}, :local)
                    local_path = LurCGT.sqlite_db_path(SU{2}, :local)
                    global_path = LurCGT.sqlite_db_path(SU{2}, :global)
                    finalize_key = "sqlite_local_finalize"
                    @test isfile(local_path)
                    @test !isfile(global_path)

                    DBInterface.execute(local_db,
                        "INSERT OR REPLACE INTO irreps (key, data) VALUES (?, ?)",
                        [finalize_key, payload])
                    result = LurCGT.finalize_sqlite!(SU{2}; tables=("irreps",), cleanup=true, verbose=0)

                    @test result.merged == 1
                    @test isfile(global_path)
                    @test !isfile(local_path)
                    global_db = SQLite.DB(global_path)
                    try
                        @test first(DBInterface.execute(global_db,
                            "SELECT COUNT(*) AS n FROM irreps WHERE key = ?",
                            [finalize_key])).n == 1
                    finally
                        SQLite.close(global_db)
                    end
                finally
                    LurCGT.close_all_sqlite_dbs()
                end
            end

            withenv("LURCGT_RUN_MODE" => "server",
                    "LURCGT_GLOBALDB_DIR" => joinpath(tmp, "finalize-source-global"),
                    "LURCGT_LOCALDB_DIR_NODE" => joinpath(tmp, "finalize-node-local"),
                    "LURCGT_GLOBALDB_DIR_NODE" => joinpath(tmp, "finalize-node-global")) do
                try
                    LurCGT.close_all_sqlite_dbs()
                    source_global_path = LurCGT.sqlite_db_source_path(SU{2}, :global)
                    node_global_path = LurCGT.sqlite_db_path(SU{2}, :global)
                    local_path = LurCGT.sqlite_db_path(SU{2}, :local)
                    finalize_key = "sqlite_server_finalize"

                    mkpath(dirname(source_global_path))
                    source_db = SQLite.DB(source_global_path)
                    try
                        LurCGT._init_sqlite_db(source_db)
                        DBInterface.execute(source_db,
                            "INSERT OR REPLACE INTO irreps (key, data) VALUES (?, ?)",
                            [global_copy_key, payload])
                    finally
                        SQLite.close(source_db)
                    end

                    global_db = LurCGT.get_sqlite_db(SU{2}, :global)
                    @test isfile(node_global_path)
                    @test first(DBInterface.execute(global_db,
                        "SELECT COUNT(*) AS n FROM irreps WHERE key = ?",
                        [global_copy_key])).n == 1

                    local_db = LurCGT.get_sqlite_db(SU{2}, :local)
                    DBInterface.execute(local_db,
                        "INSERT OR REPLACE INTO irreps (key, data) VALUES (?, ?)",
                        [finalize_key, payload])
                    result = LurCGT.finalize_sqlite!(SU{2}; tables=("irreps",), cleanup=true, verbose=0)

                    @test result.merged == 1
                    @test isfile(source_global_path)
                    @test !isfile(local_path)
                    @test !isfile(node_global_path)
                    source_db_check = SQLite.DB(source_global_path)
                    try
                        @test first(DBInterface.execute(source_db_check,
                            "SELECT COUNT(*) AS n FROM irreps WHERE key = ?",
                            [finalize_key])).n == 1
                    finally
                        SQLite.close(source_db_check)
                    end
                finally
                    LurCGT.close_all_sqlite_dbs()
                end
            end
        finally
            LurCGT.close_all_sqlite_dbs()
        end
    end
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
