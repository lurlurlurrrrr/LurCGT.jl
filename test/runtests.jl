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
    expected_exports = Set(Symbol[
        :LurCGT, :Z, :U1, :SU, :SO, :Sp, :G2,
        :Symmetry, :AbelianSymm, :NonabelianSymm,
        :add_qn, :decompose_irop, :decompose_space, :dimension,
        :get_CGTom, :get_IROP, :get_conj_perm, :get_dualq,
        :getNsave_CGTperm, :getNsave_CGTSVD, :getNsave_CGTQR,
        :getNsave_Xsymbol, :getNsave_Conjperm,
        :getNsave_omlist, :getNsave_validout, :isabelian,
        :nzops, :remove_zeros, :totxt, :transf_basis!, :to_float,
    ])
    @test Set(names(LurCGT)) == expected_exports
    @test isabelian(U1)
    exact_cgt = LurCGT.getNsave_cg3(SU{2}, BigInt, ((1,), (1,)), [(0,)])[(0,)]
    cgt, qlabels, directions = to_float(exact_cgt, Float64)
    default_cgt, _, _ = to_float(exact_cgt)
    @test size(cgt) == (2, 2, 1, 1)
    @test eltype(default_cgt) == Float64
    @test qlabels == ((1,), (1,), (0,))
    @test directions == ('+', '+', '-')
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
        LurCGT.IREP_CACHE[(LurCGT.CRYSTAL_BASIS_VERSION, cache_hit_key)] = cache_hit_obj
    end
    @test LurCGT._cached_load(
        () -> error("SQLite string key should not be generated on cache hit"),
        SU{2}, "irreps", cache_hit_key, LurCGT.IREP_CACHE, LurCGT.IREP_CACHE_LOCK) === cache_hit_obj
    conjperm = getNsave_Conjperm(SU{2}, ((1,), (1,)))
    @test conjperm isa LurCGT.Conjperm
    @test conjperm === getNsave_Conjperm(SU{2}, ((1,), (1,)))
    @test conjperm.perm == get_conj_perm(get_CGTom(SU{2}, ((1,), (1,)), ((1,), (1,))))
end

@testset "deferred virtual-crystal irrep backends" begin
    for (S,q) in ((SO{5},(1,0)),(SO{7},(0,0,1)),(G2,(1,0)))
        @test isempty(LurCGT.fundamental_recipes(S))
        @test_throws ArgumentError LurCGT.getNsave_irep(S,BigInt,q)
    end
end

@testset "exceptional scaled-root coordinates" begin
    for S in (E6,E7)
        r = nzops(S); d = S == E6 ? 3 : 2
        A = LurCGT._exceptional_cartan(S)
        metric_numerator = reverse(reverse(A; dims=1); dims=2)
        steps = LurCGT.getdzs(S)
        labels = [LurCGT._unit_label(S,k) for k in 1:r]
        push!(labels,ntuple(i -> i,r))
        for q in labels
            z = LurCGT.qlab2mwz(S,q)
            @test LurCGT.getqlabel(S,z) == q
            for i in 1:r
                lowered = z .- Tuple(steps[i])
                @test LurCGT.getqlabel(S,lowered) == q .- Tuple(A[:,i])
                @test LurCGT.rev_less(lowered,z)
            end
        end
        for i in 1:r
            @test steps[i] == [j == r+1-i ? d : 0 for j in 1:r]
            @test LurCGT.qlab2mwz(S,Tuple(A[:,i])) == Tuple(steps[i])
            for j in 1:r
                @test dot(steps[i],metric_numerator*steps[j]) == d^2*A[i,j]
            end
        end
        # Coordinate-unit steps for adjacent nodes are not physically orthogonal.
        @test dot(steps[1],metric_numerator*steps[3]) == -d^2
    end
end

@testset "minuscule seeds and integer crystal words" begin
    for (S,nodes,dims) in ((SU{3},[1,2],[3,3]),(Sp{6},[1],[6]),
                           (SO{4},[1,2],[2,2]),(SO{6},[1,2,3],[6,4,4]),
                           (SO{8},[1,3,4],[8,8,8]),(E6,[1,6],[27,27]),(E7,[7],[56]))
        symbols = Set{Int}()
        for (k,d) in zip(nodes,dims)
            q = LurCGT._unit_label(S,k)
            rep = LurCGT.get_minuscule_irep(S,BigInt,k)
            @test rep.qlabel == q
            @test dimension(rep) == d
            highest = LurCGT.mw_crystal_qlabel(S,q)
            @test length(highest.word) == 1
            vertices = LurCGT._component_vertices(S,highest)
            @test length(vertices) == d
            for v in vertices
                @test length(v.word) == 1
                @test !(only(v.word) in symbols)
                push!(symbols,only(v.word))
                @test haskey(rep.Sz,LurCGT.get_weight(v))
            end
        end
    end
end

@testset "fundamental recipe metadata and recursive requests" begin
    for S in (SU{3},Sp{6},SO{4},SO{6},SO{8},E6,E7)
        recipes = LurCGT.fundamental_recipes(S)
        @test recipes isa Vector{LurCGT.FundamentalRecipe}
        @test LurCGT.validate_recipe_catalogue(S,recipes) === nothing
        @test LurCGT.irrep_validation_rules(S) isa Vector{LurCGT.IrrepValidationRule}
    end
    bad = LurCGT.fundamental_recipes(SO{8})
    bad[2] = LurCGT.ProductRecipe([0,1,0,0],([1,0,0,0],[1,0,0,0]),0,28)
    @test_throws ArgumentError LurCGT.validate_recipe_catalogue(SO{8},bad)
    bad[2] = LurCGT.ProductRecipe([0,1,0,0],([1,0,0],[1,0,0,0]),1,28)
    @test_throws ArgumentError LurCGT.validate_recipe_catalogue(SO{8},bad)
    bad[2] = LurCGT.ProductRecipe([0,1,0,0],([-1,0,0,0],[1,0,0,0]),1,28)
    @test_throws ArgumentError LurCGT.validate_recipe_catalogue(SO{8},bad)
    @test_throws ArgumentError LurCGT.construct_recipe(SU{3},BigInt,LurCGT.AdjointRecipe([1,0],3))
    mktempdir() do tmp
        withenv("LURCGT_LOCALDB_DIR"=>joinpath(tmp,"local"),
                "LURCGT_GLOBALDB_DIR"=>joinpath(tmp,"global")) do
            LurCGT.close_all_sqlite_dbs()
            lock(LurCGT.IREP_CACHE_LOCK) do
                empty!(LurCGT.IREP_CACHE)
            end
            try
                @test dimension(LurCGT.getNsave_irep(SO{8},BigInt,(0,1,0,0))) == 28
                @test dimension(LurCGT.load_irep_sqlite(SO{8},BigInt,(1,0,0,0))) == 8
                @test isnothing(LurCGT.load_irep_sqlite(SO{8},BigInt,(0,0,1,0)))
                @test isnothing(LurCGT.load_irep_sqlite(SO{8},BigInt,(0,0,0,1)))
            finally
                LurCGT.close_all_sqlite_dbs()
                lock(LurCGT.IREP_CACHE_LOCK) do
                    empty!(LurCGT.IREP_CACHE)
                end
            end
        end
    end
end

@testset "canonical nonfundamental copies after fundamental bootstrap" begin
    for S in (Sp{4},Sp{6},SO{8})
        @test LurCGT.validate_fundamental_ireps(S) > 0
    end
    @test dimension(LurCGT.getNsave_irep(E6,BigInt,(0,1,0,0,0,0))) == 78
    @test dimension(LurCGT.getNsave_irep(E7,BigInt,(1,0,0,0,0,0,0))) == 133
    rules = LurCGT.irrep_validation_rules(SO{8})
    @test LurCGT.validate_fundamental_products(SO{8}; rules) == length(rules)
    rep = LurCGT.getNsave_irep(SU{3},BigInt,(1,1))
    Sl,Sz,_,wi = LurCGT.tensor_prod_reps(rep,rep)
    mw = LurCGT.qlab2mwz(SU{3},(1,1))
    repeated = LurCGT.IrrepValidationRule([1,1],([1,1],[1,1]))
    @test LurCGT.validate_fundamental_products(SU{3};rules=[repeated]) == 2
    absent = LurCGT.IrrepValidationRule([1,0],([1,0],[1,0]))
    @test_throws ArgumentError LurCGT.validate_fundamental_products(SU{3};rules=[absent])
    for copy_index in 1:2
        other = LurCGT.get_irep_frommw(SU{3},Sl,Sz,rep.innerprod,rep.innerprod,mw,wi;copy_index)
        @test LurCGT.same_irep_fields(rep,other)
    end
    @test_throws ArgumentError LurCGT.get_irep_frommw(SU{3},Sl,Sz,rep.innerprod,rep.innerprod,mw,wi;copy_index=3)
    @test_throws ArgumentError LurCGT.get_irep_frommw(SU{3},Sl,Sz,rep.innerprod,rep.innerprod,(100,100),wi)
    @test_throws ErrorException LurCGT._exact_divide(BigInt(1),2)
    @test all(gcd(B)==1 for (B,_) in values(rep.inv_innerprod))
    lock(LurCGT.IREP_CACHE_LOCK) do
        delete!(LurCGT.IREP_CACHE,(LurCGT.CRYSTAL_BASIS_VERSION,
            LurCGT.irep_cache_key(SU{3},BigInt,(1,1))))
    end
    @test LurCGT.same_irep_fields(rep,LurCGT.load_irep_sqlite(SU{3},BigInt,(1,1)))
end

@testset "highest tensor words" begin
    for (S, q) in ((SU{3}, (2,2)), (Sp{6}, (1,1,1)),
                    (SO{8}, (1,1,2,0)), (SO{8}, (1,1,0,2)),
                    (SO{8}, (0,0,1,1)))
        word = LurCGT.mw_crystal_qlabel(S,q)
        @test word.word isa Vector{Int}
        @test LurCGT.get_weight(word) == LurCGT.qlab2mwz(S,q)
    end
    @test LurCGT.mw_crystal_qlabel(G2, (2, 1)).word == [1, 1, 2, 1]
    @test LurCGT.mw_crystal_qlabel(SO{7}, (1, 1, 2)).word == [1, 2, 1, 3, 2, 1]
    @test_throws AssertionError LurCGT.mw_crystal_qlabel(SO{7}, (0, 0, 1))
end

@testset "SO half-spin labels round trip" begin
    for q in ((0, 0, 1, 0), (0, 0, 0, 1))
        @test LurCGT.getqlabel(SO{8}, LurCGT.qlab2mwz(SO{8}, q)) == q
    end
end

@testset "chapter 2 tensor product rule" begin
    for S in (SU{3}, SO{5}, Sp{4}, SO{8}, G2)
        cops = LurCGT.get_crystal_ops(S)
        # All two-factor states, including zero symbols and local strings of length two.
        for x in LurCGT.charlist(S), y in LurCGT.charlist(S), i in 1:LurCGT.nlops(S)
            vertex = LurCGT.CrystalVertex{S}([x, y])
            ex, ey = get(cops.ϵ, (i, x), 0), get(cops.ϵ, (i, y), 0)
            px, py = get(cops.ϕ, (i, x), 0), get(cops.ϕ, (i, y), 0)
            @test LurCGT.epsilon(vertex, cops, i) == ey + max(0, ex - py)
            @test LurCGT.phi(vertex, cops, i) == px + max(0, py - ex)
            for (op, local_ops, index) in (
                (LurCGT.apply_lowering, cops.f[i], py <= ex ? 1 : 2),
                (LurCGT.apply_raising, cops.e[i], py < ex ? 1 : 2),
            )
                result = op(vertex, cops, i)
                if haskey(local_ops, vertex.word[index])
                    expected = copy(vertex.word)
                    expected[index] = local_ops[expected[index]]
                    @test result.word == expected
                else
                    @test result === nothing
                end
            end
            @test vertex.word == [x, y]
        end
    end
end

@testset "longest Weyl words" begin
    for (S, r, count_roots) in ((SU{5}, 4, 10), (Sp{8}, 4, 16),
                               (SO{9}, 4, 16), (SO{8}, 4, 12),
                               (SO{6}, 3, 6), (SO{4}, 2, 2), (G2, 2, 6))
        word = LurCGT.crystal_lowering_word(S)
        @test length(word) == count_roots
        @test all(i -> 1 <= i <= r, word)
        # Reflect simple roots in their Euclidean realization; w0 makes them negative.
        if S != G2
            n = S == SU{5} ? r + 1 : r
            roots = [zeros(Int, n) for _ in 1:r]
            for i in 1:r-1
                roots[i][i] = 1
                roots[i][i+1] = -1
            end
            roots[r][r] = 1
            if S == SU{5}
                roots[r][r+1] = -1
            elseif S in (SO{4}, SO{6}, SO{8})
                roots[r][r-1] = 1
            elseif S == Sp{8}
                roots[r][r] = 2
            end
            for root in roots
                v = copy(root)
                for i in word
                    a = roots[i]
                    v -= div(2dot(v, a), dot(a, a)) * a
                end
                @test first(filter(!iszero, v)) < 0
            end
        end
    end
end

@testset "crystal component strings and published words" begin
    for (S, q, dim) in (
        (SU{3}, (1, 1), 8), (SU{3}, (2, 1), 15),
        (SO{5}, (1, 2), 35), (Sp{4}, (1, 1), 16),
        (SO{8}, (0, 1, 0, 0), 28), (SO{8}, (0, 0, 1, 1), 56),
        (SO{8}, (0, 0, 2, 0), 35), (SO{8}, (0, 0, 0, 2), 35),
        (G2, (0, 1), 14), (G2, (1, 1), 64),
    )
        cops = LurCGT.get_crystal_ops(S)
        highest = LurCGT.mw_crystal_qlabel(S, q)
        @test LurCGT.getqlabel(S, LurCGT.get_weight(highest)) == q
        @test all(i -> LurCGT.epsilon(highest, cops, i) == 0, 1:LurCGT.nlops(S))
        seen = Dict(highest => 1)
        vertices = [highest]
        for vertex in vertices
            saved_word, saved_hash = copy(vertex.word), hash(vertex)
            for i in 1:LurCGT.nlops(S)
                eps, ph = LurCGT.epsilon(vertex, cops, i), LurCGT.phi(vertex, cops, i)
                for (op, count) in ((LurCGT.apply_raising, eps), (LurCGT.apply_lowering, ph))
                    current = vertex
                    for _ in 1:count
                        current = op(current, cops, i)
                        @test current !== nothing
                    end
                    @test op(current, cops, i) === nothing
                end
                lowered = LurCGT.apply_lowering(vertex, cops, i)
                if lowered !== nothing
                    @test LurCGT.apply_raising(lowered, cops, i) == vertex
                    @test LurCGT.get_weight(lowered) == LurCGT.get_weight(vertex) .- LurCGT.getdz(S, i)
                    @test LurCGT.epsilon(lowered, cops, i) == eps + 1
                    @test LurCGT.phi(lowered, cops, i) == ph - 1
                    if !haskey(seen, lowered)
                        seen[lowered] = length(vertices) + 1
                        push!(vertices, lowered)
                    end
                end
                scratch = LurCGT.CrystalVertex{S}(copy(vertex.word))
                @test LurCGT.apply_lowering!(scratch, cops, i) == (lowered !== nothing)
                @test scratch == (lowered === nothing ? vertex : lowered)
                raised = LurCGT.apply_raising(vertex, cops, i)
                scratch = LurCGT.CrystalVertex{S}(copy(vertex.word))
                @test LurCGT.apply_raising!(scratch, cops, i) == (raised !== nothing)
                @test scratch == (raised === nothing ? vertex : raised)
            end
            @test vertex.word == saved_word
            @test hash(vertex) == saved_hash
            @test seen[LurCGT.CrystalVertex{S}(copy(vertex.word))] == seen[vertex]
        end
        @test length(vertices) == dim
        staged = [highest]
        staged_seen = Set(staged)
        for i in LurCGT.crystal_lowering_word(S)
            nsource = length(staged)
            for k in 1:nsource
                current = staged[k]
                while true
                    current = LurCGT.apply_lowering(current, cops, i)
                    current === nothing && break
                    current in staged_seen && break
                    push!(staged_seen, current)
                    push!(staged, current)
                end
            end
        end
        @test staged_seen == Set(vertices)
        if S == SU{3} && q == (1, 1)
            @test count(v -> LurCGT.get_weight(v) == (0, 0), vertices) == 2
        end
    end
end

@testset "mixed highest-weight numerical construction" begin
    for (S, q, dim) in ((SU{3}, (2, 1), 15),
                       (Sp{4}, (1, 1), 16), (SO{8}, (1, 1, 0, 0), 160))
        rep = LurCGT.getNsave_irep(S, BigInt, q)
        @test LurCGT.dimension(rep) == dim
        @test rep.qlabel == q
        @test sum(size(g, 1) for g in values(rep.innerprod)) == dim
    end
end

@testset "crystal basis cache isolation" begin
    mktempdir() do tmp
        env = Dict("LURCGT_LOCALDB_DIR" => joinpath(tmp, "local"),
                   "LURCGT_GLOBALDB_DIR" => joinpath(tmp, "global"))
        @test LurCGT.sqlite_db_path(SU{3}, :global; env) != joinpath(tmp, "global", "SU3.db")
        @test endswith(LurCGT.sqlite_db_path(SU{3}, :local; env, process_id="test"),
                       "test_$(LurCGT.CRYSTAL_BASIS_VERSION).db")
        @test LurCGT.sqlite_db_source_path(SU{3}, :global; env) ==
              LurCGT.sqlite_db_path(SU{3}, :global; env) ==
              LurCGT.sqlite_merge_global_db_path(SU{3}; env)
        withenv("LURCGT_LOCALDB_DIR" => env["LURCGT_LOCALDB_DIR"],
                "LURCGT_GLOBALDB_DIR" => env["LURCGT_GLOBALDB_DIR"]) do
            try
                LurCGT.close_all_sqlite_dbs()
                LurCGT.clear_sqlite_cache!()
                key = LurCGT.irep_cache_key(SU{3}, BigInt, (1, 1))
                legacy = (size_byte=1, value=:legacy)
                LurCGT.IREP_CACHE[key] = legacy
                @test LurCGT.load_irep_sqlite(SU{3}, BigInt, (1, 1)) === nothing
                rep = LurCGT.getNsave_irep(SU{3}, BigInt, (1, 1))
                @test LurCGT.dimension(rep) == 8
                @test LurCGT.IREP_CACHE[(LurCGT.CRYSTAL_BASIS_VERSION, key)] === rep
                LurCGT.clear_sqlite_cache!()
                loaded = LurCGT.load_irep_sqlite(SU{3}, BigInt, (1, 1))
                @test loaded.Sl == rep.Sl
                @test loaded.innerprod == rep.innerprod
            finally
                LurCGT.close_all_sqlite_dbs()
                LurCGT.clear_sqlite_cache!()
            end
        end
    end
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

@testset "CGTQR tests" begin
    test_CGTQR_unitary_examples(SU{2}; verbose=0)
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
        rep = LurCGT.getNsave_irep(Sp{6}, BigInt, qlabel)
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
    test_Xsym(SU{3}, 5000, 3, 20; verbose=0)
    test_Xsym(SU{2}, 500, 8, 100; verbose=0)
    LurCGT.merge_all_to_global(SU{2})
    LurCGT.merge_all_to_global(SU{3})
end
