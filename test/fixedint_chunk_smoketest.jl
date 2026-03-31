using Test
using Serialization
using LurCGT

include("fixedint_catalog_read_driver.jl")

@testset "fixedint dimension chunks" begin
    @test LurCGT.fixedint_dimension_chunks(1, 10, 4) == [(1, 2), (3, 5), (6, 7), (8, 10)]
    @test LurCGT.fixedint_dimension_chunks(5, 5, 3) == [(5, 5)]
end

@testset "fixedint symmetry parser accepts G2" begin
    @test parse_nonabelian_symmetry("G2") == G2
end

@testset "fixedint canonical pairs" begin
    entries1 = [
        (qlabel=(0, 1), dim=3),
        (qlabel=(1, 0), dim=3),
        (qlabel=(2, 0), dim=6),
    ]
    entries2 = [
        (qlabel=(0, 1), dim=3),
        (qlabel=(1, 0), dim=3),
        (qlabel=(0, 2), dim=6),
    ]
    pairs = LurCGT.fixedint_canonical_pairs(entries1, entries2)
    @test pairs == [
        (q1=(0, 1), dim1=3, q2=(0, 1), dim2=3),
        (q1=(0, 1), dim1=3, q2=(1, 0), dim2=3),
        (q1=(1, 0), dim1=3, q2=(1, 0), dim2=3),
        (q1=(0, 1), dim1=3, q2=(0, 2), dim2=6),
        (q1=(1, 0), dim1=3, q2=(0, 2), dim2=6),
    ]
    @test !any(pair -> pair.q1 == (2, 0) && pair.q2 == (0, 2), pairs)
end

@testset "fixedint chunk path tokens" begin
    mktempdir() do tmp
        path = LurCGT.fixedint_chunk_result_path(SU{2}, Int64, (1, 2), (3, 4); base_dir=tmp)
        @test endswith(path, joinpath("Int64", "SU2", "chunks", "1_2__3_4.jls"))
    end
end

@testset "fixedint irrep merge helper targets irreps table" begin
    seen = Ref{Any}(nothing)
    result = LurCGT.merge_fixedint_ireps_to_global(
        SU{2};
        clear_local_after=false,
        verbose=0,
        merge_fn=(S, table_name; clear_local_after=true, verbose=1) -> begin
            seen[] = (
                symmetry=S,
                table_name=table_name,
                clear_local_after=clear_local_after,
                verbose=verbose,
            )
            return (merged=3, skipped=1)
        end,
    )

    @test seen[] == (
        symmetry=SU{2},
        table_name="irreps",
        clear_local_after=false,
        verbose=0,
    )
    @test result == (merged=3, skipped=1)
end

@testset "fixedint catalog and chunk smoke" begin
    mktempdir() do tmp
        for RT in (Int64, Int128)
            catalog = LurCGT.update_fixedint_irrep_catalog(SU{2}, RT; maxdim=4, base_dir=tmp, save=true)
            @test catalog.numtype == LurCGT.totxt(RT)
            @test catalog.symmetry == LurCGT.totxt(SU{2})
            @test catalog.accepted == [
                (qlabel=(0,), dim=1),
                (qlabel=(1,), dim=2),
                (qlabel=(2,), dim=3),
                (qlabel=(3,), dim=4),
            ]

            summary = LurCGT.run_fixedint_cgt_chunk(
                SU{2},
                RT,
                1,
                4,
                1,
                4,
                2,
                2,
                1,
                2;
                base_dir=tmp,
                save=true,
                verbose=0,
            )

            @test summary.failed_pairs == 0
            @test summary.dim_range1 == (1, 2)
            @test summary.dim_range2 == (3, 4)

            chunk_path = LurCGT.fixedint_chunk_result_path(
                SU{2},
                RT,
                summary.dim_range1,
                summary.dim_range2;
                base_dir=tmp,
            )
            @test isfile(chunk_path)

            cells = LurCGT.collect_fixedint_plot_cells(SU{2}, RT; base_dir=tmp)
            @test any(cell -> cell.dim1 == 1 && cell.dim2 == 3 && cell.status == :passed, cells)
        end
    end
end

@testset "fixedint catalog supports G2" begin
    mktempdir() do tmp
        catalog = LurCGT.update_fixedint_irrep_catalog(G2, Int64; maxdim=14, base_dir=tmp, save=false)
        @test catalog.numtype == "Int64"
        @test catalog.symmetry == "G2"
        @test catalog.accepted[1:3] == [
            (qlabel=(0, 0), dim=1),
            (qlabel=(1, 0), dim=7),
            (qlabel=(0, 1), dim=14),
        ]
    end
end

@testset "fixedint chunk canonicalizes qlabels for CGT generation" begin
    mktempdir() do tmp
        LurCGT.update_fixedint_irrep_catalog(G2, Int64; maxdim=14, base_dir=tmp, save=true)
        seen_pairs = Tuple{Tuple{Int, Int}, Tuple{Int, Int}}[]

        summary = LurCGT.run_fixedint_cgt_chunk(
            G2,
            Int64,
            1,
            14,
            1,
            14,
            1,
            1,
            1,
            1;
            base_dir=tmp,
            save=false,
            update_catalog=false,
            merge_local_ireps=false,
            generate_cgt_fn=(S, CT1, CT2, qpair, out; assertlev=1, save=false) -> begin
                push!(seen_pairs, qpair)
                return nothing
            end,
            verbose=0,
        )

        @test summary.failed_pairs == 0
        @test ((0, 1), (1, 0)) in seen_pairs
        @test !(((1, 0), (0, 1)) in seen_pairs)
        @test all(qpair -> qpair[1] <= qpair[2], seen_pairs)
    end
end

@testset "fixedint catalog can merge local irreps after search" begin
    mktempdir() do tmp
        seen = Ref{Any}(nothing)
        catalog = LurCGT.update_fixedint_irrep_catalog(
            SU{2},
            Int64;
            maxdim=3,
            base_dir=tmp,
            save=false,
            merge_local_ireps=true,
            merge_ireps_fn=(S; clear_local_after=true, verbose=1) -> begin
                seen[] = (
                    symmetry=S,
                    clear_local_after=clear_local_after,
                    verbose=verbose,
                )
                return (merged=2, skipped=0)
            end,
            verbose=0,
        )

        @test catalog.accepted == [
            (qlabel=(0,), dim=1),
            (qlabel=(1,), dim=2),
            (qlabel=(2,), dim=3),
        ]
        @test seen[] == (
            symmetry=SU{2},
            clear_local_after=true,
            verbose=0,
        )
    end
end

@testset "fixedint chunk requires prebuilt catalog in read-only mode" begin
    mktempdir() do tmp
        err = try
            LurCGT.run_fixedint_cgt_chunk(
                SU{2},
                Int64,
                1,
                4,
                1,
                4,
                2,
                2,
                1,
                2;
                base_dir=tmp,
                save=false,
                update_catalog=false,
                verbose=0,
            )
            nothing
        catch caught
            caught
        end

        @test err isa ArgumentError
        @test occursin("Run test/fixedint_irep_catalog_driver.jl first.", sprint(showerror, err))
    end
end

@testset "fixedint chunk reads catalog without updating or merging" begin
    mktempdir() do tmp
        LurCGT.update_fixedint_irrep_catalog(SU{2}, Int64; maxdim=4, base_dir=tmp, save=true)

        summary = LurCGT.run_fixedint_cgt_chunk(
            SU{2},
            Int64,
            1,
            4,
            1,
            4,
            2,
            2,
            1,
            2;
            base_dir=tmp,
            save=false,
            update_catalog=false,
            merge_local_ireps=false,
            merge_ireps_fn=(S; clear_local_after=true, verbose=1) -> error("merge should not run for chunk workers"),
            verbose=0,
        )

        @test summary.failed_pairs == 0
        @test summary.dim_range1 == (1, 2)
        @test summary.dim_range2 == (3, 4)
    end
end

@testset "fixedint plot cells stay separated by numtype" begin
    mktempdir() do tmp
        summary64 = (
            numtype=LurCGT.totxt(Int64),
            symmetry=LurCGT.totxt(SU{2}),
            dim_range1=(1, 4),
            dim_range2=(1, 4),
            total_pairs=1,
            passed_pairs=1,
            failed_pairs=0,
            results=[(q1=(0,), dim1=1, q2=(1,), dim2=2, status=:passed, reason="passed")],
        )
        summary128 = (
            numtype=LurCGT.totxt(Int128),
            symmetry=LurCGT.totxt(SU{2}),
            dim_range1=(1, 4),
            dim_range2=(1, 4),
            total_pairs=1,
            passed_pairs=1,
            failed_pairs=0,
            results=[(q1=(0,), dim1=1, q2=(2,), dim2=3, status=:passed, reason="passed")],
        )

        serialize(LurCGT.fixedint_chunk_result_path(SU{2}, Int64, (1, 4), (1, 4); base_dir=tmp), summary64)
        serialize(LurCGT.fixedint_chunk_result_path(SU{2}, Int128, (1, 4), (1, 4); base_dir=tmp), summary128)

        @test LurCGT.collect_fixedint_plot_cells(SU{2}, Int64; base_dir=tmp) == [(dim1=1, dim2=2, status=:passed)]
        @test LurCGT.collect_fixedint_plot_cells(SU{2}, Int128; base_dir=tmp) == [(dim1=1, dim2=3, status=:passed)]
    end
end

@testset "fixedint catalog read driver prints accepted entries" begin
    mktempdir() do tmp
        LurCGT.update_fixedint_irrep_catalog(SU{2}, Int64; maxdim=3, base_dir=tmp, save=true)

        output = withenv("FIXEDINT_DATA_ROOT" => tmp) do
            sprint() do io
                run_fixedint_catalog_read_cli("Int64", "SU2"; mode="accepted", io=io)
            end
        end

        @test occursin("Path:", output)
        @test occursin("Accepted: 3", output)
        @test occursin("accepted dim=1 q=(0,)", output)
        @test occursin("accepted dim=3 q=(2,)", output)
    end
end
