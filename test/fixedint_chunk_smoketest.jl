using Test
using Serialization
using LurCGT

@testset "fixedint dimension chunks" begin
    @test LurCGT.fixedint_dimension_chunks(1, 10, 4) == [(1, 2), (3, 5), (6, 7), (8, 10)]
    @test LurCGT.fixedint_dimension_chunks(5, 5, 3) == [(5, 5)]
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
