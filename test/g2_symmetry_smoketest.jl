using Test
using LurCGT

@testset "G2 symmetry basics" begin
    @test G2 == LurCGT.G2
    @test LurCGT.isvalidsymm(G2)
    @test LurCGT.totxt(G2) == "G2"
    @test LurCGT.defirepdim(G2) == 7
    @test LurCGT.nlops(G2) == 2
    @test LurCGT.nzops(G2) == 2

    @test LurCGT.getsz_def(G2) == Dict(
        (1, 1) => (1, 1),
        (-1, 1) => (2, 2),
        (2, 0) => (3, 3),
        (0, 0) => (4, 4),
        (-2, 0) => (5, 5),
        (1, -1) => (6, 6),
        (-1, -1) => (7, 7),
    )
    @test LurCGT.getsz_def_vec(G2) == [
        [1, 1], [-1, 1], [2, 0], [0, 0], [-2, 0], [1, -1], [-1, -1],
    ]
    @test Tuple(LurCGT.getdz(G2, 1)) == (2, 0)
    @test Tuple(LurCGT.getdz(G2, 2)) == (-3, 1)

    @test LurCGT.charlist(G2) == [1, 2, 3, 0, -3, -2, -1]
    @test LurCGT.crystal_chars_map(G2) == Dict(1 => 1, 2 => 2, 3 => 3, 0 => 4, -3 => 5, -2 => 6, -1 => 7)
    @test LurCGT.get_fops_std(G2) == [
        Dict(1 => 2, 3 => 0, 0 => -3, -2 => -1),
        Dict(2 => 3, -3 => -2),
    ]

    @test LurCGT.qlab2mwz(G2, (1, 0)) == (1, 1)
    @test LurCGT.qlab2mwz(G2, (0, 1)) == (0, 2)
    @test LurCGT.getqlabel(G2, (1, 1)) == (1, 0)
    @test LurCGT.getqlabel(G2, (0, 2)) == (0, 1)

    cops = LurCGT.get_crystal_ops(G2)
    fund2 = LurCGT.mw_tableau_qlabel(G2, (0, 1))
    @test isnothing(LurCGT.apply_lowering(fund2, cops, 1))
    @test LurCGT.apply_lowering(fund2, cops, 2).colread == [3, 1]

    tab = LurCGT.mw_tableau_qlabel(G2, (2, 1))
    @test tab.shape == [3, 1]
    @test tab.colread == [2, 1, 1, 1]
end

@testset "G2 alternating square matrices" begin
    defrep = LurCGT.getdefirep(G2, BigInt)
    sl_alt, _, mw_alt, _ = LurCGT.def_altprod(defrep, 2)

    @test mw_alt == (0, 2)
    @test Matrix(sl_alt[1][(1, 1)]) == BigInt[
        2 0
        1 1
    ]
end

@testset "G2 alternating square norms" begin
    defrep = LurCGT.getdefirep(G2, BigInt)
    custom_inner = copy(defrep.innerprod)
    custom_inner[(1, 1)] = BigInt[2;;]
    custom_inner[(-1, 1)] = BigInt[3;;]
    custom_inner[(2, 0)] = BigInt[5;;]
    custom_inner[(0, 0)] = BigInt[7;;]
    custom_inner[(-2, 0)] = BigInt[11;;]
    custom_inner[(1, -1)] = BigInt[13;;]
    custom_inner[(-1, -1)] = BigInt[17;;]

    custom_rep = LurCGT.Irep(
        G2,
        BigInt,
        defrep.Sl,
        defrep.Sz,
        custom_inner,
        defrep.inv_innerprod,
        LurCGT.dimension(defrep),
    )

    _, _, _, alt_inner = LurCGT.def_altprod(custom_rep, 2)
    @test alt_inner[(1, 1)] == BigInt[
        14 0
        0 15
    ]
end

g2_expected_dims() = [
    1, 7, 14, 27, 64, 77, 77, 182, 189, 273, 286, 378, 448, 714, 729, 748,
    896, 924, 1254, 1547, 1728, 1729, 2079, 2079, 2261, 2926, 3003, 3289,
    3542, 4096, 4914, 4928, 4928, 5005, 5103, 6630, 7293, 7371, 7722, 8372,
    9177, 9660, 10206, 10556, 11571, 11648, 12096, 13090,
]

g2_weyl_dim(a::Int, b::Int) =
    ((a + 1) * (b + 1) * (a + b + 2) * (a + 2b + 3) * (a + 3b + 4) * (2a + 3b + 5)) ÷ 120

@testset "G2 irrep dimensions" begin
    @test LurCGT.dimension(LurCGT.getNsave_irep(G2, BigInt, (0, 0))) == 1
    @test LurCGT.dimension(LurCGT.getNsave_irep(G2, BigInt, (1, 0))) == 7
    @test LurCGT.dimension(LurCGT.getNsave_irep(G2, BigInt, (0, 1))) == 14

    dims = Int[]
    for a in 0:12, b in 0:12
        g2_weyl_dim(a, b) <= 13090 || continue
        println("checking G2 qlabel ", (a, b))
        push!(dims, LurCGT.dimension(LurCGT.getNsave_irep(G2, BigInt, (a, b))))
    end
    @test sort(dims) == g2_expected_dims()
end
