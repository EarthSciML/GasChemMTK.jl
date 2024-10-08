using GasChem
using DifferentialEquations, ModelingToolkit, DynamicQuantities

tspan = (0.0, 360.0)

@testset "Base case" begin
    answer = 17.204446308176035
    rs = structural_simplify(SuperFast())
    sol = solve(ODEProblem(rs, [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)

    @test sol[rs.O3][end] ≈ answer
end

@testset "DMS sensitivity" begin
    u_dms = 0.8231258833472341

    rs1 = structural_simplify(SuperFast())
    o1 = solve(ODEProblem(rs1, [rs1.DMS => 76], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    rs2 = structural_simplify(SuperFast())
    o2 = solve(ODEProblem(rs2, [rs2.DMS => 46], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    test1 = o1[rs1.SO2][end] - o2[rs2.SO2][end]

    @test test1 ≈ u_dms
end

@testset "ISOP sensitivity" begin
    u_isop = 0.15008086216045768

    rs1 = structural_simplify(SuperFast())
    o1 = solve(ODEProblem(rs1, [rs1.ISOP => 0.54], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    rs2 = structural_simplify(SuperFast())
    o2 = solve(ODEProblem(rs2, [rs2.ISOP => 0.13], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    test2 = o1[rs1.O3][end] - o2[rs2.O3][end]

    @test test2 ≈ u_isop
end

@testset "NO2 sensitivity" begin
    u_no2 = 35.75701073785234

    rs1 = structural_simplify(SuperFast())
    o1 = solve(ODEProblem(rs1, [rs1.NO2 => 100.0, rs1.DMS => 0.1], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    rs2 = structural_simplify(SuperFast())
    o2 = solve(ODEProblem(rs2, [rs2.DMS => 0.1], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    test3 = o1[rs1.O3][end] - o2[rs2.O3][end]

    @test test3 ≈ u_no2
end

@testset "CO sensitivity" begin
    u_co = -0.1680113594146384

    rs1 = structural_simplify(SuperFast())
    o1 = solve(ODEProblem(rs1, [rs1.CO => 50.0], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    rs2 = structural_simplify(SuperFast())
    o2 = solve(ODEProblem(rs2, [rs2.CO => 500.0], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    test4 = o1[rs1.O3][end] - o2[rs2.O3][end]

    @test test4 ≈ u_co
end

@testset "CH4 sensitivity" begin
    u_ch4 = 0.00463593346606217

    rs1 = structural_simplify(SuperFast())
    o1 = solve(ODEProblem(rs1, [rs1.CH4 => 1900.0], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    rs2 = structural_simplify(SuperFast())
    o2 = solve(ODEProblem(rs2, [rs2.CH4 => 1600.0], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    test5 = o1[rs1.O3][end] - o2[rs1.O3][end]

    @test test5 ≈ u_ch4
end
