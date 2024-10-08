using Main.GasChem
using EarthSciMLBase
using DifferentialEquations, ModelingToolkit, DynamicQuantities
using ModelingToolkit:t

@testset "2wayCoupling" begin
    sol_middle = 9.9476146376358

    sf = couple(SuperFast(), FastJX())
    combined_mtk = convert(ODESystem, sf)
    sys = structural_simplify(combined_mtk)
    tspan = (0.0, 3600 * 24)

    sol = solve(ODEProblem(sys, [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    @test sol[sys.SuperFast.O3][4320] ≈ sol_middle
end
