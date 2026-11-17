using ExactRiemannProblemSolver
using Test

@testset "ExactRiemannProblemSolver" begin
    # Sod's problem, Test 1 of Toro (2009, Chapter 4), whose star region he gives as
    # p* = 0.30313 and u* = 0.92745
    left = HydroStatus(1.0, 0.0, 1.0, 1.4)
    right = HydroStatus(0.125, 0.0, 0.1, 1.4)
    p_star = p_star_calc(left, right, PrimitiveValue(), 1.0e-12)
    @test isapprox(p_star, 0.30313; atol = 5.0e-6)
    @test isapprox(u_star_calc(p_star, left, right), 0.92745; atol = 5.0e-6)

    # A stiffened gas is an ideal gas in p + p_inf: with the same gamma and p_inf on both
    # sides, its solution is that of the ideal gas at p + p_inf, shifted back by p_inf
    p_inf = 3.0
    x = collect(range(-0.5, 0.5; length = 101))
    stiffened = sample_riemann(x, 0.2, HydroStatus(1.0, 0.0, 1.0, 1.4, p_inf),
                               HydroStatus(0.125, 0.0, 0.1, 1.4, p_inf), PrimitiveValue(),
                               1.0e-12)
    ideal = sample_riemann(x, 0.2, HydroStatus(1.0, 0.0, 1.0 + p_inf, 1.4),
                           HydroStatus(0.125, 0.0, 0.1 + p_inf, 1.4), PrimitiveValue(),
                           1.0e-12)
    @test all(isapprox(a.rho, b.rho) && isapprox(a.u, b.u) && isapprox(a.p + p_inf, b.p)
              for (a, b) in zip(stiffened, ideal))
end
