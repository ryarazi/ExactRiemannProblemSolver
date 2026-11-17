module ExactRiemannProblemSolver

using Roots

include("TypeDefine.jl")
include("PressureGuess.jl")
include("RiemannSolver.jl")

export HydroStatus, HydroStatusSimple, sample_riemann, p_star_calc, u_star_calc
export TwoRarefaction, PrimitiveValue, TwoShock, MeanPressure

end
