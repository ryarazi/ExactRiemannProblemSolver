# ExactRiemannProblemSolver

A Julia program which generates exact solutions to Riemann problems for the hydrodynamics Euler equations based on the algorithm of Toro which is described in:
Toro, E.F., Riemann Solvers and Numerical Methods for Fluid Dynamics, A Practical Introduction, 3nd ed., Springer, Berlin, 2009.

This program also knows how to handle case of vaccum (both vaccum in the initial conditions and generation of vaccum through thr dynamics).
To see how to use the program you can read the documentation of the main function sample_riemann inside "RiemannSolver.jl".
In addition to the solver code, the file RiemannTests.jl have set of test problems from different sources to check the code.

## Examples

An example of the pressure, density, energy and velocity profiles printed by the "SodShockTube" test, which solves the famous Sod Shock Tube problem:

![Sod Shock Tube Graphs](sod_example.png?raw=true)

An example of the pressure, density, energy and velocity profiles printed by the some test (which is called "RCVCR" in the tests file) in which vaccum is generated:

![RCVCR](vaccum_example.png?raw=true)

## This fork

This fork makes the program a Julia package, with the solver in `src/`, and extends it to
stiffened gases, `p = (gamma - 1) rho e - gamma p_inf`, with `gamma` and `p_inf` on each
side: `HydroStatus(rho, u, p, gamma, p_inf)`, where `p_inf = 0` (the default) is the ideal
gas. A stiffened gas is an ideal gas in `p + p_inf`, so the changes put `p + p_inf` in place
of `p`; each changed line is marked `stiffened gas`. As before, the iteration for `p_star` is
in `log p`, so `p_star` must be positive.

```julia
using ExactRiemannProblemSolver
left = HydroStatus(rho_l, u_l, p_l, gamma_l, p_inf_l)
right = HydroStatus(rho_r, u_r, p_r, gamma_r, p_inf_r)
profile = sample_riemann(x, t, left, right)    # rho, u and p at the points x at the time t
```
