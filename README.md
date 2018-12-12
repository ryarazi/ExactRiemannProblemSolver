# ExactRiemannProblemSolver
A Julia program which generates exact solutions to Riemann problems for the hydrodynamics Euler equations based on the algorithm of Toro which is described in:
Toro, E.F., Riemann Solvers and Numerical Methods for Fluid Dynamics, A Practical Introduction, 3nd ed., Springer, Berlin, 2009.

This program also knows how to handle case of vaccum (both vaccum in the initial conditions and generation of vaccum through thr dynamics).
To see how to use the program you can read the documentation of the main function sample_riemann inside "RiemannSolver.jl".
In addition to the solver code, the file RiemannTests.jl have set of test problems from different sources to check the code.