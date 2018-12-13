using PyPlot
include("RiemannSolver.jl")

N = 1000
xmin = -0.5
xmax = 0.5
x = collect(LinRange(xmin, xmax, N))
save_figs = true #if this is false the figures will use show(), otherwise they will be saved locally in 'figs' dir

#create directory for figures
if save_figs && !isdir("figs")
	mkdir("figs")
end

struct TestCase
	name::String
	t::Float64
	left::HydroStatus
	right::HydroStatus
end

function calculate_and_plot(test::TestCase)
	name = test.name
	
	println("Time for the test $name:")
	
	#dummy run to compile all code
	sample_riemann(x, test.t, test.left, test.right)
	
	@time profiles = sample_riemann(x, test.t, test.left, test.right)
	density = [status.rho for status in profiles]
	velocity = [status.u for status in profiles]
	pressure = [status.p for status in profiles]
	energy = pressure ./ density ./ (test.left.gamma - 1.)
	
	figure()
	suptitle(name)
	subplot(221)
	plot(x, density)
	title("density")
	xlabel("x")
	grid(true)
	subplot(222)
	plot(x, velocity)
	title("velocity")
	xlabel("x")
	grid(true)
	subplot(223)
	plot(x, pressure)
	title("pressure")
	xlabel("x")
	grid(true)
	subplot(224)
	plot(x, energy)
	title("energy")
	xlabel("x")
	grid(true)
	tight_layout()
	
	if !save_figs
		show()
	else
		
		savefig("figs/$name.pdf",format="pdf",dpi=200)
		savefig("figs/$name.png",format="png",dpi=200)
	end
end

gamma = 7. / 5.

#The first 5 case are taken from Toro's book: (Table 4.1)
#Toro, E.F., Riemann Solvers and Numerical Methods for Fluid Dynamics, A Practical Introduction, 3nd ed., Springer, Berlin, 2009.
#cases 6-13 are takeb from:
#Kamm et. al. 'Enhanced Verification Test Suite for Physics Simulation Codes' SAND2008-7813 2009

#case 1: 
sod = TestCase("SodShockTube",  0.25, HydroStatus(1., 0., 1., gamma),
									  HydroStatus(0.125, 0., 0.1, gamma))
calculate_and_plot(sod)

#case 2: 
einfeldt = TestCase("Einfeldt123",  0.15, HydroStatus(1., -2., 0.4, gamma),
									      HydroStatus(1., 2., 0.4, gamma))
calculate_and_plot(einfeldt)

#case 3:
toro3 =  TestCase("ToroTest3",  0.012, HydroStatus(1., 0., 1000., gamma),
									   HydroStatus(1., 0., 0.1, gamma))
calculate_and_plot(toro3)

#case 4: 
toro4 =  TestCase("ToroTest4",  0.035, HydroStatus(1., 0., 0.01, gamma),
									   HydroStatus(1., 0., 100., gamma))
calculate_and_plot(toro4)

#case 5: 
toro5 =  TestCase("ToroTest5",  0.035, HydroStatus(5.99924, 19.5975, 460.894, gamma),
									   HydroStatus(5.99242, -6.19633, 46.0950, gamma))
calculate_and_plot(toro5)

#case 6: 
vaccum_left =  TestCase("VaccumExpansionLeft",  0.75, HydroStatus(0., 0., 0., gamma),
													  HydroStatus(1., 0., 1., gamma))
calculate_and_plot(vaccum_left)

#case 7: 
vaccum_right =  TestCase("VaccumExpansionRight",  0.75, HydroStatus(1., 0., 1., gamma),
														HydroStatus(0., 0., 0., gamma))
calculate_and_plot(vaccum_right)

#case 8: 
rcvcr =  TestCase("RCVCR",  0.75, HydroStatus(1., -4., 0.4, gamma),
								  HydroStatus(1., 4., 0.4, gamma))
calculate_and_plot(rcvcr)

#case 9:
modified_sod = TestCase("ModifiedSod",  0.2, HydroStatus(1., 0.75, 1., gamma),
											 HydroStatus(0.125, 0., 0.1, gamma))
calculate_and_plot(modified_sod)

#case 9:
stream = TestCase("StreamCollision",  0.8, HydroStatus(1., 2., 0.1, gamma),
										   HydroStatus(1., -2., 0.1, gamma))
calculate_and_plot(stream)

#case 10:
leblanc = TestCase("LeBlanc",  0.5, HydroStatus(1., 0., (2. / 3.)*1.e-1, gamma),
									HydroStatus(1.e-3, 0., (2. / 3.)*1.e-10, gamma))
calculate_and_plot(leblanc)

#case 11:
peak = TestCase("PeakProblem",  3.9e-3, HydroStatus(0.1261192, 8.9047029, 782.92899, gamma),
									    HydroStatus(6.591493, 2.2654207, 3.1544874, gamma))
calculate_and_plot(peak)

#case 12:
slow = TestCase("SlowShock",  2., HydroStatus(3.857143, -0.810631, 10.33333, gamma),
								  HydroStatus(1.0, -3.44, 1.0, gamma))
calculate_and_plot(slow)

#case 13:
stationary = TestCase("StationaryContact",  0.012, HydroStatus(1.0, -19.59745, 1.e3, gamma),
												   HydroStatus(1.0, -19.59745, 1.e-2, gamma))
calculate_and_plot(stationary)

##if you want to profile the tests - remove the following comments
##you will also need a larger N than used before
#
#using Profile
#N = 1000000
#x = collect(LinRange(xmin, xmax, N))
#chosen_test = rcvcr
#sample_riemann(x, test.t, test.left, test.right)
#Profile.clear()
#Profile.clear_malloc_data()
#@time sample_riemann(x, test.t, test.left, test.right)
#@profile sample_riemann(x, test.t, test.left, test.right)
#Profile.print()