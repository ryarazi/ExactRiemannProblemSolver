#struct for saving the state of the system in some half
struct HydroStatus
	rho::Float64 #density
	u::Float64 #velocity
	p::Float64 #pressure
	gamma::Float64 #adiabatic index
	c::Float64 #speed of cound
	A::Float64 #Toro constant (4.8)
	B::Float64 #Toro constant (4.8)
	
	function HydroStatus(rho, u, p, gamma)
		c = sqrt(gamma * p / rho)
		A = 2. / (gamma + 1.) / rho
		B = (gamma - 1.) / (gamma + 1.) * p
		new(rho, u, p, gamma, c, A, B)
	end
end

#like the HydroStatus but saves only the important parameters - density, pressure and velocity
struct HydroStatusSimple
	rho::Float64 #density
	u::Float64 #velocity
	p::Float64 #pressure
end
HydroStatusSimple(status::HydroStatus) = HydroStatusSimple(status.rho, status.u, status.p)

#different types to dispatch the first guess of p_star calculation
#see "PressureGuess.jl"
abstract type FirstGuessScheme end
struct TwoRarefaction <: FirstGuessScheme end
struct PrimitiveValue <: FirstGuessScheme end
struct TwoShock <: FirstGuessScheme end	
struct MeanPressure <: FirstGuessScheme end	

#the type of the wave that will be created on some side of the contact discontinuity
abstract type WaveType end
struct Rarefaction <: WaveType end
struct Shock <: WaveType end
struct Vaccum <: WaveType end #for vaccum generated in the middle

#the side in which the calculation is happening
abstract type Side end
struct LeftSide <: Side end
struct RightSide <: Side end

#in some equations there is factor difference between the left and right side calculations
side_factor(::LeftSide) = -1
side_factor(::RightSide) = 1

#choose the some variable based on the side given
choose_side(left, right, ::LeftSide) = left
choose_side(left, right, ::RightSide) = right

