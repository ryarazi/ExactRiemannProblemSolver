using Roots
include("TypeDefine.jl")
include("PressureGuess.jl")

#function (4.6) and (4.7) from Toro
f(p, status::HydroStatus, ::Shock) = (p - status.p) * sqrt(status.A / (p + status.B))
f(p, status::HydroStatus, ::Rarefaction) = 2. * status.c / (status.gamma - 1.) * ((p / status.p)^((status.gamma - 1.) / 2. /status.gamma) - 1.)
function f(p, status::HydroStatus)
	wave_type = p > status.p ? Shock() : Rarefaction()
	return f(p, status, wave_type)
end

#solver for p_star based on the first guess of the pressure(see PressureGuess.jl)
#uses secant method to find solve the equation
function p_star_calc(left::HydroStatus, right::HydroStatus, guess_scheme::T, TOL::Float64) where {T<:FirstGuessScheme}	
	p0 = max(TOL, p_guess(left, right, guess_scheme))
	
	#we define p = exp(u) so inside the interation it won't get negative value
	u0 = log(p0) #p = exp(u) => u = log(p)
	f_star(u) = f(exp(u), left) + f(exp(u), right) + right.u - left.u #function (4.5) from Toro
	u = find_zero(f_star, u0, Order1(), rtol = TOL) #Secant method
	return exp(u) #p = exp(u) 
end

#equation (4.9) from Toro
u_star_calc(p_star, left::HydroStatus, right::HydroStatus) = 0.5 * (left.u + right.u) +
															 0.5 * (f(p_star, right) - f(p_star, left))


#equations (4.23) and (4.32) from Toro
rho_star_calc(p_star, status::HydroStatus, ::Rarefaction) = status.rho * (p_star / status.p) ^ (1. / status.gamma)

#equation (4.25) and (4.34) from Toro
c_star_calc(p_star, status::HydroStatus, ::Rarefaction) = status.c * (p_star / status.p) ^ ((status.gamma - 1.) / 2. / status.gamma)

#equations (4.50) and (4.57) from Toro
function rho_star_calc(p_star, status::HydroStatus, ::Shock)
	gamma_ratio = (status.gamma - 1.) / (status.gamma + 1.)
	pressure_ratio = p_star / status.p
	return status.rho * (pressure_ratio + gamma_ratio) / (pressure_ratio * gamma_ratio + 1.)
end

#equation (4.55) and (4.62) from Toro 
head_speed_calc(p_star, u_star, status::HydroStatus,
				side::T, ::Rarefaction) where {T<:Side} = status.u + side_factor(side) * status.c

#equation (4.55) and (4.62) from Toro 
function tail_speed_calc(p_star, u_star, status::HydroStatus, side::T, ::Rarefaction) where {T<:Side}
	c_star = c_star_calc(p_star, status, Rarefaction())
	return u_star + side_factor(side) * c_star
end

#equations (4.52) and (4.59) from Toro
function head_speed_calc(p_star, u_star, status::HydroStatus, side::T, ::Shock) where {T<:Side}
	gamma = status.gamma
	gamma_ratio1 = (gamma + 1.) / 2. / gamma
	gamma_ratio2 = (gamma - 1.) / 2. / gamma
	return status.u + side_factor(side) * status.c * sqrt(gamma_ratio1 * p_star / status.p + gamma_ratio2)
end

#for shock the head and the tail are the same
tail_speed_calc(p_star, u_star, status::HydroStatus, side::T, ::Shock) where {T<:Side} = head_speed_calc(p_star, u_star, status, side, Shock())


#for vaccum generation the density in the middle is zero (= vaccum)
rho_star_calc(p_star, status::HydroStatus, ::Vaccum) = 0.

#those are the velocities of the head and the tail of the rarefaciton wave when a vaccum is generated
#see (4.76), (4.77), (4.79) and (4.80) from Toro
head_speed_calc(p_star, u_star, status::HydroStatus,
				side::T, ::Vaccum) where {T<:Side} = status.u + side_factor(side) * status.c
tail_speed_calc(p_star, u_star, status::HydroStatus, side::T, ::Vaccum) where {T<:Side} = status.u - side_factor(side) * 2 * status.c / (status.gamma - 1.)

#this is the density, velocity and pressure profile of a rarefaction wave in some specific x/t	
#see (4.56) and (4.63) from Toro
function rarefaction_profile(x_over_t, status::HydroStatus, side::T) where {T<:Side}
	gamma = status.gamma
	gamma_ratio = (gamma - 1.) / (gamma + 1.)
	gamma_plus = 2. / (gamma + 1.)
	gamma_minus = 2. / (gamma - 1.)
	
	rarefaction_factor = (gamma_plus - side_factor(side) * gamma_ratio / status.c * (status.u - x_over_t)) ^ gamma_minus
	
	rho = status.rho * rarefaction_factor
	u = gamma_plus * (-side_factor(side)*status.c + (1. /gamma_minus) * status.u + x_over_t)
	p = status.p * rarefaction_factor ^ gamma
	
	return HydroStatusSimple(rho, u, p)
end

#sample the riemann problem in the situation where both sides are not vaccum (but vaccum can be generated through the dynamics)
function sample_riemann_regular(x, t::Float64, left::HydroStatus, right::HydroStatus,
								  guess_scheme::T, TOL::Float64) where {T<:FirstGuessScheme}
	@assert left.rho != 0.
	@assert right.rho != 0.
	
	profile = similar(x, HydroStatusSimple) #we return the values of rho, u and p in each point x in space
	
	#the status "far" from the center = like the starting condition
	status_left = HydroStatusSimple(left)
	status_right = HydroStatusSimple(right)
	
	x_over_t = x ./ t
	
	if 2. / (left.gamma - 1.) * (left.c + right.c) <= right.u - left.u #vaccum generation
		#if vaccum is generated - the density in the middle will be zero
		p_star = 0. 
		u_star = 0.
		wave_type_left = Vaccum()
		wave_type_right = Vaccum()
	else
		#if vaccum is not generated we can calculate the star profile in the middle
		p_star = p_star_calc(left, right, guess_scheme, TOL)
		u_star = u_star_calc(p_star, left, right)
		#check what kind of waves are in every direction (shock or rarefaction)
		wave_type_left = p_star > left.p ? Shock() : Rarefaction()
		wave_type_right = p_star > right.p ? Shock() : Rarefaction()
	end

	rho_star_left = rho_star_calc(p_star, left, wave_type_left)
	status_left_star = HydroStatusSimple(rho_star_left, u_star, p_star) #the profile in the left near the contact discontinuity
	head_speed_left = head_speed_calc(p_star, u_star, left, LeftSide(), wave_type_left)
	tail_speed_left = tail_speed_calc(p_star, u_star, left, LeftSide(), wave_type_left)
	
	rho_star_right = rho_star_calc(p_star, right, wave_type_right)
	status_right_star = HydroStatusSimple(rho_star_right, u_star, p_star) #the profile in the right near the contact discontinuity
	head_speed_right = head_speed_calc(p_star, u_star, right, RightSide(), wave_type_right)
	tail_speed_right = tail_speed_calc(p_star, u_star, right, RightSide(), wave_type_right)
	
	
	for i = 1:length(x)
		S = x_over_t[i] #this is like the S which Toro use in Section 4.5
		
		#see Figure 4.14 in Toro for the flow of the following lines
		side = S < u_star ? LeftSide() : RightSide()
		status = choose_side(left, right, side)
		head_speed = choose_side(head_speed_left, head_speed_right, side)
		tail_speed = choose_side(tail_speed_left, tail_speed_right, side)
		
		#this is used to flip the direction of the inequallity in the branching between the left and the right side
		#the xor which will be in the following lines uses this boolean like some "controlled not"
		#when right_condition is "true" the inequallity will be flipped
		#when right_condition is "false" the inequallity will be stay the same
		right_condition = isa(side, RightSide)

		if xor(S < head_speed, right_condition)
			profile[i] = choose_side(status_left, status_right, side)
		elseif xor(S < tail_speed, right_condition) #can only happen in Rarefaction, because in Shock  head_speed == tail_speed
			profile[i] = rarefaction_profile(S, status, side)
		else
			profile[i] = choose_side(status_left_star, status_right_star, side)
		end
	end
	
	return profile							   
end

#sample the riemann problem in the situation where the gas in "side" is actually a vaccum
function sample_riemann_side_vaccum(x, t::Float64, status::HydroStatus, side::T) where {T<:Side}
	
	profile = similar(x, HydroStatusSimple) #we return the values of rho, u and p in each point x in space
	vaccum = HydroStatusSimple(0., 0., 0.) #vaccum
	simple = HydroStatusSimple(status)
	
	x_over_t = x ./ t
	
	tail_speed = tail_speed_calc(0., 0., status, side, Vaccum()) #the head of the rarefaciton wave
	head_speed = head_speed_calc(0., 0., status, side, Vaccum()) #the tail of the rarefaciton wave (=the boundary with the vaccum)
	
	#see "sample_riemann_regular" for explnation about this boolean
	right_condition = isa(side, RightSide)
	
	for i = 1:length(x)
		S = x_over_t[i]
		
		if xor(S < head_speed, right_condition)
			profile[i] = simple
		elseif xor(S < tail_speed, right_condition)
			profile[i] = rarefaction_profile(S, status, side)
		else
			profile[i] = vaccum
		end
	end
	
	return profile							   
end

#this is the main function that samples the Riemann problem and returns the density, pressure and velocity in the space points
#specified by the array x. It gets the following parameters:
# x - an array type which saves all the points in space in which we want to sample the problem (pay attention - the contact discontinuity begins at x=0.)
# t - a floating number which tells us in what time we sample to problem
# left + right - struct from type HydroStatus (see "TypeDefine.jl") to get the initial state of the system
# guess_scheme - the type of guess which will be used to solve the p_star equation (see "PressureGuess.jl")
# TOL - the tolerance for the convergance of the p_star finding algorithm
function sample_riemann(x, t::Float64, left::HydroStatus =  HydroStatus(1.,    0., 1.,  7. / 5.),
									   right::HydroStatus = HydroStatus(0.125, 0., 0.1, 7. / 5.),
									   guess_scheme::T = PrimitiveValue(), TOL::Float64 = 1.e-6) where {T<:FirstGuessScheme}
	@assert left.gamma == right.gamma
	
	if left.rho == 0 #material is on the right, vaccum on the left
		return sample_riemann_side_vaccum(x, t, right, RightSide())
	elseif right.rho == 0 #material is on the left, vaccum on the right
		return sample_riemann_side_vaccum(x, t, left, LeftSide())
	else
		return sample_riemann_regular(x, t, left, right, guess_scheme, TOL)
	end

end
