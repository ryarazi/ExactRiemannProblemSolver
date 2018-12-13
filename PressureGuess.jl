include("TypeDefine.jl")

#see (4.46) from Toro
function p_guess(left::HydroStatus, right::HydroStatus, ::TwoRarefaction)
		@assert left.gamma == right.gamma
		gamma = left.gamma
		gamma_power = (gamma-1) / 2. / gamma
		delta_u = right.u - left.u
		return ( (left.c + right.c - 0.5 * (gamma - 1) * delta_u) /
			   (left.c/left.p^gamma_power + right.c/right.p^gamma_power) ) ^
		       (1. / gamma_power)
end

#see (4.47) from Toro
function p_guess(left::HydroStatus, right::HydroStatus, ::PrimitiveValue)
	delta_u = right.u - left.u
	mean_pressure = 0.5 * (left.p + right.p)
	mean_density = 0.5 * (left.rho + right.rho)
	mean_speed_of_sound = 0.5 * (left.c + right.c)
	return mean_pressure - delta_u * mean_density * mean_speed_of_sound
end

#see (4.48) from Toro
g(p, status::HydroStatus) = sqrt(status.A / (p + status.B))	
function p_guess(left::HydroStatus, right::HydroStatus, ::TwoShock)
	p_hat = p_guess(left, right, PrimitiveValue)
	g_left = g(p_hat, left)
	g_right = g(p_hat, right)
	delta_u = right.u - left.u
	return (g_right * left.p + g_left * right.p - delta_u)/(g_right + g_left)
end				  

#see (4.49) from Toro
p_guess(left::HydroStatus, right::HydroStatus, ::MeanPressure) = 0.5*(left.p + right.p)