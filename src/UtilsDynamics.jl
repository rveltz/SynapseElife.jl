# nspeedfactor = 1, speed up Kdr channel rates
@inline alpha_n(V, nspeedfactor = 1) = nspeedfactor * exp(-0.11*(V - 13) )
@inline beta_n(V, nspeedfactor = 1)  = nspeedfactor * exp(-0.08*(V - 13) )
@inline alpha_m(V) = 0.4   * (V + 30) / (1 - exp(-(V + 30)/7.2))
@inline beta_m(V)  = 0.124 * (V + 30) / (exp((V + 30)/7.2) - 1)
@inline alpha_h(V) = 0.01  * (V + 45) / (exp((V + 45)/1.5) - 1) #error!! this is beta_h
@inline beta_h(V)  = 0.03  * (V + 45) / (1 - exp(-(V + 45)/1.5) ) #error!! this is alpha_h

##########################################
function mollifier(t::T, duration; pw = 20) where T
	if abs(t/duration) > 10
		return zero(T)
	else
		return one(T) / (one(T) + (t/duration)^pw)
	end
end

"""
$(SIGNATURES)

Compute the current from the BAPs. It is a smooth version of a sum of dirac with fixed amplitude `amp`.
"""
function inputBaP(t, bapTimes::Vector, duration, amp::T) where T
	if isempty(bapTimes); return zero(T); end
	res = zero(T)
	Δ = duration / 2
	for ts in bapTimes
		res += mollifier(t - (ts + Δ), Δ)
	end
	return res * amp
end

##########  ltd-ltp regions (rectangular or elliptical)  ##########
@inline function plasticityRate(p, nhill, K)
	Pmax = 1
	r = p^nhill
	return Pmax * r / (r + K^nhill)
end

@inline function sigmoid(p, nhill, K)
	Pmax = 1
	r = p^nhill
	return Pmax * r / (r + K^nhill)
end

##########  VGCC  ##########
function rates_m_r(Vsp)
	beta_m_r_star = 1/(4e-1) # /ms
	minf_m_r_star = 1/(1+exp((3-10)/8))
	alpha_m_r_star = beta_m_r_star*minf_m_r_star/(1-minf_m_r_star)
	tau_m_r = 1/(alpha_m_r_star + beta_m_r_star)
	minf_r = 1/(1+exp((3-Vsp)/8))
	alpha_m_r = minf_r/tau_m_r
	beta_m_r = (1-minf_r)/tau_m_r
	return alpha_m_r, beta_m_r
end

function rates_h_r(Vsp)
	tau_h_r = 100 # ms
	hinf_r = 1/(1 + exp((Vsp+39)/9.2))
	alpha_h_r = hinf_r/tau_h_r
	beta_h_r = (1-hinf_r)/tau_h_r
	return alpha_h_r, beta_h_r
end

function rates_m_t(Vsp)
	beta_m_t_star = 1 # /ms
	minf_m_t_star = 1/(1+exp((-32+20)/7))
	alpha_m_t_star = beta_m_t_star*minf_m_t_star/(1-minf_m_t_star)
	tau_m_t = 1/(alpha_m_t_star + beta_m_t_star)
	minf_t = 1/(1+exp((-32-Vsp)/7))
	alpha_m_t = minf_t/tau_m_t
	beta_m_t = (1-minf_t)/tau_m_t
	return alpha_m_t, beta_m_t
end

function rates_h_t(Vsp)
	tau_h_t = 50 # ms
	hinf_t = 1/(1 + exp((Vsp+70)/6.5))
	alpha_h_t = hinf_t/tau_h_t
	beta_h_t = (1-hinf_t)/tau_h_t
	return alpha_h_t, beta_h_t
end

@inline function rates_l(Vsp)
	return 0.83/(1+exp((13.7-Vsp)/6.1)),
			0.53/(1+exp((Vsp-11.5)/6.4)),
			1.86/(1+exp((Vsp-18.8)/6.17))
end

@inline CaNnorm(s, mu, x) = 1/(s*sqrt(2*pi))*exp(-0.5*((x-mu)/s)^2)
@inline CaT_inf_m(v) = 1/(1+exp((-32-v)/7))
@inline CaT_inf_h(v) = 1/(1+exp((v-(-67))/6.5))
@inline CaR_inf_m(v) = 1/(1+exp((3-v)/8.3))
@inline CaR_inf_h(v) = 1/(1+exp((v-(-32))/9.2))

@inline function CaL_inf_m(v)
	theta_m = -10
	kappa_m = -6
	return  1 / (1 + (exp((v - theta_m)/ kappa_m)))
end

########### NMDA and Ca flux  ##########
# Term representing NMDA voltage-dependent Mg unblock
@inline B(v, Mg) = 1 / (1 + exp(- 0.062 * v) * (Mg / 3.57 ) ) #Jahr Stevens

@inline function ghk(V, Ca_int, Ca_ext, p_synapse)
	# Kelvin physiological temp.
	@unpack_SynapseParams p_synapse
	x = z * V * faraday / (gas*(temp_rates + 273.15))
	v = z * x * faraday * ((Ca_int)- (Ca_ext) * exp(-x)) / (1 - exp(-x))
	return v
end

########### SK  ##########
@inline SK_chnnl(Ca) = Ca^6/(Ca^6 + 0.333^6) #0.333
