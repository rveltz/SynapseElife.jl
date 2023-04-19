"""
$(SIGNATURES)

Calcium-dependent sigmoid for the vesicle release test.

# Arguments
- `Ca_pre` phenomenological presynaptic calcium
- `nhill` steepness of the sigmoid
- `K` h-activation value
"""

@inline function releaseProbaSTP(Ca_pre, nhill, K)
	return Ca_pre^nhill / (Ca_pre^nhill + K^nhill)
end

"""
$(SIGNATURES)

Deterministic part of the vesicle release. The `Ca_pre` (phenomenological presynaptic calcium)
has exponential decay constant. The variable `Ca_jump` is the adaptation of the deterministic jump of `Ca_pre`.
`V_evoke` has exponential decay it is used to determine whether an EPSP will trigger an AP.

# Arguments
- `xdot` derivative vector
- `pop_c` continuos variables
- `discrete_var` discrete variables, unnused here
- `p_pre::PreSynapseParams` list of parameters for the presynaptic side
- `t` time
"""
function stp_F_synapse!(xdot, pop_c, discrete_var, p_pre::PreSynapseParams, t)
	Ca_pre, Ca_jump, V_evoke  = pop_c
	@unpack_PreSynapseParams p_pre

	# Calcium
	xdot[1]   =  -Ca_pre / τ_pre
	# Decay term
	xdot[2]   =  (1-Ca_jump) / τ_rec - δ_ca * Ca_jump * Ca_pre
	# AP induced EPSP
	xdot[3]   =  -V_evoke / τ_V
	xdot
end

"""
$(SIGNATURES)

Discrete part of the vesicle release. Estimate the transitions rate of transition between the
different pools respecting their limit capacity.

- `τ_R` is the mixing constant from docked (xd[2]) to the recovery pool (xd[3]).
- `τ_D` is the recycling constant from recovery (xd[3]) to the docked pool (xd[2]).
- `τ_R_ref` is the recycling constant of recovery pool (xd[3]), representing a reffiling
- `R_0` is initial capacity of the recovery pool (xd[3])
- `D_0` is initial capacity of the docked pool (xd[3])

# Arguments
- `rate` rate of each transition
- `xc` continuous variables, unnused here
- `xd` discrete variables
- `p_pre::PreSynapseParams` list of parameters for the presynaptic side
- `t` time
- `issum::Bool` variable used for the rejection algorithm
"""
function stp_R_synapse!(rate, xc, xd, p_pre::PreSynapseParams, t, issum::Bool)
	# the discrete variables are [nPrint, nDocked, nReserve]
	# unpack the parameters
	@unpack_PreSynapseParams p_pre

	rate[1]  = sampling_rate			   # sampling rate for plotting
	rate[2]  = (R_0 - xd[3]) * xd[2] / τ_R
	rate[3]  = (D_0 - xd[2]) * xd[3] / τ_D
	rate[4]  = (R_0 - xd[3]) / τ_R_ref

	if issum == false
		return 0., 0. # we do not use rejection algorithm
	else
		return sum(rate), 0.
	end
end

"""
$(SIGNATURES)
Definition of the PDMP problem for the presynaptic side.

# Arguments
- `xc0, xd0` initial conditions for continuous (3d) and discrete variables (3d)
- `parms` list of parameters for the presynaptic side
- `nustp` stp_build_transition_matrix()` transition matrix
- `ti` intial time of the segment
- `tf` final time of the segment

# Keyword arguments
`algo = PDMP.CHV(:lsoda)` LSODA solver, used to solve this PDMP
"""
function stpPDMP(xc0, xd0, parms, nu_stp, ti, tf; algo = PDMP.CHV(:lsoda), kwargs...)
	problem = PDMP.PDMPProblem(stp_F_synapse!, stp_R_synapse!, nu_stp, xc0, xd0, parms, (ti, tf))
	return solve(problem, algo; kwargs...)
end

"""
$(SIGNATURES)

Evolution rule for presynaptic side PDMP. It first solves the PDMP before the presynaptic stimulation, then
a deterministic jump is applied to `Ca_pre` (xd[1]) and `V_evoke`, the variable used to estimate the AP induced by EPSP (xd[3]).
Then a test using `Ca_pre` (xd[2]) is done to evaluate if a vesicle will be released (given the availability in the docked pool), once it is released,
the docked pool (xd[2]) loses one vesicle. If the test is successful the variable `event`, which marks the presynaptic stimulation being tested,
will be added to the `release time`, the list saving the successful glutamate releases. Similarly `V_evoke` (xd[3]) is used to evaluate if
an AP was triggered by an EPSP. However, the successful "AP induced by EPSP" are saved 15 ms after the EPSP originating it, and saved in the list
 `release_time_auxbap`. After the cycle of tests, the simulation continues until the end or until the next presynaptic stimulation.

# Arguments
- `t_end` end simulation time
- `xc0, xd0` initial conditions for continuous (3d) and discrete variables (3d)
- `par_pre::PreSynapseParams` parameters for the presynaptic side
- `prespike::Vector{Float64}` times at which the presynaptic stimulation occur which will determine released glutamate times
- `nu` transition matrix for the discrete reactions of the presynaptic side

# Output
- `tt` time vector of the simulation
- `XC` continuous variables for every time in the simulation
- `XD` discrete variables for every time in the simulation
- `release_time` list of successful presynaptic stimulation
- `release_time_auxbap` list of the AP induced by EPSPs
"""
function stp_evolve_synapse(t_end,
			xc0, xd0,
			par_pre::PreSynapseParams,
			prespike::Vector{Float64},
			nu = stp_build_transition_matrix();
			kwargs...)
	@assert findfirst(prespike .== 1) == nothing "Remove BaP from the list!!"

	XC = VectorOfArray([xc0])
	XD = VectorOfArray([xd0])
	tt = [0.0]
	res = PDMP.PDMPResult([0.,0.], copy(XC), copy(XD))

	release_time = Vector{Float64}(undef, 0)
	release_time_auxbap = Vector{Float64}(undef, 0)

	for (countloop, event) in enumerate(prespike)
		################### Glutamate  ###################
		res = stpPDMP(res.xc[:,end], res.xd[:,end], par_pre, nu, tt[end], event; kwargs...)

		append!(XC, res.xc); append!(XD, res.xd); append!(tt, res.time)

		# we perform a deterministic jump
		res.xc[1,end] += res.xc[2,end] # Ca_pre -> Ca_pre + Ca_jump
		res.xc[3,end] += 1. 		   # V_evoke -> V_evoke + 1

		if rand() < releaseProbaSTP(res.xc[1,end], par_pre.s, par_pre.h)
			# we may have a  Glutamate release
			if res.xd[2, end] > 0
				# we have a docked vesicule
				# we have a Glutamate release
				res.xd[2, end] -= 1
				# we save the Glu release time
				push!(release_time, event)
			end
		end

		# this is for Bap
		if res.xd[2, end] > 0
			auxbap = sum(rand(par_pre.D_0) .< releaseProbaSTP(res.xc[3,end], par_pre.s, par_pre.h))
			if auxbap > .8 * par_pre.D_0 					# 80%
				push!(release_time_auxbap, event + par_pre.δ_delay_AP)
			end
		end

		res = stpPDMP(res.xc[:,end],res.xd[:,end], par_pre, nu, event, event + .1; kwargs...)
		append!(XC, res.xc); append!(XD, res.xd); append!(tt, res.time)
	end
	res = stpPDMP(res.xc[:,end],res.xd[:,end], par_pre, nu, tt[end], t_end; kwargs...)
	append!(XC, res.xc); append!(XD, res.xd); append!(tt, res.time)
	@assert res.time[end] == t_end "Error in PDMP. Did not reach requested simulated time"
	return tt, XC, XD, release_time, release_time_auxbap
end

"""
$(SIGNATURES)

This function performs the simulation of the presynaptic side.

# Arguments:
- `t_end` end of simulation time
- `param` named tuple with parameters, example: `(τ_rec = 20000, δ_decay	= .0004, tau_pre = 20, tau_soma = 40)`
- `all_events_times` sorted list of times of external events
- `is_pre_or_post_index` whether the external events are pre (`true`) or post (`false`)

# Keyword arguments
- `_plot = false` whether to plot the result or not
- `nu_stp = stp_build_transition_matrix()` transition matrix

# Output
- `is_glu_release` `true` (success) or `false` (failure), list of presynaptic stimuli in `all_events_index` that led to a vesicle (glutamate) release
- `D`, also (XD[2,:]), the evolution of the docked pool
- `R`, also (XD[3,:]), the evolution of the reserve pool
- `time` simulation times
- `glu_release_times` the vector of times in which a successful releases (glutamate) occurred
- `bap_by_epsp_times` the vector of times in which an AP was triggered by an EPSP
"""
function stp(t_end, param,
			all_events_times,
			is_pre_or_post_index;
			_plot = false,
			nu_stp = stp_build_transition_matrix(),
			kwargs...)
	# presynaptic spikes
	_prespike = all_events_times[is_pre_or_post_index .== true ]
	prespike = _prespike[_prespike.<t_end]

	tt, XC, XD, glu_release_times, bap_by_epsp_times = stp_evolve_synapse(
				t_end,
				[0., 1. , 0.],
				[0, param.D_0, param.R_0],
				param,
				prespike,
				nu_stp;
				kwargs...) # model function

	is_glu_release = zeros(Bool, length(all_events_times))
	idx = findall(x -> x ∈ glu_release_times, all_events_times)
	is_glu_release[idx] .= true

	return is_glu_release, XD[2,:], XD[3,:], tt, glu_release_times, bap_by_epsp_times, XC[1,:]
end
