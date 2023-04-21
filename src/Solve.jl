"""
$(SIGNATURES)


"""

"""
$(SIGNATURES)

Perform a simulation of the synapse model. Among other things, you need to provide the external events impacting the synapse: Glutamate releases and BaPs.

# Arguments
- `xc0` initial condition for the continuous variables. Example `xc0 = Synapse.initial_conditions_continuous_temp(p_synapse)`
- `xd0` initial condition for the discrete variables. Example `xd0 = Synapse.initial_conditions_discrete(p_synapse)`
- `p_synapse::SynapseParams` synapse parameters. Example `p_synapse = SynapseParams()`.
- `events_sorted_times` sorted list of times (ms) for external events (Glutamate / BaP)
- `is_pre_or_post_event::Vector{Bool}` whether the corresponding event in `events_sorted_times` is a Glutamate event. If `false`, it corresponds to a BaP event.
- `bap_by_epsp::Vector{<:Real}` Additionnal BaPs time events, these are evoked by EPSPs. There are added to the ones indexed in `is_events_glu`.
- `is_glu_released::Vector{Bool}` variable to shut down the Glutamate event, i.e. make the glutamate amplitude be zero. This proves useful to have this variable for reproducing some experiments. If equals to `false`, then glutamate amplitude is set to zero.
- `algos` simulation algorithms from `PiecewiseDeterministicMarkovProcesses`. For example `(PDMP.CHV(:lsoda), PDMP.CHV(:lsoda))`

# Optional arguments
- `verbose = false` display information during simulation
- `abstol = 1e-8` absolute tolerance for ODE time stepper
- `reltol = 1e-7` relative tolerance for ODE time stepper
- `progress = false` show a progressbar during simulation
- `save_positions = (false, true)` save the values (before, after) the jumps (transitions)
- `nu` transition matrix. It is initialised with `buildTransitionMatrix()`.
"""
function evolveSynapse(xc0::Vector{T}, xd0, p_synapse::SynapseParams,
	events_sorted_times, is_pre_or_post_event, bap_by_epsp,
	is_glu_released, nu, algos, agg = nothing;
	verbose = false, progress = false, abstol = 1e-8, reltol = 1e-7,
	save_positions = (false, true), kwargs...) where T

	tt, XC, XD = evolveSynapse_noformat(xc0, xd0, p_synapse,
		events_sorted_times, is_pre_or_post_event, bap_by_epsp,
		is_glu_released, nu, algos, agg; verbose = verbose, progress = progress,
		abstol = abstol, reltol = reltol, save_positions = save_positions, kwargs...)

	# format the output to make it convenient to parse
	# this is wasting a lot of ressources but is convenient for plotting
	verbose && @printf("=> done! parsing results")

	out = formatSynapseResult(tt, XC, XD)
end


"""
Same as `evolveSynapse` but do not format the output because it takes time.
"""
function evolveSynapse_noformat(xc0::Vector{T}, xd0, p_synapse::SynapseParams,
	events_sorted_times, is_pre_or_post_event, bap_by_epsp,
	is_glu_released, nu, algos, agg = nothing;
	verbose = false, progress = false, abstol = 1e-8, reltol = 1e-7,
	save_positions = (false, true), kwargs...) where T

	if save_positions isa Tuple{Bool, Bool}
		save_positionsON = save_positions
		save_positionsOFF = save_positions
	else
		save_positionsON = save_positions[1]
		save_positionsOFF = save_positions[2]
	end

	@assert eltype(is_pre_or_post_event) == Bool "Provide booleans for glutamate releases."
	@assert eltype(is_glu_released) == Bool "Provide booleans for glutamate indices."

	verbose && printstyled(color=:red,"\n+++++++++++++++++++++++++++++\n")
	verbose && printstyled(color=:red,"Synapse simulation")

	XC = VectorOfArray([xc0]) # vector to hold continuous variables
	if isnothing(agg)
		XD = VectorOfArray([xd0]) # vector to hold discrete variables
	else
		XD = VectorOfArray([typeof(xc0)(xd0)])
	end
	tt = [0.0] # vector of times

	# we collect which external events correspond to BaPs
	events_bap = events_sorted_times[is_pre_or_post_event .== false]

	# function to simulate the synapse when Glutamate is ON
	SimGluON = (xc, xd, t1, t2, glu) -> SynapseProblem(xc, xd, t1, t2, events_bap, bap_by_epsp, glu, p_synapse, nu, algos[1], agg; save_positions = save_positionsON, reltol = reltol, abstol = abstol, kwargs...)

	# function to simulate the synapse when Glutamate is OFF
	SimGluOFF = (xc, xd, t1, t2) -> SynapseProblem(xc, xd, t1, t2, events_bap, bap_by_epsp, zero(T), p_synapse, nu, algos[2], agg; save_positions = save_positionsOFF, reltol = reltol, abstol = abstol, kwargs...)

	# variable to display progressbar during simulation
	# +1 for the last big till p_synapse.t_end
	pbar = progress ? Progress(length(events_sorted_times) + 1, 1) : nothing

	# random variable for Glutamate concentration
	gluDist = Gamma(1/p_synapse.glu_cv^2, p_synapse.glu_cv^2)

	# we loop over the external events, simulate them and append to res
	for (eveindex, eve) in enumerate(events_sorted_times)
		verbose && printstyled(color=:red,"\n$(eveindex) / $(length(events_sorted_times)) +++++++++++++++++++++++\n")
		if is_pre_or_post_event[eveindex] == true # it is a pre-synaptic event
			# we simulate the synapse with Glutamate OFF until event time
			# then we put  Glutamate ON for dt = p_synapse.glu_width with variable amplitude (concentration)
			verbose && @printf("=> Glu Off,%4d, t ∈ [%9.4e, %9.4e]\n", eveindex, tt[end], eve)

			# simulate the event with Glutamate OFF
			res = SimGluOFF(XC[:,end], XD[:,end], tt[end], eve)
			formatSimResult!(res, XC, XD, tt)
			gluamp = rand(gluDist)
			verbose && @printf("=> Glu on, %4d, t ∈ [%9.4e, %9.4e]\n", eveindex, eve, eve+ p_synapse.glu_width )

			# simulate the event with Glutamate ON
			# variability here
			res = SimGluON(XC[:,end], XD[:,end], eve, eve + p_synapse.glu_width,  ifelse(is_glu_released[eveindex], gluamp, zero(T)))
			formatSimResult!(res, XC, XD, tt)
		end
		# update the progress bar
		progress && next!(pbar; showvalues = [(:steps, eveindex), (:t, tt[end])])
	end

	# reaching tend: we simulate the synapse with Glutamate OFF until simulation end time required
	# by the user. In  most protocol, this is taking most of the time.
	verbose && @printf("=> Reaching the end, t ∈ [%9.4e, %9.4e]\n",tt[end], p_synapse.t_end)
	res = @time SimGluOFF(XC[:,end], XD[:,end], tt[end], p_synapse.t_end)
	formatSimResult!(res, XC, XD, tt)
	if isnothing(agg)
		@info "last bit" length(res.time) tt[end] p_synapse.t_end
	else
		@info "last bit" length(res.t) tt[end] p_synapse.t_end
	end

	# update the progress bar
	progress && next!(pbar; showvalues = [(:steps, length(events_sorted_times) + 1), (:t, p_synapse.t_end)])

	if tt[end] != p_synapse.t_end 
            @warn "The simulation did not reach requested simulated time."
        end

	return (t = tt, XC = XC, XD = XD)
end

function formatSimResult!(res::PDMP.PDMPResult, XC, XD, tt)
	append!(XC, res.xc)
	append!(XD, res.xd)
	append!(tt, res.time)
	nothing
end

function formatSimResult!(res::ODESolution, XC, XD, tt)
	if res.u isa ExtendedJumpArray
		u = VectorOfArray([u.u for u in res.u]) 
	else
		u = VectorOfArray(res.u)
	end
	append!(XC, VectorOfArray([i[1:34] for i in u]))
	append!(XD, VectorOfArray([i[35:84] for i in u]))
	append!(tt, res.t)
	nothing
end

function formatSynapseResult(tt, XC, XD)
	namesC = (:Vsp, :Vdend, :Vsoma, :λ, :ImbufCa, :Ca, :Dye, :CaM0, :CaM2C,
		:CaM2N, :CaM4, :mCaN, :CaN4, :mKCaM, :KCaM0, :KCaM2N, :KCaM2C, :KCaM4,
		:PCaM0, :PCaM2C, :PCaM2N, :PCaM4, :P, :P2, :LTD, :LTP, :act_D, :act_P,
		:m, :h, :n, :SK ,:λ_age, :λ_aux)
	values = (XC[i, :] for i in 1:length(namesC))
	return (t = tt, XD = XD, XC = XC, zip(namesC, values)...)
end

function indexOfVariable(name::Symbol)
	names = (:Vsp, :Vdend, :Vsoma, :λ, :ImbufCa, :Ca, :Dye, :CaM0, :CaM2C,
		:CaM2N, :CaM4, :mCaN, :CaN4, :mKCaM, :KCaM0, :KCaM2N, :KCaM2C, :KCaM4,
		:PCaM0, :PCaM2C, :PCaM2N, :PCaM4, :P, :P2, :LTD, :LTP, :act_D, :act_P,
		:m, :h, :n, :SK ,:λ_age, :λ_aux)
	return findfirst(isequal(name), names)
end

function getCaM(t, XC, XD)
	XC[indexOfVariable(:CaM2C), :] .+
	XC[indexOfVariable(:CaM2N), :] .+
	XC[indexOfVariable(:CaM4), :]
end

getCaN(t, XC, XD) = XC[indexOfVariable(:CaN4), :]

function getCamKII(t, XC, XD)
	return XC[indexOfVariable(:KCaM0), :]   .+
		XC[indexOfVariable(:KCaM2C), :] .+
		XC[indexOfVariable(:KCaM2N), :] .+
		XC[indexOfVariable(:KCaM4), :]  .+
		XC[indexOfVariable(:PCaM0), :]  .+
		XC[indexOfVariable(:PCaM2C), :] .+
		XC[indexOfVariable(:PCaM2N), :] .+
		XC[indexOfVariable(:PCaM4), :]  .+
		XC[indexOfVariable(:P), :]      .+
		XC[indexOfVariable(:P2), :]
end
