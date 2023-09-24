using Revise, SynapseElife,
	Random,
	Plots,
	PiecewiseDeterministicMarkovProcesses,
	ColorSchemes, Parameters, Sundials

############# Initials ##################################################


function makeFigure2(;detailed::Bool = true)
	# extract protocols
	data_protocol = dataProtocol("TigaretMellor16")

	# make plot layout
	l = @layout [[a{.5w} b{.5w}]; [c{.5w}  d{.5w}]; [e{.5w} f{.5w}]]
	Plots.plot(windowsize = (1.5*700 * .6,350*1.5), layout = l, grid = false)

	colorss = ColorSchemes.viridis

	counter = 1
	K = 0
	for s in [[7, 8], [6, 4], [3, 2]]
		counter = counter + K
		colorrr = rand()
		K = 0
		for k in s
			@info data_protocol[k,:protocol]
			K = K + 1
			pls = data_protocol[!,:pulse][k]
			start = .5e3
			#####EVENT REPRESENTATION
			events_times, is_pre_or_post_event = firingPattern(
						start_time = start,
						n_pos  		 = data_protocol[!,:n_pos][k],
						delay_pos	 = data_protocol[!,:delay_pos][k],
						n_pre  		 = data_protocol[!,:n_pre][k],
						delay_pre	 = data_protocol[!,:delay_pre][k],
						delay  		 = data_protocol[!,:delay][k],
						pulse  		 = pls,#data_protocol[:pulse][k],		# PULSE
						freq   		 = data_protocol[!,:freq][k],			# FREQUENCY
						causal 		 = data_protocol[!,:causal][k],
						repeat_times = data_protocol[!,:repeat_times][k],	# FREQUENCY
						repeat_after = data_protocol[!,:repeat_after][k])
			afterend = 1.5e5
			param_synapse = SynapseParams(
						t_end		   = start+(data_protocol[!,:repeat_times][k]+1)*pls*1000/data_protocol[!,:freq][k]+afterend,
						soma_dist 		= 200. ,
						temp_rates 		= data_protocol[!,:temp][k],
						Ca_ext			= data_protocol[!,:exca][k],
						Mg				= data_protocol[!,:exmg][k],
						age				= data_protocol[!,:age][k],
						injbap			= data_protocol[!,:inj_time][k],
						I_clamp			= data_protocol[!,:injection][k],
						sampling_rate   = .01)

			p = param_synapse.p_release
			pre_synapse = PreSynapseParams(h = (p[4]+ p[3]/(1 + exp(p[1] * (data_protocol[!,:exca][k] - p[2])))))

			xc0 = initial_conditions_continuous_temp(param_synapse)
			xd0 = initial_conditions_discrete(param_synapse)

			#####RUN MODEL
			is_glu_release, Docked, Reserve, t_stp, glu_release_times, bap_by_epsp_times = stp(param_synapse.t_end, pre_synapse, events_times, is_pre_or_post_event, _plot = false, algo = CHV(CVODE_BDF()))
			@show "number of releases $(sum(is_glu_release))"

			# this is slightly optimised to limit allocations and
			result = @time evolveSynapse_noformat(
					xc0, xd0,
					param_synapse,
					events_times[events_times .< param_synapse.t_end],
					is_pre_or_post_event[1:length(events_times)],
					data_protocol[!,:AP_by_EPSP][k] == "yes" ? bap_by_epsp_times : Float64[], #optional BaP induced by EPSP
					is_glu_release,
					# (CHV(:lsoda), CHV(:lsoda));
					(CHV(CVODE_BDF(linear_solver=:GMRES)), CHV(CVODE_BDF(linear_solver=:GMRES)));
					abstol = 1e-6, reltol = 1e-5,
					save_positions = ((false, false), (false, detailed)),
					verbose = false, progress = true) # model function

			@info "Extracting data..." length(result.t)
			tt = result.t
			out = SynapseElife.get_names(result.XC, result.XD)

			CaMKII = getCamKII(result...)
			CaM = getCaM(result...)
			CaN = getCaN(result...)

			@info "Plotting..."
			args = (color = get(colorss, colorrr), w = 2)
			plot!(tt/1000, CaMKII; subplot = counter,xlabel="time (s)", ylabel = "enzymes (μM)",ylim=[0,32], label = "$(data_protocol[k,:protocol])", args...)
			plot!(tt/1000,  CaN;  subplot = counter,label = "",linestyle = :dash,args...)
			plot!(CaN,  CaMKII;  subplot = counter+1,ylabel="CaMKII (μM)", xlabel="CaN (μM)", label = "$(data_protocol[k,:protocol])",ylim=[0,32],xlim=[0,11.5], args...) |> display
			colorrr=1-colorrr #invert color
		end
	end
end
@time makeFigure2(detailed = true)
