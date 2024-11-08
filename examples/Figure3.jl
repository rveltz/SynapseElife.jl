using Revise, SynapseElife,
	Random,
	Plots,
	PiecewiseDeterministicMarkovProcesses,
	ColorSchemes, Sundials

############# Initials ##################################################
data_protocol = dataProtocol("TigaretMellor16")

colorss = ColorSchemes.coolwarm
# pyplot()
l = @layout [a{.5w} b{.5w}]
Plots.plot(windowsize=(0.9*1100*2/3,0.9*250),layout=l,grid=false)
plot!(w=0,SynapseParams().LTP_region,color="red",alpha=.1)
plot!(w=0,SynapseParams().LTD_region,color="blue",alpha=.1); plot!(grid=false)

for k in 2:8
	pls = data_protocol[!,:pulse][k]
	start = 3e3
		#####EVENT REPRESENTATION
		events_times, is_pre_or_post_event = firingPattern(
					start_time=start,
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

		afterend=1.5e5

		param_synapse = SynapseParams(
					t_end		   = start+(data_protocol[!,:repeat_times][k]+1)*pls*1000/data_protocol[!,:freq][k]+afterend,
					soma_dist 		= 200. ,
					temp_rates 		= data_protocol[!,:temp][k],
					Ca_ext			= data_protocol[!,:exca][k],
					Mg				= data_protocol[!,:exmg][k],
					age				= data_protocol[!,:age][k],
					injbap			= data_protocol[!,:inj_time][k],
					I_clamp			= data_protocol[!,:injection][k],
					sampling_rate   = .1)

		p = param_synapse.p_release
		pre_synapse = PreSynapseParams(h = (p[4]+ p[3]/(1+exp(p[1]* (data_protocol[!,:exca][k]-p[2])))))

		xc0 = initial_conditions_continuous_temp(param_synapse)
		xd0 = initial_conditions_discrete(param_synapse)

		#####RUN MODEL
		is_glu_release, Docked, Reserve, t_stp, glu_release_times, bap_by_epsp_times = stp(param_synapse.t_end, pre_synapse, events_times, is_pre_or_post_event, _plot = false, algo = CHV(CVODE_BDF()))

		@info data_protocol[k,:protocol]
		@show "number of releases $(sum(is_glu_release))"

		ode = :lsoda
		ode = CVODE_BDF(linear_solver=:GMRES)

		result = @time evolveSynapse(
				xc0,
				xd0,
				param_synapse,
				events_times[events_times .< param_synapse.t_end],
				is_pre_or_post_event,
				ifelse(data_protocol[!,:AP_by_EPSP][k] == "yes",bap_by_epsp_times,Float64[]), #optional BaP induced by EPSP
				is_glu_release,
				(CHV(ode), CHV(ode), CHV(ode));
				save_positions = (false, true),
				verbose = false, progress = true) # model function

		@info "Extracting data..."
		tt = result.t
		out = SynapseElife.get_names(result.XC, result.XD)

		CaMKII = out[:KCaM0] .+ out[:KCaM2C] .+ out[:KCaM2N] .+ out[:KCaM4] .+ out[:PCaM0] .+ out[:PCaM2C] .+ out[:PCaM2N] .+ out[:PCaM4] .+ out[:P] .+ out[:P2]
		CaM = out[:CaM2C] .+ out[:CaM2N] .+ out[:CaM4]
		CaN = out[:CaN4]

		@info "Plotting..."
		args = (color = get(colorss, k/8), w = 2, grid= false)
		plot!(CaN,  CaMKII; subplot = 1, ylabel="CaMKII (μM)", xlabel="CaN (μM)", label = "", ylim=[0,32], xlim=[0,11.5], args...) |> display
		plot!(tt/1e3, [out[:lt]][1:end][1][1:end,3]-[out[:lt]][1:end][1][1:end,2]; subplot = 2,ylabel="weight change (%)", xlabel="Time (s)", label = "$(data_protocol[k,:protocol])", ylim = [-80,90], args...) |> display
end

plot!(subplot = 1,legend = :none) |> display
