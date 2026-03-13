using Revise, SynapseElife,
	Random,
	Plots,
	PiecewiseDeterministicMarkovProcesses,
	ColorSchemes, Sundials

############# Initials ##################################################
data_protocol = dataProtocol("TigaretMellor16")
############################  Simulation 1 - panel b ################################
k = 8
pls = 1
start = 0.5e3

#####EVENT REPRESENTATION
events_times, is_pre_or_post_event = firingPattern(
			start_time=start,
			n_pos  		 = data_protocol[!,:n_pos][k],
			delay_pos	 = data_protocol[!,:delay_pos][k],
			n_pre  		 = data_protocol[!,:n_pre][k],
			delay_pre	 = data_protocol[!,:delay_pre][k],
			delay  		 = data_protocol[!,:delay][k],
			pulse  		 = pls, #data_protocol[:pulse][k],		# PULSE
			freq   		 = data_protocol[!,:freq][k],			# FREQUENCY
			causal 		 = data_protocol[!,:causal][k],
			repeat_times = data_protocol[!,:repeat_times][k],	# FREQUENCY
			repeat_after = data_protocol[!,:repeat_after][k])

if data_protocol[!,:post_poisson_rate][k] > 0. #for Poisson, mandatory positive rate, the regular rate in the parameters is used to set the event_time end
	post = Synapse.PoissonProcess(data_protocol[!,:post_poisson_rate][k], start, events_times[end] )
		post_ = repeat([0],inner=length(post))
	pre = Synapse.PoissonProcess(data_protocol[!,:pre_poisson_rate][k], start, events_times[end] )
		pre_ = repeat([0],inner=length(pre))
	events_times=[pre;post]
	p=sortperm(events_times)
	global is_pre_or_post_event=[pre_;post_][p]
	global events_times=[pre;post][p]
end

param_synapse = SynapseParams(
			t_end		   = start+(data_protocol[!,:repeat_times][k]+1)*pls*1000/data_protocol[!,:freq][k]+.7e3,
			soma_dist 		= 200. ,
			temp_rates 		= data_protocol[!,:temp][k],
			Ca_ext			= data_protocol[!,:exca][k],
			Mg				= data_protocol[!,:exmg][k],
			age				= data_protocol[!,:age][k],
			injbap			= data_protocol[!,:inj_time][k],
			I_clamp			= data_protocol[!,:injection][k],
			sampling_rate   = 10.)

p = param_synapse.p_release
pre_synapse = PreSynapseParams(h = 0) # A zero is placed here since we don't want to have a release failure when showing panel B in Figure 1.

xc0 = initial_conditions_continuous_temp(param_synapse) # initial conditions deterministic vars
xd0 = initial_conditions_discrete(param_synapse)		# initial conditions stochastic channels

# RUN PRESYNAPTIC MODEL
is_glu_release, Docked, Reserve, t_stp, glu_release_times, bap_by_epsp_times = stp(param_synapse.t_end, 
						pre_synapse, 
						events_times, 
						is_pre_or_post_event, 
						_plot = false, 
						algo = CHV(CVODE_BDF()))
@show "number of releases $(sum(is_glu_release))"

Random.seed!(7)
result = @time evolveSynapse(
	xc0, xd0,
	param_synapse,
	events_times[events_times .< param_synapse.t_end],
	is_pre_or_post_event,
	ifelse(data_protocol[!,:AP_by_EPSP][k] == "yes", bap_by_epsp_times, Float64[]), #optional BaP induced by EPSP
	is_glu_release,
	(CHV(CVODE_BDF(linear_solver=:GMRES)), CHV(CVODE_BDF(linear_solver=:GMRES)));
	# (CHV(:lsoda), CHV(:lsoda));
	abstol = 1e-6, reltol = 1e-5,
	save_positions = (false, true),
	verbose = true); # model function 0.25s


tt = result.t
XD = result.XD

colorss = ColorSchemes.viridis
limits_d(x) = (minimum(x) ;maximum(x))

gr()
begin	
	l = @layout [[a{.5w} b{.3w}]; [c{.5w}  d{.3w}]; [e{.5w} f{.3w}]; [g{.5w} h{.3w}]]
	Plots.plot(windowsize=(1.5*700 * .6,350*1.5),layout=l,grid=false)

	#AMPAr
	lim=[(start-15*.02),(start+15*.28)]
	val = XD[14,:] .+ XD[15,:] .+ XD[16,:]
	plot!(tt,val,xlabel=" ",label="",ylabel="AMPAr",w=2.5,subplot=1,alpha=1,color=get(colorss, 0.) ,linetype=:steppost,layout=l,axis=false,grid=false,xlim=[start-20 ,start .+ 250],xticks=(collect(500:100:700),collect(0:100:200))
	,yticks=(collect(0:20:80)))
	plot!(tt,val,xlabel="time (ms)",label="",ylabel="AMPAr",w=2.5,subplot=2,alpha=1,color=get(colorss, 0.) ,linetype=:steppost,layout=l,xlim=lim,
	xticks=(collect(500:1:504),collect(0:1:4)),
	yticks=(collect(0:20:80)))

	#VGCC
	lim = [(start+11.2),(start+20)]
	val = XD[28,:]
	plot!(tt,XD[32,:],xlabel="",label="T-type",ylabel="",w=2.5,subplot=3,alpha=1,color=get(colorss, .5) ,linetype=:steppost,layout=l)
	plot!(tt,XD[28,:],xlabel="",label="R-type",ylabel="",w=2.5,subplot=3,alpha=1,color=get(colorss, 1.) ,linetype=:steppost,layout=l,axis=false,grid=false)
	plot!(tt,XD[34,:] .+ XD[35,:],xlabel="",label="L-type",ylabel="VGCC",w=2.5,subplot=3,alpha=1,color=get(colorss, .0) ,linetype=:steppost,layout=l,
	xlim=[start-20 ,start .+ 250],xticks=(collect(500:100:700),collect(0:100:200)))

	plot!(tt,XD[28,:],xlabel="",label="",ylabel="",w=2.5,subplot=4,alpha=1,color=get(colorss, 1.) ,linetype=:steppost,layout=l)
	plot!(tt,XD[32,:],xlabel="",label="",ylabel="",w=2.5,subplot=4,alpha=1,color=get(colorss, .5) ,linetype=:steppost,layout=l)
	plot!(tt,XD[34,:] .+ XD[35,:],xlabel="time (ms)",label="",ylabel="VGCC",w=2.5,subplot=4,alpha=1,color=get(colorss, .0) ,linetype=:steppost,layout=l,xlim=lim,
	xticks=(collect(512:2:520),collect(0:2:8)))

	#GABAr
	lim = [(start-100*.015),(start+100*.3)]
	val = XD[49,:] .+ XD[50,:]
	plot!(tt,val,xlabel="",label="",ylabel="GABAr",w=2.5,subplot=5,alpha=1,color=get(colorss, 0.5) ,linetype=:steppost,layout=l,axis=false,grid=false,
	xlim=[start-20 ,start .+ 250],xticks=(collect(500:100:700),collect(0:100:200)))
	plot!(tt,val,xlabel="time (ms)",label="",ylabel="GABAr",w=2.5,subplot=6,alpha=1,color=get(colorss, 0.5) ,linetype=:steppost,layout=l,
	xlim=lim,
	xticks=(collect(500:10:530),collect(0:10:30)),yticks=(collect(0:10:30)))

	#NMDAr
	lim = [(start-550*.06 ),(start+550 )]
	val = XD[22,:] .+ XD[23,:]
	plot!(tt,val,xlabel="",label="GluN2A",ylabel="",w=2.5,subplot=7,alpha=1,color=get(colorss, 0.25) ,linetype=:steppost,layout=l,yaxis=false,grid=false,xlim=[start-20,start .+ 250])
	plot!(tt,val,xlabel="time (ms)",label="",ylabel="NMDAr",w=2.5,subplot=8,alpha=1,color=get(colorss, 0.25) ,linetype=:steppost,layout=l,xlim=[500-5,600],
	xticks=(collect(500:25:600),collect(0:25:100)),yticks=(collect(0:2:8)))
	val = XD[44,:] .+ XD[45,:]
	plot!(tt,val,xlabel="time (ms)",label="GluN2B",ylabel="NMDAr",w=2.5,subplot=7,alpha=1,color=get(colorss, 0.75) ,linetype=:steppost,layout=l,
	xticks=(collect(500:100:700),collect(0:100:200)),yticks=(collect(0:2:8))) |> display
end
############################  Simulation 1 - panel c and d ################################
begin
	l = @layout [[a{.5w} b{.5w}]; [c{.5w}  d{.5w}]]
	Plots.plot(windowsize=(1.5*700 * .6,350*1),layout=l,grid=false)

	k = 8
	pls = 1
	start = .5e3

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

		param_synapse = SynapseParams(
					t_end		   = start+(data_protocol[!,:repeat_times][k]+1)*pls*1000/data_protocol[!,:freq][k]+2e3,
					soma_dist 		= 200. ,
					temp_rates 		= data_protocol[!,:temp][k],
					Ca_ext			= data_protocol[!,:exca][k],
					Mg				= data_protocol[!,:exmg][k],
					age				= data_protocol[!,:age][k],
					injbap			= data_protocol[!,:inj_time][k],
					I_clamp			= data_protocol[!,:injection][k],
					sampling_rate   = 0.1)

		p = param_synapse.p_release
		pre_synapse = PreSynapseParams(h = 0)


		xc0 = initial_conditions_continuous_temp(param_synapse)
		xd0 = initial_conditions_discrete(param_synapse)

		#####RUN MODEL
		is_glu_release, Docked, Reserve, t_stp, glu_release_times, bap_by_epsp_times = stp(param_synapse.t_end, pre_synapse, events_times, is_pre_or_post_event, algo = CHV(CVODE_BDF()))
		@show "number of releases $(sum(is_glu_release))"


		for i in 1:5
		result = @time evolveSynapse(
				xc0,
				xd0,
				param_synapse,
				events_times[events_times .< param_synapse.t_end],
				is_pre_or_post_event,
				ifelse(data_protocol[!,:AP_by_EPSP][k] == "yes",bap_by_epsp_times,Float64[]), #optional BaP induced by EPSP
				is_glu_release,
				# (CHV(:lsoda), CHV(:lsoda));
				(CHV(CVODE_BDF(linear_solver=:GMRES)), CHV(CVODE_BDF(linear_solver=:GMRES)));
				abstol = 1e-6, reltol = 1e-5,
				save_positions = (false, true),
				verbose = false) # model function


			tt = result.t
			out=SynapseElife.get_names(result.XC, result.XD)

		args = (color = :black, label = "", xlabel="time (ms)", xlim = [470, 570],alpha=.3,w=2,xticks=(collect(500:25:550),collect(0:25:50)))
		plot!(tt, out[:Vsp]; subplot = 1, ylabel = "Voltage (mv)", args...)
		plot!(tt,  out[:Ca];  subplot = 2, ylabel = "Ca2+ (μM)", args...) |> display
	end

	k = 8
	pls = 30
	start = .5e3

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

	param_synapse = SynapseParams(
				t_end		   = start+(data_protocol[!,:repeat_times][k]+1)*pls*1000/data_protocol[!,:freq][k]+2e3,
				soma_dist 		= 200. ,
				temp_rates 		= data_protocol[!,:temp][k],
				Ca_ext			= data_protocol[!,:exca][k],
				Mg				= data_protocol[!,:exmg][k],
				age				= data_protocol[!,:age][k],
				injbap			= data_protocol[!,:inj_time][k],
				I_clamp			= data_protocol[!,:injection][k],
				sampling_rate   = 0.1)

	p = param_synapse.p_release
	pre_synapse = PreSynapseParams(h = (p[4]+ p[3]/(1+exp(p[1]* (data_protocol[!,:exca][k]-p[2])))))

	xc0 = initial_conditions_continuous_temp(param_synapse)
	xd0 = initial_conditions_discrete(param_synapse)

	#####RUN MODEL
	is_glu_release, Docked, Reserve, t_stp, glu_release_times, bap_by_epsp_times = stp(param_synapse.t_end, pre_synapse, events_times, is_pre_or_post_event, _plot = false, algo = CHV(CVODE_BDF()))
	@show "number of releases $(sum(is_glu_release))"

	plot!(t_stp/1000, Docked; label = "docked",subplot = 3,xlabel="time(s)", ylabel = "Vesicles",linetype=:steppost ,w=2)
	plot!(t_stp/1000, Reserve; label = "reserve",subplot = 3,xlabel="time(s)", ylabel = "Vesicles",linetype=:steppost ,w=2)
	scatter!(glu_release_times/1000, 31 .* ones(length(glu_release_times)); label = "releases",subplot = 3,xlabel="time (s)", ylabel = "Vesicles",w=2)

	result = @time evolveSynapse(
			xc0,
			xd0,
			param_synapse,
			events_times[events_times .< param_synapse.t_end],
			is_pre_or_post_event,
			ifelse(data_protocol[!,:AP_by_EPSP][k] == "yes", bap_by_epsp_times,Float64[]), #optional BaP induced by EPSP
			is_glu_release,
			(CHV(CVODE_BDF(linear_solver=:GMRES)), CHV(CVODE_BDF(linear_solver=:GMRES)));
			abstol = 1e-6, reltol = 1e-5,
			save_positions = (false, true),
			verbose = false, progress = true) # model function

	tt = result.t
	out = SynapseElife.get_names(result.XC, result.XD)

	args = (color = :black, label = "", xlabel="time (s)", w=2 )
	plot!(tt/1000, out[:λ]; subplot = 4, ylabel = "BaP efficiency ", args...) |>display
end
############################  Simulation 1 - panel c and d ################################
begin
	l = @layout [a ; b ; c]
	Plots.plot(windowsize=(1*500, 1*400),layout=l,grid=false)

	k = 8
	pls = 30
	start = .5e3

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

		param_synapse = SynapseParams(
					t_end		   = start+(data_protocol[!,:repeat_times][k]+1)*pls*1000/data_protocol[!,:freq][k]+1.5e5,
					soma_dist 		= 200. ,
					temp_rates 		= data_protocol[!,:temp][k],
					Ca_ext			= data_protocol[!,:exca][k],
					Mg				= data_protocol[!,:exmg][k],
					age				= data_protocol[!,:age][k],
					injbap			= data_protocol[!,:inj_time][k],
					I_clamp			= data_protocol[!,:injection][k],
					sampling_rate   = 0.1)

		p = param_synapse.p_release
		pre_synapse = PreSynapseParams(h = (p[4]+ p[3]/(1+exp(p[1]* (data_protocol[!,:exca][k]-p[2])))))


		xc0 = initial_conditions_continuous_temp(param_synapse)
		xd0 = initial_conditions_discrete(param_synapse)

		#####RUN MODEL
		is_glu_release, Docked, Reserve, t_stp, glu_release_times, bap_by_epsp_times = stp(param_synapse.t_end, pre_synapse, events_times, is_pre_or_post_event, _plot = false, algo = CHV(CVODE_BDF()))
		@show "number of releases $(sum(is_glu_release))"

		result = @time evolveSynapse(
				xc0,
				xd0,
				param_synapse,
				events_times[events_times .< param_synapse.t_end],
				is_pre_or_post_event,
				ifelse(data_protocol[!,:AP_by_EPSP][k] == "yes",bap_by_epsp_times,Float64[]), #optional BaP induced by EPSP
				is_glu_release,
				# (CHV(:lsoda), CHV(:lsoda));
				(CHV(CVODE_BDF(linear_solver=:GMRES)), CHV(CVODE_BDF(linear_solver=:GMRES)));
				abstol = 1e-6, reltol = 1e-5,
				save_positions = (false, true),
				verbose = false, progress = true) # model function


		tt = result.t
		out = SynapseElife.get_names(result.XC, result.XD)

		CaMKII  =  out[:KCaM0] .+ out[:KCaM2C] .+ out[:KCaM2N] .+ out[:KCaM4] .+ out[:PCaM0] .+ out[:PCaM2C] .+ out[:PCaM2N] .+ out[:PCaM4] .+ out[:P] .+ out[:P2]
		CaM	 =   out[:CaM2C] .+ out[:CaM2N] .+ out[:CaM4]
		CaN = out[:CaN4]


		args = (color = :black, label = "", xlabel="time (s)" ,alpha=1,w=2)
		plot!(tt/1000, CaM; subplot = 1, ylabel = "CaM (μM)", args...)
		plot!(tt/1000,  CaMKII;  subplot = 2, ylabel = "CaMKII (μM)", args...)
		plot!(tt/1000,  CaN;  subplot = 3, ylabel = "CaN (μM)", args...)  |> display
end