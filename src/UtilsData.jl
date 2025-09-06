"""
Returns the number of jumps for each component and `t ∈ stpan`. Usefull to determine who jumps the most. The syntax is as follows

$(SIGNATURES)

where `t` is time, `xd` the discrete component and `It` is an inferior bound on `t`. You can use it like

	statistics_jumps(t, out[:ampa])

"""
function statistics_jumps(t, xd; tspan = (0., Inf64))
	n = size(xd)
	@assert n[1] < n[2] "Perhaps you want to pass xd'?"
	# sort xd by the number of jumps
	Ind = findall((t .>= tspan[1]) .* (t .<= tspan[2]))
	jps = div.(sum(abs, diff(xd[:, Ind], dims = 2), dims = 2), 2) |> vec
	return 1:n[1], jps
end

"""
Returns a dictionary `out` from a simulation with associated names like `out[:ampa]`.
Use it like

$(SIGNATURES)

Another use is for plotting.
"""
function get_names(xc, xd)
	out = Dict{Symbol, Union{Array{Float64,1}, Adjoint{Int64,Array{Int64,2}}}}()
	out[:ampa]	   = xd[1:16,:]'
	out[:nmda]	   = xd[17:23,:]'
	out[:sampling] = xd[24:24,:]'
	out[:vgcc_r]   = xd[25:28,:]'
	out[:vgcc_t]   = xd[29:32,:]'
	out[:vgcc_l]   = xd[33:35,:]'
	out[:lt]	   = xd[36:38,:]'

	ModelNamesContinuous = [:Vsp,:Vdend,:Vsoma,:λ,:ImbufCa,:Ca,:Dye,:CaM0,:CaM2C,:CaM2N,:CaM4,:mCaN,:CaN4,:mKCaM,:KCaM0,:KCaM2N,:KCaM2C,:KCaM4,:PCaM0,:PCaM2C,:PCaM2N,:PCaM4,:P,:P2,:LTD,:LTP,:LTD_act,:LTP_act,:m,:h,:n,:SK,:λ_age,:λ_aux]
	for (key,val) in enumerate(ModelNamesContinuous)
		out[val] = xc[key,:]
	end
	return out
end
####################################################################################################

_protocols = ["oconnor", "Bittner", "Goldings01", "Buchenan", "Tigaret_jitter_timespent", "Tigaret_jitter_double_1", "Poisson", "Tigaret_jitter", "Dudek_jitter", "TigaretMellor16", "TigaretMellor16_poisson", "Sleep_age", "Poisson_physiological_range", "YannisDebanne20_freq", "YannisDebanne20_freq_delay", "YannisDebanne20_ratio", "YannisDebanne20_inv", "Tigaret_burst", "Tigaret_burst_temp", "Tigaret_burst_age", "Tigaret_burst_freq", "Tigaret_burst_ca", "Tigaret_burst_dist", "Tigaret_freq_1200", "Tigaret_tripost", "Tigaret_preburst", "Tigaret_single", "TigaretMellor_sparse", "TigaretMellor_jitter_sparse", "DudekBear92-BCM-Ca", "DudekBear92-BCM", "DudekBear92-BCM-900", "DudekBear92-BCM-37", "DudekBear92-BCM-33", "DudekBear92-BCM-priming", "DudekBear92_timespend", "DudekBear92-100", "DudekBear_short", "FujiBito", "DudekBear92-sliding", "DudekBear92_temp", "DudekBear92_Ca", "DudekBear92_Mg", "DudekBear92_dist", "DudekBear92-Age", "Blocking_age_control", "Blocking_yNMDA", "Blocking_oNMDA", "Blocking_yBaP", "Blocking_oBaP", "Blocking_yGABA", "Blocking_oGABA", "DudekBear92_BCM_recovery", "DudekBear93-LFS", "DudekBear93-TBS", "Cao-TBS", "RecoverLTD", "Chang19", "Fujii_CaN", "YannisDebanne20", "YannisDebanne_temp", "YannisDebanne_age", "Meredith03-GABA", "TigaretvsMeredith", "YannisvsMeredith", "WittenbergWang06_D", "WittenbergWang06_B", "WittenbergWang06_P", "Mizuno01-LTP-Mg", "Mizuno01LTP", "Mizuno01LTD"]

"""
$(SIGNATURES)

Structure to describe a plasticity protocol "conoc"

# Arguments
- `paper` can be `$_protocols`
- `n_pre` number of presynaptic pulses
- `delay_pre` delay between presynaptic pulses
- `n_pos` number of postsynaptic pulses
- `delay_pos` delay between postsynaptic pulses
- `delay` delay between pre and postsynaptic spikes (used in STDP)
- `pulse` number of pulses/pairing repetitions
- `freq` frequency
- `causal` causal inverts the order of pre and post, uses true or false
- `protocol` name of the protocol
- `weight` info entry with weight change value (not used in the model)
- `outcome` `String` to show the qualitative outcome (not used in the model)
- `paper` paper is a `String` to choose a predefined protocol (e.g. paper = TigaretMellor16)
- `temp` temperature by which all temperature-sensitive mechanisms will be adapted
- `injection` (nA) current injected in the soma, used to make postsynaptic spikes (BaPs)
- `exca` (μM) extracellular calcium, we expressed it in μM to be used in the ghk function
- `exmg` (mM) extracellular magnesium, we expressed it in mM, however, it is converted to μM in the magnesium relief function
- `repeat_times` number of epochs additional epochs to include, for instance, TBS is usually expressed in epochs
- `repeat_after` time difference between the epochs (ms)
- `testing_freq` some protocols use test frequencies, this can be useful to evaluate the effect or absence of it, or to add a regular background presynaptic stimuli. 0 Hz will turn it off
- `inj_time` duration (ms) of the current injection to elicit a BaP
- `age` animal age (rat)
- `AP_by_EPSP` `yes` or `no`, to let additional BaPS induced by EPSPs to be included in the stimulation
- `GABA_block` `yes` or `no` GABAr conductance is set zero
- `jitter	` standard deviation of a Gaussian used to jitter the spikes
- `sparse` percentage of spikes skipped
- `dista` distance from the soma in μm
"""
function dataProtocol(paper)
	 	data_protocol = DataFrame(
		n_pre         = Int64[],
		delay_pre     = Float64[],
		n_pos         = Int64[],
		delay_pos     = Float64[],
		delay         = Float64[],
		pulse         = Int64[],
		freq          = Float64[],
		causal        = Bool[],
		protocol      = String[],
		weight        = Float64[],
		outcome       = String[],
		paper         = String[],
		temp          = Float64[],
		injection     = Float64[],
		exca          = Float64[],
		exmg          = Float64[],
		repeat_times  = Int64[],
		repeat_after  = Float64[],
		testing_freq  = Float64[],
		inj_time      = Float64[],
		age           = Float64[],
		AP_by_EPSP    = String[],
		GABA_block	  = String[],
		jitter		  = Float64[],
		sparse	 	  = Float64[],
		dista		  = Float64[],
		pre_poisson_rate = Float64[],
		post_poisson_rate = Float64[])

		########## not shown in the paper
		if paper == "oconnor"
			for freq in [ .5 ; collect(1. : 2 : 9.); collect(10. : 10. : 100.)  ]
				for temp in collect(25. : 1. : 34)
					push!(data_protocol,[1 0. 0 0. 0. 100 freq true "oconnor_$(freq)_$(temp)" -0.5 "LTD" "$paper" temp 1e3 2.2e3 1.2 3 3e5 0. 2. 21 "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## not shown in the paper
		if paper == "Bittner"
			for delay in  collect(50:100:3000)
				push!(data_protocol,[10 50.  12 40.    delay   5 1000/15e3	false    "Bittner_-$(delay)" 0.   "NC"  "$paper" 35. 2e3  2e3 .7 0 0. 0. 3.  50. "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[10 50.  12 40.    delay   5 1000/15e3	true     "Bittner_$(delay)" 0.    "NC"  "$paper" 35. 2e3  2e3 .7 0 0. 0. 3.  50. "yes" "no" 0. 0. 200. 0. 0.])
			end
		end

		########## used in validation (methods)
		if paper == "Goldings01"
				push!(data_protocol,[0 0. 1 0. 0.	 300 5. false "1Post" 0. "ND" "$paper" 35. 1.5e3 2.5e3 1.3 0 0. 0. 2. 60. "no" "yes" 0. 0. 200. 0. 0.])
		end

		########## used in validation (methods)
		if paper == "Buchenan"
				push!(data_protocol,[0 0.  1 0.    0.	  5 100.    false	 "1Post"        0. 	 "ND"  "$paper" 25. 2.5e3 2.5e3 1.3 0 0. 0. 2.  21. "no" "no"  0. 0. 200. 0. 0.])
				push!(data_protocol,[0 0.  1 0.    0.	  5 100.    false	 "1Post"        0. 	 "ND"  "$paper" 25. 2.5e3 2.5e3 1.3 0 0. 0. 2.  55. "no" "no"  0. 0. 200. 0. 0.])
		end

		########## Figure 7
		if paper == "Tigaret_jitter_timespent"
			for s in  [50]
				push!(data_protocol,[1 0.  2 10    50.0   300 5.	false    "2Post-1Pre50_$(s)" 0.   "NC"  "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" s 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    10.    300 5.    true     "1Pre1Post10_$(s)" 0.1  "NC"  "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" s 0. 200. 0. 0.])
			end
		end

		########## Figure 7
		if paper == "Tigaret_jitter_double_1"
			for s in  collect(0:2.5:50)
				push!(data_protocol,[1 0.  2 10    50.    300 5.	false    "2Post-1Pre50_$(s)" 0.   "NC"  "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" s 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    50.	  300 5.    true     "1Pre2Post50_$(s)" 0.5  "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" s 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    10.	  300 5.    true     "1Pre2Post10_$(s)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" s 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    10.    300 5.    true     "1Pre1Post10_$(s)" 0.1  "NC"  "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" s 0. 200. 0. 0.])
			end
		end

		########## Figure 7
		if paper == "Poisson"
			for pre_poisson in [ collect(1. : .5 : 3.) ; collect(4 : 1. : 9.); collect(10. : 3. : 20.);  collect(25. : 5. : 50.) ]
				for post_poisson in [ collect(1. : .5 : 3.) ; collect(4 : 1. : 9.); collect(10. : 3. : 20.);  collect(25. : 5. : 50.) ]
					push!(data_protocol,[1 0. 0  0.  0.   10 1.    true     "Poisson_$(pre_poisson)_$(post_poisson)" 0.1  "NC"  "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2.  35. "yes" "no" 0. 0. 200. pre_poisson post_poisson])
				end
			end
		end


		########## Figure 7
		if paper == "Tigaret_jitter"
			for s in  [0. ; 50. ; 500. ;1000.]
				push!(data_protocol,[2 10.  0 0.   0.	  300 5.    true 	 "2Pre10_$(s)"      0.05 "NC"  "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" s 0. 200. 0. 0.])
				push!(data_protocol,[2 50. 0 0.    0.     900 3.    true     "2Pre50_$(s)"      -0.5 "LTD" "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" s 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    50.    300 5.	false    "2Post1Pre50_$(s)" 0.   "NC"  "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" s 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    20.    300 5.	false    "2Post-1Pre20_$(s)" 0.25 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" s 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    50.	  300 5.    true     "1Pre2Post50_$(s)" 0.5  "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" s 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    10.	  300 5.    true     "1Pre2Post10_$(s)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" s 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    10.    300 5.    true     "1Pre1Post10_$(s)" 0.1  "NC"  "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" s 0. 200. 0. 0.])
			end
		end

		########## not shown in the paper
		if paper == "Dudek_jitter"
			for freq in [ collect(1. : .5 : 3.) ; collect(4 : 1. : 9.); collect(10. : 3. : 20.);  collect(25. : 5. : 50.) ]
				for s in  [collect(2.5: :2.5 : 7.5) ; collect(10 : 10. : 50) ;  collect(50 : 50. : 300);  collect(400 : 100. : 900)]
					push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "BCM_$(freq)_$s"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. 35. "yes" "no" s 0. 200. 0. 0.])
				end
			end
		end

		########## Used in Figure 1, 2 and 3
		if paper == "TigaretMellor16"
				push!(data_protocol,[1 0.  0 0.    0.	  300 5.    true 	 "1Pre"         0.05 "NC"  "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[2 10.  0 0.   0.	  300 5.    true 	 "2Pre10"      0.05 "NC"  "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[2 50. 0 0.    0.     900 3.    true     "2Pre50"      -0.5 "LTD" "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    50.0   300 5.	false    "2Post1Pre50" 0.   "NC"  "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    20.0   300 5.	false    "2Post1Pre20" 0.25 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    50.	  300 5.    true     "1Pre2Post50" 0.5  "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    10.	  300 5.    true     "1Pre2Post10" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    10.    300 5.    true     "1Pre1Post10" 0.1  "NC"  "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
		end

		########## not shown in the paper
		if paper == "TigaretMellor16_poisson"
			for delay in [collect(10:10:200);collect(225:25:500)]
				push!(data_protocol,[1 0.  2 10    delay	300 5.	false	 "2Post1Pre_$delay" 0.   "NC"  "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "yes" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    delay	300 5.	true	 "1Pre2Post_$delay" 0.5  "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "yes" "yes" 0. 0. 200. 0. 0.])
			end
		end

		########## not shown in the paper
		if paper == "Sleep_age"
			for age in collect(5:2.5:45)
				for T in collect(32.0 : .5 :37.5)
					push!(data_protocol,[1 0.  1 0.    0.	  100   20.  true 	 "a_sleep_$(age)_$(T)"         0.05 "NC"  "$paper" T 1e3  1.5e3 1.2 12 1e4 0. 2.  age "yes" "no" 0. 0. 200. 0. 0.])
					push!(data_protocol,[1 0.  1 0.    0.	  100   20.  true 	 "b_sleep_$(age)_$(T)"         0.05 "NC"  "$paper" T 1e3  1.5e3 1.2 12 1e3 0. 2.  age "yes" "no" 0. 0. 200. 0. 0.])
					push!(data_protocol,[1 0.  1 0.    0.	  100   20.  true 	 "c_sleep_$(age)_$(T)"         0.05 "NC"  "$paper" T 1e3  1.5e3 1.2 12 250 0. 2.  age "yes" "no" 0. 0. 200. 0. 0.])

				end
			end
		end

		########## Supp files
		if paper == "Poisson_physiological_range"
				push!(data_protocol,[1 0.  1 0.    10.    900 5.    true     "Poisson_34_5"  0.1    "NC"  "$paper" 34.  1e3  2.5e3 1.3 0 0. 0. 2.  5. "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    10.    900 5.    true     "Poisson_34.5_10" 0.1  "NC"  "$paper" 34.5 1e3  2.5e3 1.3 0 0. 0. 2. 10. "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    10.    900 5.    true     "Poisson_35_15" 0.1  "NC"  "$paper"   34.5 1e3  2.5e3 1.3 0 0. 0. 2. 15. "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    10.    900 5.    true     "Poisson_35_20" 0.1    "NC"  "$paper" 35.  1e3  2.5e3 1.3 0 0. 0. 2. 20. "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    10.    900 5.    true     "Poisson_35_35" 0.1    "NC"  "$paper" 35.  1e3  2.5e3 1.3 0 0. 0. 2. 35. "yes" "no" 0. 0. 200. 0. 0.])
		end

		########## Figure 6
		if paper == "YannisDebanne20_freq"
			for Ca in collect(1.1 : .1 : 3.1)
				for freq in [collect(.2 : .1 : 9.); collect(1. : 1 : 10.)]
					  push!(data_protocol,[1 0.  1 0. abs(10)  100  freq	true   string("$(paper)_10_$(freq)_$Ca") 0. "LTP" "$(paper)" 30.45 .5e3 Ca*1e3 Ca/1.5 0 0. 0. 6. 21 "no" "no" 0.  0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		if paper == "YannisDebanne20_freq_delay"
			for freq in collect(.1: .1 : 1.)
				for delay in [2.5;collect(10.:10.:100.)]
					  push!(data_protocol,[1 0.  1 0. abs(delay)  100  freq	true   string("$(paper)_$(delay)_$(freq)") 0. "LTP" "$(paper)"   30.45 .5e3 2.5*1e3 2.5/1.5 0 0. 0. 6. 21 "no" "no" 0. 0. 200. 0. 0.])
					  push!(data_protocol,[1 0.  1 0. abs(delay)  150  freq	false  string("$(paper)_-$(delay)_$(freq)") 0. "LTD" "$(paper)"  30.   .5e3 2.5*1e3 2.5/1.5 0 0. 0. 6. 14 "no" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		if paper == "YannisDebanne20_ratio"
			for ratio in collect(.5 : .1 : 2.5)
				for Ca in collect(1. : .25 : 3.5)
					  push!(data_protocol,[1 0.  1 0. abs(10)  100  .3	true   string("$(paper)_10_$(Ca)_$(ratio)") 0. "LTP" "$(paper)" 30.5 .5e3 Ca*1e3 Ca/ratio 0 0. 0. 6. 21 "no" "no" 0.  0. 200. 0. 0.])
					  push!(data_protocol,[1 0.  1 0. abs(10)  150  .3	false  string("$(paper)_-10_$(Ca)_$(ratio)") 0. "LTD" "$(paper)" 30.  .5e3 Ca*1e3 Ca/ratio 0 0. 0. 6. 14 "no" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Figure 6
		if paper == "YannisDebanne20_inv"
			for temp in collect(28 : .25 : 32)
				for age in collect(14 : 1. : 21)
					  push!(data_protocol,[1 0.  1 0. abs(10)  100  .3	true   string("$(paper)_$(temp)_$(age)") 0. "LTP" "$(paper)" temp .5e3 1.8*1e3 1.8/1.5 0 0. 0. 6. age "yes" "no" 0.  0. 200. 0. 0.])
				end
			end
		end

		########## Figure 3
		if paper == "Tigaret_burst"
			for dt in collect(10. : 10. : 100.)
			  for pulse in collect(100. : 100. : 1200.)
				push!(data_protocol,[1 0.  2 10    dt	  pulse 5.    true    "burst_$(dt)_$(pulse)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    dt	  pulse 5.    false   "burst_-$(dt)_$(pulse)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
			  end
			end
		end

		########## Supp files
		if paper == "Tigaret_burst_temp"
			for dt in collect(10. : 10. : 100.)
				for temp in collect(30:.25:37)
				  	push!(data_protocol,[1 0.  2 10    dt	  300 5.    true    "burst_$(dt)_$(temp)" 0.75 "LTP" "$paper" temp 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
					push!(data_protocol,[1 0.  2 10    dt	  300 5.    false   "burst_-$(dt)_$(temp)" 0.75 "LTP" "$paper" temp 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		if paper == "Tigaret_burst_age"
			ages = collect(2.5:2.5:80)
			for dt in collect(10. : 10. : 100.)
				for age in ages
				  	push!(data_protocol,[1 0.  2 10    dt	  300 5.    true    "burst_age_$(dt)_$(age)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2. age  "no" "yes" 0. 0. 200. 0. 0.])
					push!(data_protocol,[1 0.  2 10    dt	  300 5.    false   "burst_age_-$(dt)_$(age)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2. age  "no" "yes" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		if paper == "Tigaret_burst_freq"
			for dt in collect(10. : 10. : 100.)
				for freq in [.5; .75 ; collect(1. : 1. : 9.); collect(10. :  5. : 50.) ]
					push!(data_protocol,[1 0.  2 10    dt	  300 freq   true    "burst_freq_$(dt)_$(freq)"  0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
					push!(data_protocol,[1 0.  2 10    dt	  300 freq   false   "burst_freq_-$(dt)_$(freq)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		if paper == "Tigaret_burst_ca"
			for dt in collect(10. : 10. : 100.)
				for Ca in collect(1.1 : .1 : 3.1)
					push!(data_protocol,[1 0.  2  10.  dt 300 5.  true    "burst_ca_$(dt)_$(Ca)"  0.75 "LTP" "$paper" 35. 2e3  Ca*1e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
					push!(data_protocol,[1 0.  2  10.  dt 300 5.  false   "burst_ca_-$(dt)_$(Ca)" 0.75 "LTP" "$paper" 35. 2e3  Ca*1e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		if paper == "Tigaret_burst_dist"
			for dt in collect(10. : 10. : 100.)
				for dista in collect(50. : 50. : 400.)
					push!(data_protocol,[1 0.  2  10.  dt 300 5.  true    "burst_dista_$(dt)_$(dista)"  0.75 "LTP" "$paper" 35. 2e3  2.5*1e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. dista 0. 0.])
					push!(data_protocol,[1 0.  2  10.  dt 300 5.  false   "burst_dista_-$(dt)_$(dista)" 0.75 "LTP" "$paper" 35. 2e3  2.5*1e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. dista 0. 0.])
				end
			end
		end

		########## not shown in the paper
		if paper == "Tigaret_freq_1200"
			for freq in [collect(1. : 1 : 9.); collect(10. : 5. : 50.) ]
				push!(data_protocol,[1 0.  2 10    50	  1200 freq    true    "burstfreq_50_$(freq)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    50	  1200 freq    false   "burstfreq_-50_$(freq)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
			end
		end


		########## not shown in the paper
		if paper == "Tigaret_tripost"
			for dt in collect(10. : 10. : 100.)
			  for pulse in collect(100. : 100. : 900.)
				push!(data_protocol,[1 0.  3 10    dt	  pulse 5.    true    "tripost_$(dt)_$(pulse)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  3 10    dt	  pulse 5.    false   "tripost_-$(dt)_$(pulse)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
			  end
			end
		end

		########## not shown in the paper
		if paper == "Tigaret_preburst"
			for dt in collect(10. : 10. : 100.)
			  for pulse in collect(100. : 100. : 900.)
				push!(data_protocol,[2 10.  1 0    dt	  pulse 5.    true    "preburst_$(dt)_$(pulse)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0.  0. 200. 0. 0.])
				push!(data_protocol,[2 10.  1 0    dt	  pulse 5.    false   "preburst_-$(dt)_$(pulse)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0.  0. 200. 0. 0.])
			  end
			end
		end

		########## not shown in the paper
		if paper == "Tigaret_single"
			for dt in collect(10. : 10. : 100.)
			  for pulse in collect(100. : 100. : 900.)
				push!(data_protocol,[1 0.  1 0.    dt	  pulse 5.    true    "single_$(dt)_$(pulse)" 0.75 "LTP" "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    dt	  pulse 5.    false   "single_-$(dt)_$(pulse)" 0.75 "LTP" "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. 0. 200. 0. 0.])
			  end
			end
		end

		########## Figure 7
		if paper == "TigaretMellor_sparse"
			for i in [collect(0. : .05 : .5);.53;.58;collect(0.65 : .05 : .8)]
				push!(data_protocol,[1 0.  2 10    50.0   300 5.	false    "2Post1Pre50_$i" 0.   "NC"  "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. i 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    50.	  300 5.    true     "1Pre2Post50_$i" 0.5  "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. i 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    10.	  300 5.    true     "1Pre2Post10_$i" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. i 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    10.    300 5.    true     "1Pre1Post10_$i" 0.1  "NC"  "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" 0. i 200. 0. 0.])
			end
		end

		if paper == "TigaretMellor_jitter_sparse"
				for i in collect(0. : 2 : 50)
					for j in collect(0.2 : .1 : 1.)
						push!(data_protocol,[1 0.  2 10    50.0   300 5.	false    "2Post1Pre50_$(i)_$(j)" 0.   "NC"  "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2.  60. "no" "yes" i j 200. 0. 0.])
					end
				end
		end


		########## Supp files
		ages = collect(5:5:60)
		if paper == "DudekBear92-BCM-Ca"
			for freq in [1. ; 3. ; 5. ; 10 ; 30; 50; 70 ; 100.]
				for Ca in collect(1. : .25 : 3.0)
						push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "BCM_$(freq)_$(Ca)"  -0.5 "LTD" "$paper" 35. 2e3 Ca*1e3 1.5 0 0. 0. 3. 35 "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Figure 4
		if paper == "DudekBear92-BCM"
			for freq in [.5; .75 ; collect(1. : .25 : 6.) ; collect(6.5 : .5 : 9.); collect(10. :  5. : 50.) ]
				for pulses in collect(50: 50: 1200)
					for age in [35.]
						push!(data_protocol,[1 0. 0 0.    0.     pulses  freq    true     "BCM_$(age)_$(freq)_$(pulses)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
					end
				end
			end
		end

		########## Figure 4
		if paper == "DudekBear92-BCM-900"
			for freq in [.5; .75 ; collect(1. : .25 : 6.) ; collect(6.5 : .5 : 9.); collect(10. :  5. : 50.) ]
				for pulses in [900]
					for age in [35.]
						push!(data_protocol,[1 0. 0 0.    0.     pulses  freq    true     "BCM_$(age)_$(freq)_$(pulses)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
					end
				end
			end
		end


		########## Supp files
		if paper == "DudekBear92-BCM-37"
			for freq in [.5; .75 ; collect(1. : .25 : 6.) ; collect(6.5 : .5 : 9.); collect(10. :  5. : 50.) ]
				for pulses in collect(50: 50: 1200)
					for age in [35.]
						push!(data_protocol,[1 0. 0 0.    0.     pulses  freq    true     "BCM37_$(age)_$(freq)_$(pulses)"  -0.5 "LTD" "$paper" 37. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
					end
				end
			end
		end

		########## Supp files
		if paper == "DudekBear92-BCM-33"
			for freq in [.5; .75 ; collect(1. : .25 : 6.) ; collect(6.5 : .5 : 9.); collect(10. :  5. : 50.) ]
				for pulses in collect(50: 50: 1200)
					for age in [35.]
						push!(data_protocol,[1 0. 0 0.    0.     pulses  freq    true     "BCM33_$(age)_$(freq)_$(pulses)"  -0.5 "LTD" "$paper" 33. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
					end
				end
			end
		end

		########## not shown in the paper
		if paper == "DudekBear92-BCM-priming"
			for freq in [.5; .75 ; collect(1. : 1 : 9); collect(10. :  5. : 50.) ]
				for priming in [.5; .75 ; collect(1. : 1 : 9.); collect(10. :  5. : 50.) ]
						push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "BCM_priming_$(freq)_$(priming)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. priming 2. 35. "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Figure 4
		if paper == "DudekBear92_timespend"
			for freq in [1.;3.;5;10.;50.]
						push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "BCMtime_$(freq)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. 35. "yes" "no" 0. 0. 200. 0. 0.])
			end
		end

		########## Figure 4
		if paper == "DudekBear92-100"
			for freq in [3. ; 30.]
					push!(data_protocol,[1 0. 0 0.    0.  900  freq    true     "BCM_$(freq)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. 35. "yes" "no" 0. 0. 200. 0. 0.])
			end
		end

		########## not shown in the paper
		if paper == "DudekBear_short"
			for freq in [.5; .75 ; collect(1. : .25 : 6.) ; collect(6.5 : .5 : 9.); collect(10. : 2. : 20.);  collect(25. : 5. : 50.) ]
				for pulses in collect(1: 1: 30)
					for age in [35.]
						# missing one arg
						push!(data_protocol,[1 0. 0 0.    0.     pulses  freq    true     "BCM_$(age)_$(freq)_$(pulses)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
					end
				end
			end
		end

		########## not shown in the paper
		if paper == "FujiBito"
			for freq in [collect(.1 : .2 : 1.) ; collect(1. : 3 : 9.); collect(10. : 5. : 20.) ]
				for pulses in collect(1: 3: 30)
						push!(data_protocol,[1 0.  0 0.    0.     pulses  freq    true     "FujiBito_$(freq)_$(pulses)"  -0.5 "LTD" "$paper"  25. 2e3  2e3  .0001 0 0. 0. 2.  13. "no" "yes" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## not shown in the paper
		if paper == "DudekBear92-sliding"
			for freq in [.8 ; 1.0 ; 1.5 ; collect(2. : 1. : 9.) ; collect(10. : 5. : 50.) ; collect(60. : 10 : 100.) ]
				for pulses in collect(50: 50: 1200)
					for age in [35.]
						# missing one arg
						push!(data_protocol,[1 0. 0 0.    0.     pulses  freq    true     "BCM_$(age)_$(freq)_$(pulses)"  -0.5 "LTD" "$paper" 25. 1e3 2.5e3 1.5 0 0. 3. age "yes" "no"  0. 0. 200. 0. 0.])
					end
				end
			end
		end

		########## Supp files
		if paper == "DudekBear92_temp"
			for freq in [.5; .75 ; collect(1. : .5 : 9.); collect(10. : 5. : 50.) ]
				for temp in collect(33:.5:37)
					push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "BCM_temp_$(freq)_$(temp)"  -0.5 "LTD" "$paper" temp 1e3 2.5e3 1.5 0 0. 0. 2. 35. "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		if paper == "DudekBear92_Ca"
			for freq in [.5;.75;collect(1. : .5 : 9.); collect(10. : 5. : 50.) ]
				for Ca in collect(1:.25:3.5)
					push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "BCM_Ca_$(freq)_$(Ca)"  -0.5 "LTD" "$paper" 35. 1e3 Ca*1e3 1.5 0 0. 0. 2. 35. "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		if paper == "DudekBear92_Mg"
			for freq in [.5;.75;collect(1. : .5 : 9.); collect(10. : 5. : 50.) ]
				for Mg in collect(.5:.1:2.)
					push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "BCM_Mg_$(freq)_$(Mg)"  -0.5 "LTD" "$paper" 35. 1e3 2.5*1e3 Mg 0 0. 0. 2. 35. "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		if paper == "DudekBear92_dist"
			for freq in [.5;.75;collect(1. : .5 : 9.); collect(10. : 5. : 50.) ]
				for dist in collect(50: 25. : 400)
					push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "BCM_dist_$(freq)_$(dist)"  -0.5 "LTD" "$paper" 35. 1e3 2.5*1e3 1.5 0 0. 0. 2. 35. "yes" "no" 0. 0. dist 0. 0.])
				end
			end
		end

		########## Supp files
		ages = collect(2.5:2.5:80)
		if paper == "DudekBear92-Age"
			for freq in [.5; .75 ; collect(1. : .25 : 6.) ; collect(6.5 : .5 : 9.); collect(10. : 5. : 50.) ]
				for pulses in [900]
					for age in ages
						push!(data_protocol,[1 0. 0 0.    0.     pulses  freq    true     "BCM_$(age)_$(freq)_$(pulses)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
					end
				end
			end
		end

		########## Supp files
		ages = collect(5. : 5. : 80.)
		if paper == "Blocking_age_control"
			for freq in [ 3.;30.]
				for age in ages
					push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "BCMAGEBLOCK_$(age)_$(freq)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		ages = collect(5. : 5. : 80.)
		if paper == "Blocking_yNMDA"
			for freq in [3. ; 30.]
				for age in ages
					push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "yNMDA_$(age)_$(freq)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		ages = collect(5. : 5. : 80.)
		if paper == "Blocking_oNMDA"
			for freq in [3. ; 30.]
				for age in ages
					push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "oNMDA_$(age)_$(freq)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		ages = collect(5. : 5. : 80.)
		if paper == "Blocking_yBaP"
			for freq in [3. ; 30.]
				for age in ages
					push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "yBaP_$(age)_$(freq)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		ages = collect(5. : 5. : 80.)
		if paper == "Blocking_oBaP"
			for freq in [3. ; 30.]
				for age in ages
					push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "oBaP_$(age)_$(freq)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		ages = collect(5. : 5. : 80.)
		if paper == "Blocking_yGABA"
			for freq in [3. ; 30.]
				for age in ages
					push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "yGABA_$(age)_$(freq)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		########## Supp files
		ages = collect(5. : 5. : 80.)
		if paper == "Blocking_oGABA"
			for freq in [3. ; 30.]
				for age in ages
					push!(data_protocol,[1 0. 0 0.    0.     900  freq    true     "oGABA_$(age)_$(freq)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		##########  Figure 5
		if paper == "DudekBear92_BCM_recovery"
			for freq in [.1;.5;.75;collect(1. : 3. : 9.); collect(10. : 5. : 25.); collect(30. : 10. : 50.) ]
				for age in collect(5: 5. : 80) #reviewed
						push!(data_protocol,[2 48. 0 0.    0.     900  freq    true     "BCMduplet_$(age)_$(freq)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
						push!(data_protocol,[3 48. 0 0.    0.     900  freq    true     "BCMtriplet_$(age)_$(freq)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
						push!(data_protocol,[4 48. 0 0.    0.     900  freq    true     "BCMquplet_$(age)_$(freq)"  -0.5 "LTD" "$paper" 35. 1e3 2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		##########  Figure 5
		ages = collect(7:4:56)
		if paper == "DudekBear93-LFS"
			for age in ages
				push!(data_protocol,[1 0. 0 0.    0.     900  1.    true     "LFS_$age"	 -0.5 "LTD" "$paper" 35. 2e3  2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
			end
		end

		##########  Figure 5
		if paper == "DudekBear93-TBS"
			for age in collect(2.5:2.5:80)
                # one arg too many
				push!(data_protocol,[1 0. 0 0.    0.     4  100.    true     "TBS_$(age)"	 -0.5 "LTD" "$paper" 34. 1e3  2.5e3 1.5 10 200. 0. 0. 5. age "yes" "no" 0. 0. 200. 0. 0.])
# suggestion #	push!(data_protocol,[1 0. 0 0.    0.     4  100.    true     "TBS_$(age)"	 -0.5 "LTD" "$paper" 34. 1e3  2.5e3 1.5 10 200. 0. 5. age "yes" "no" 0. 0. 200. 0. 0.])
			end
		end

		##########  not shown in the paper
		if paper == "Cao-TBS"
			for age in collect(5:5:65)
				push!(data_protocol,[1 0. 0 0.    0.     4  100.    true     "TBS_$(age)"	 -0.5 "LTD" "$paper" 33 1e3  2.5e3 1.3 10 200. 0. 5. age "yes" "no" 0. 0. 200. 0. 0.])
			end
		end

		##########  Figure 5
		if paper == "RecoverLTD"
			for age in collect(30:5:70)
				push!(data_protocol,[1 0.  0 0.    0.     300  1.    true      "1Pre"	 -0.5 "LTD" "$paper" 35. 2e3  2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[2 50. 0 0.    0.     300  1.    true     "2Pre50"	 -0.5 "LTD" "$paper" 35. 2e3  2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[3 50. 0 0.    0.     300  1.    true     "3Pre50"	 -0.5 "LTD" "$paper" 35. 2e3  2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[4 50. 0 0.    0.     300  1.    true     "4Pre50"	 -0.5 "LTD" "$paper" 35. 2e3  2.5e3 1.5 0 0. 0. 2. age "yes" "no" 0. 0. 200. 0. 0.])
			end
		end

		##########  methods
		if paper == "Chang19"
			push!(data_protocol,[1 0. 0 0. 0.  30  0.49    true     "Chang_19"	 	 	 	-0.5 "LTD" "$paper" 25. 0.  1e3   2.  0 0. 0. 0. 4.  "no"  "no" 0. 0. 200. 0. 0.])
			push!(data_protocol,[1 0. 0 0. 0.  30  0.49    true     "Chang_19"	 	 	 	-0.5 "LTD" "$paper" 35. 0.  1e3   2.  0 0. 0. 0. 4.  "no"  "no" 0. 0. 200. 0. 0.])
		end

		##########  methods
		if paper == "Fujii_CaN"
			push!(data_protocol,[1 0. 0 0.    0.   100  20.    true     "FujiBito"  -0.5 "LTD" "$paper" 25. 0e3 2e3 .0001 0 0. 0.  0. 12. "no" "no"  0. 0. 200. 0. 0.])
			push!(data_protocol,[1 0. 0 0.    0.   100  20.    true     "FujiBito"  -0.5 "LTD" "$paper" 35. 0e3 2e3 .0001 0 0. 0.  0. 12. "no" "no"  0. 0. 200. 0. 0.])
		end

		##########  Figure 6
		if paper == "YannisDebanne20"
			for Ca in collect(1.1 : .1 : 3.1)
				for delay in [2.5;collect(10.:10.:100.)]
					  push!(data_protocol,[1 0.  1 0. abs(delay)  100  0.3	true   string("$(paper)_$(delay)_$(Ca)") 0. "LTP" "$(paper)"  30.45 .5e3 Ca*1e3 Ca/1.5 0 0. 0. 6. 21 "no" "no" 0. 0. 200. 0. 0.])
					  push!(data_protocol,[1 0.  1 0. abs(delay)  150  0.3	false  string("$(paper)_-$(delay)_$(Ca)") 0. "LTD" "$(paper)" 30.  .5e3 Ca*1e3 Ca/1.5 0 0. 0. 6. 14 "no" "no" 0. 0. 200. 0. 0.])
				end
					  push!(data_protocol,[1 0.  1 0. abs(25)  150  0.3	false  string("$(paper)_-25_$(Ca)") 0. "LTD" "$(paper)" 30.  .5e3 Ca*1e3 Ca/1.5 0 0. 0. 6. 14 "no" "no" 0. 0. 200. 0. 0.])
			end
		end

		##########  Supp files
		if paper == "YannisDebanne_temp"
			for temp in collect(29:.5:36)
				for delay in [2.5;collect(10.:10.:100.)]
					  push!(data_protocol,[1 0.  1 0. abs(delay)  100  0.3	true   string("$(paper)_$(delay)_$(temp)") 0. "LTP" "$(paper)"   temp .5e3 2.5*1e3 2.5/1.5 0 0. 0. 6. 21 "no" "no" 0. 0. 200. 0. 0.])
					  push!(data_protocol,[1 0.  1 0. abs(delay)  150  0.3	false  string("$(paper)_-$(delay)_$(temp)") 0. "LTD" "$(paper)"  temp .5e3 2.5*1e3 2.5/1.5 0 0. 0. 6. 14 "no" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		##########  Supp files
		if paper == "YannisDebanne_age"
			for age in collect(14:.5:21)
				for delay in [2.5;collect(10.:10.:100.)]
					  push!(data_protocol,[1 0.  1 0. abs(delay)  100  0.3	true   string("$(paper)_$(delay)_$(age)") 0. "LTP" "$(paper)"   30.45 .5e3 2.5*1e3 2.5/1.5 0 0. 0. 6. age "no" "no" 0. 0. 200. 0. 0.])
					  push!(data_protocol,[1 0.  1 0. abs(delay)  150  0.3	false  string("$(paper)_-$(delay)_$(age)") 0. "LTD" "$(paper)"  30.   .5e3 2.5*1e3 2.5/1.5 0 0. 0. 6. age "no" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		##########  Supp files
		if paper == "Meredith03-GABA"
			for age in  collect(5: 5. : 20)
				push!(data_protocol,[1 0.  4 10.  10.	 30  .2    true     "1Pre4Post10_$(age)" 0.85 "LTP" "$paper"  28.  1e3   2e3 2. 0 0. 0. 4.   age  "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  3 10.  10.	 30  .2    true     "1Pre3Post10_$(age)" 0.85 "LTP" "$paper"  28.  1e3   2e3 2. 0 0. 0. 4.   age  "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10.  10.	 30  .2    true     "1Pre2Post10_$(age)" 0.85 "LTP" "$paper"  28.  1e3   2e3 2. 0 0. 0. 4.   age  "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.   10.    30  .2    true     "1Pre1Post10_$(age)" 0.85  "NC" "$paper"  24.5  1e3   2e3 2. 0 0. 0. 2.  age  "yes" "no" 0. 0. 200. 0. 0.])
			end
			for age in  collect(25. : 5 : 60)
				push!(data_protocol,[1 0.  4 10.  10.	 30  .2    true     "1Pre4Post10_$(age)" 0.85 "LTP" "$paper"  27.  1e3  2e3 2. 0 0. 0. 4.  age  "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  3 10.  10.	 30  .2    true     "1Pre3Post10_$(age)" 0.85 "LTP" "$paper"  27.  1e3  2e3 2. 0 0. 0. 4.  age  "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10.  10.	 30  .2    true     "1Pre2Post10_$(age)" 0.85 "LTP" "$paper"  27.  1e3  2e3 2. 0 0. 0. 4.  age  "yes" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.   10.    30  .2    true     "1Pre1Post10_$(age)" 0.85  "NC" "$paper"  25.  1e3  2e3 2. 0 0. 0. 2.  age  "yes" "no" 0. 0. 200. 0. 0.])
			end
		end

		##########  Supp files
		if paper == "TigaretvsMeredith"
			for age in  collect(5: 5. : 60)
				push!(data_protocol,[1 0.  4 10    10.	  300 5.    true     "1Pre4Post10_$(age)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2. age "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  3 10    10.	  300 5.    true     "1Pre3Post10_$(age)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2. age "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10    10.	  300 5.    true     "1Pre2Post10_$(age)" 0.75 "LTP" "$paper" 35. 2e3  2.5e3 1.3 0 0. 0. 2. age "no" "yes" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    10.    300 5.    true     "1Pre1Post10_$(age)" 0.1  "NC"  "$paper" 35. 1e3  2.5e3 1.3 0 0. 0. 2. age "no" "yes" 0. 0. 200. 0. 0.])
			end
		end

		##########  Supp files
		if paper == "YannisvsMeredith"
			for age in  collect(5: 5. : 60)
				push!(data_protocol,[1 0.  4 10. 10.  100  0.3	true    "1Pre4Post10_$(age)" 0. "LTP"  "$(paper)"  30. .5e3 1.8*1e3 1.8/1.5 0 0. 0. 6. age "no" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  3 10. 10.  100  0.3	true    "1Pre3Post10_$(age)" 0. "LTP"  "$(paper)"  30. .5e3 1.8*1e3 1.8/1.5 0 0. 0. 6. age "no" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10. 10.  100  0.3	true    "1Pre2Post10_$(age)" 0. "LTP"  "$(paper)"  30. .5e3 1.8*1e3 1.8/1.5 0 0. 0. 6. age "no" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.  10.  100  0.3	true    "1Pre1Post10_$(age)" 0. "LTP"  "$(paper)"  30. .5e3 1.8*1e3 1.8/1.5 0 0. 0. 6. age "no" "no" 0. 0. 200. 0. 0.])
			end
		end

		##########  Supp files
		if paper == "WittenbergWang06_D"
			sepp=10.
				push!(data_protocol,[1 0.  1 0.    	   0.        100 5.    true    "1Pre-1Post0"      0.0  "NC"  "$paper" 22.75   1.8e3 2e3 1.0 0 0. 0. 3. 21. "no" "no" 0. 0. 200. 0. 0.])
			for delay in collect(10.: sepp: 30.)
				push!(data_protocol,[1 0.  1 0.    	   delay     100 5.    true    "1Pre-1Post$(delay)" 0.0  "NC"  "$paper" 22.75   1.8e3 2e3 1.0 0 0. 0. 3. 21. "no" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    	   delay     100 5.    false   "1Post-1Pre$(delay)" 0.0  "NC"  "$paper" 22.75   1.8e3 2e3 1.0 0 0. 0. 3. 21. "no" "no" 0. 0. 200. 0. 0.])
			end
			for delay in collect(40.: sepp: 50.)
				push!(data_protocol,[1 0.  1 0.    	   delay     100 5.    true    "1Pre-1Post$(delay)" 0.0  "NC"  "$paper" 22.75   1.8e3 2e3 1.0 0 0. 0. 3. 21. "no" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    	   delay     100 5.    false   "1Post-1Pre$(delay)" 0.0  "NC"  "$paper" 22.75   1.8e3 2e3 1.0 0 0. 0. 3. 21. "no" "no" 0. 0. 200. 0. 0.])
			end
			for delay in collect(60.: sepp: 100.)
				push!(data_protocol,[1 0.  1 0.    	   delay     70 5.     true    "1Pre-1Post$(delay)" 0.0  "NC"  "$paper" 23.  1.8e3 2e3 1.0 0 0. 0. 3. 21. "no" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  1 0.    	   delay     70 5.     false   "1Post-1Pre$(delay)" 0.0  "NC"  "$paper" 23.  1.8e3 2e3 1.0 0 0. 0. 3. 21. "no" "no" 0. 0. 200. 0. 0.])
			end
		end

		##########  Supp files
		if paper == "WittenbergWang06_B"
			sepp=10.
			for delay in collect(10.: sepp: 30.)
				push!(data_protocol,[1 0.  2 10. delay-10     100 5.    true    "1Pre-2Post$delay" 0.0  "NC"  "$paper"   25. 1.8e3 2e3 1.0 0 0. 0. 3. 14. "no" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10. delay        100 5.    false   "2Post-1Pre$delay" 0.0  "NC"  "$paper"   22. 1.2e3 2e3 1.0 0 0. 0. 3. 21. "no" "no" 0. 0. 200. 0. 0.])
			end
			for delay in collect(40.: sepp: 50)
				push!(data_protocol,[1 0.  2 10. delay-10     70 5.    true    "1Pre-2Post$delay" 0.0  "NC"  "$paper"   22. 1.9e3 2e3 1.0 0 0. 0. 3. 18. "no" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10. delay        70 5.    false   "2Post-1Pre$delay" 0.0  "NC"  "$paper"   23 1.8e3 2e3 1.0 0 0. 0. 3. 18. "no" "no" 0. 0. 200. 0. 0.])
			end
			for delay in collect(60.: sepp: 100.)
				push!(data_protocol,[1 0.  2 10. delay-10     100 5.    true    "1Pre-2Post$delay" 0.0  "NC"  "$paper"   23 1.8e3 2e3 1.0 0 0. 0. 3. 18. "no" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10. delay        100 5.    false   "2Post-1Pre$delay" 0.0  "NC"  "$paper"   23 1.8e3 2e3 1.0 0 0. 0. 3. 18. "no" "no" 0. 0. 200. 0. 0.])
			end
				push!(data_protocol,[1 0.  2 10.  100.        100 5.    true    "1Pre-2Post100"    0.0  "NC"  "$paper"   23 1.8e3 2e3 1.0 0 0. 0. 3. 18. "no" "no" 0. 0. 200. 0. 0.])
		end

		##########  Supp files
		if paper == "WittenbergWang06_P"
			sepp=10.
			for delay in collect(10.: sepp: 100.)
				push!(data_protocol,[1 0.  2 10. delay-10     30 5.    true    "1Pre-2Post$delay" 0.0  "NC"  "$paper"   25.5    1.8e3 2e3 1.0 0 0. 0. 3. 14. "no" "no" 0. 0. 200. 0. 0.])
				push!(data_protocol,[1 0.  2 10. delay        30 5.    false   "2Post-1Pre$delay" 0.0  "NC"  "$paper"   25.5    1.8e3 2e3 1.0 0 0. 0. 3. 14. "no" "no" 0. 0. 200. 0. 0.])

			end
				push!(data_protocol,[1 0.  2 10.  100.        30 5.    true    "1Pre-2Post100"    0.0  "NC"  "$paper"   26. 1.8e3 2e3 1.0 0 0. 0. 3. 14. "no" "no" 0. 0. 200. 0. 0.])
		end

		##########  Supp files
		if paper == "Mizuno01-LTP-Mg"
			for Mg in [.0001]#; collect(.1:.1:1.);collect(1.2: .2:2.0)]
				for pulses in [1; 3; 5; 10; 25; 50; 100; 150]
						# missing one arg
						push!(data_protocol,[1 0. 0 0.  0.     pulses  1.    true     "LTP_$(pulses)_$(Mg)"	 -0.5 "LTD" "$paper" 26. 3e3 2.4e3 Mg 0 0. 2. 12 "yes" "no" 0. 0. 200. 0. 0.])
# suggestion #			push!(data_protocol,[1 0. 0 0.  0.     pulses  1.    true     "LTP_$(pulses)_$(Mg)"	 -0.5 "LTD" "$paper" 26. 3e3 2.4e3 Mg 0 0. 0. 2. 12 "yes" "no" 0. 0. 200. 0. 0.])
				end
			end
		end

		##########  Supp files
		if paper == "Mizuno01LTP"
			for pulses in [ 1; 3; 5; 10 ;  25; 50 ; 60 ; 100; 120 ;200; 600; 900  ]
					push!(data_protocol,[1 0. 0 0.  0.     pulses  1.    true     "LTP_$(pulses)"	 -0.5 "LTD" "$paper" 26.5 3e3 2.4e3 0.0001 0 0. 0. 2. 12 "yes" "no" 0. 0. 200. 0. 0.])
			end
		end

		##########  Supp files
		if paper == "Mizuno01LTD"
			for pulses in [ 1; 3; 5; 10 ;  25; 50 ; 60 ; 100; 120 ;200; 600; 900  ]
				push!(data_protocol,[1 0. 0 0.  0.     pulses  1.    true     "LTD_$(pulses)"	 -0.5 "LTD" "$paper" 31. 3e3 2.4e3 0.0001 0 0. 0. 2. 12 "yes" "no" 0. 0. 200. 0. 0.])
			end
		end

		return data_protocol
end
#########################################
function PoissonProcess(rate, start_time, end_time)
	post = [start_time]
		while true
			#ISI is exponential(rate)
			dt = rand(Exponential(1000/rate))
			if post[end] + dt > end_time
				break
			else
				push!(post, post[end] + dt)
			end
		end
		return post
end
#########################################
rates_adapt(a, b, c, d, Ca) = a * b /(c + d * Ca)
########################### Passing the var names to the sim outcome
function passing(s, XC)
	for i in 1:length(s)
		sy=s[i]
		value=XC[i,:]
		@eval (($sy)=($value))
	end
end
#########################################
function integrate(tt, vals)
	res = zero(eltype(tt))
	for i in 1:(length(tt)-1)
		res = res + (tt[i+1]-tt[i]) *(vals[i+1] + vals[i]) / 2
	end
	return res
end
##########  Compilation firing events  ##########
function paired_representation(start_time,n_pre,delay_pre,n_pos,delay_pos,delay,freq,repetion,causal,repeat_times,repeat_after)
	post_pulses = Float64[]
	pre_pulses = Float64[]
	step=1.0/((freq)*0.001)
  for i in collect(start_time:step:(start_time+step*(repetion-1)))
	append!(pre_pulses,ifelse(causal==true,0,delay) .+ i .+ range(0., length=n_pre, stop=delay_pre*ifelse(n_pre>1,n_pre-1,n_pre)))
	append!(post_pulses,ifelse(causal==false,0,delay) .+ i .+ range(0., length=n_pos, stop=delay_pos*ifelse(n_pos>1,n_pos-1,n_pos)))
  end
  if repeat_times>0 && repetion>0
		  	pre_pulses_aux = pre_pulses .-start_time
			for i in 2:repeat_times
		  		append!(pre_pulses,pre_pulses_aux .+ repeat_after .+ pre_pulses[end])
			end
			post_pulses_aux = post_pulses .-start_time
			for i in 2:repeat_times
				append!(post_pulses,post_pulses_aux .+ repeat_after .+ post_pulses[end])
			end
  end
  return pre_pulses, post_pulses
end

"""
$(SIGNATURES)

Generate external stimulation times and indices for pre/post stimulation. Usually used with `dataProtocol` (see folder examples in `examples/`).

# Output
- `event_times` sorted list of external event times
- `is_pre_or_post_index` whether the external events are pre (`true`) or post (`false`)
"""
function firingPattern(  ;start_time::Float64 = 1.,
						n_pre::Int64 = 1, delay_pre::Float64 = 1.,
						n_pos::Int64 = 1, delay_pos::Float64 = 1.,
						delay::Float64 = 1.,
						pulse::Int64 = 1,
						freq::Float64 = 1.,
						causal::Bool,
						repeat_times::Int64 = 0,
						repeat_after::Float64 = 0.)
	pre_pulses,post_pulses = paired_representation(start_time,n_pre,delay_pre,n_pos,delay_pos,delay,freq,pulse,causal,repeat_times,repeat_after)
	# specify the dendritic BAP times
	nglus = length(pre_pulses)
	nbaps = length(post_pulses)
	#%% ###########################################################################
	# Concatenate event trains
	all_events_times = vec([pre_pulses; post_pulses])
	all_events_sort = sort(all_events_times)
	# index of sorted events
	all_events_inds = sortperm(all_events_times)
	is_pre_or_post_index = [convert(Bool, e > length(pre_pulses)) for e in all_events_inds]
	return all_events_sort, .!(is_pre_or_post_index)
end

function initial_conditions_continuous(param_synapse)
	@unpack_SynapseParams param_synapse
		return vec([
		E_leak			#Vsp,
		E_leak			#Vdend,
		E_leak			#Vsoma,
		1.				#λ,
		0.				#ImbufCa,
		Ca_infty		#Ca,
		0.				#Dye,
		CaM_con			#CaM0
		0.				#CaM2C
		0.				#CaM2N
		0.				#CaM4
		mCaN_con		#mCaN
		0.				#CaN4
		mKCaM_con		#mKCaM
		0.				#KCaM0
		0.				#KCaM2N
		0.				#KCaM2C
		0.				#KCaM4
		0.				#PCaM0
		0.				#PCaM2C
		0.				#PCaM2N
		0.				#PCaM4
		0.				#P
		0.				#P2
		0.				#LTD,
		0.				#LTP,
		0.				#LTD_act,
		0.				#LTP_act,
		0.				#m,
		0.				#h,
		0.				#n
		0.				#SK
		1.
		1.
	  ])
	  end

function initial_conditions_deterministic(param_synapse)
	@unpack_SynapseParams param_synapse
		return vec([
		E_leak			#Vsp,
		E_leak			#Vdend,
		E_leak			#Vsoma,
		1.				#λ,
		0.				#ImbufCa,
		Ca_infty		#Ca,
		0.				#Dye,
		CaM_con			#CaM0
		0.				#CaM2C
		0.				#CaM2N
		0.				#CaM4
		mCaN_con		#mCaN
		0.				#CaN4
		mKCaM_con		#mKCaM
		0.				#KCaM0
		0.				#KCaM2N
		0.				#KCaM2C
		0.				#KCaM4
		0.				#PCaM0
		0.				#PCaM2C
		0.				#PCaM2N
		0.				#PCaM4
		0.				#P
		0.				#P2
		0.				#LTD,
		0.				#LTP,
		0.				#LTD_act,
		0.				#LTP_act,
		0.				#m,
		0.				#h,
		0.				#n
		0.				#SK
		1.
		1.
		N_ampa
		0.
		0.
		0.
		0.
		0.
		0.
		0.
		0.
		0.
		0.
		0.
		0.
		0.
		0.
		0.
	  	N_N2A
		0.
		0.
		0.
		0.
		0.
		0.
		N_N2B
		0.
		0.
		0.
		0.
		0.
		0.
	  	N_caR
		0.
		0.
		0.
		N_caT
		0.
		0.
		0.
	  	N_caL
		0.
		0.
		N_GABA
		0.
		0.
		0.
		0.
		100.
		0.
		0.
	  ])
	  end

function initial_conditions_continuous_temp(param_synapse)
	@unpack_SynapseParams param_synapse
	if temp_rates <= 25
		return vec([
		-70.10245699808998
		-70.02736715107497
		-70.01992573979436
		1.0
		5.251484030952095
		0.17942311304488254
		0.0
		18.422417385628144
		0.061182835845181506
		0.007491230194287401
		2.5601235159798556e-5
		18.60768149870139
		1.3923185018623478
		47.943467158481965
		0.33484901412981444
		0.00873078540237604
		0.26160543666353603
		0.006945031752908706
		5.28369636945184
		4.008117080305045
		0.11294478061086997
		0.09967594670655765
		4.656494469862573
		7.283473926518337
		0.0
		0.0
		0.0
		0.0
		0.012283139643655655
		0.9999998289470913
		0.00010811866849049202
		0.09878906052566663
		1.0
		1.0])
	end
	if 25 < temp_rates <= 30.
		return vec([
		-70.0140727673961
		-70.00177103943689
		-70.00007589726667
		1.0
		3.48177628683147
		0.11706193391209137
		0.0
		19.902059159951726
		0.04813075789246737
		0.0033156783870759324
		1.3058864642411135e-5
		19.473469645955362
		0.5265303542539336
		49.658388018645105
		0.39979810740476907
		0.005987535794514187
		0.22114161049072126
		0.003443459770906218
		5.486877056572574
		3.2429635679748126
		0.0967411298498974
		0.06299852278885718
		4.315666798657558
		6.505994192089851
		0.0
		0.0
		0.0
		0.0
		0.0123167984282587
		0.9999998268859568
		0.00010833059621284
		0.016112012029112062
		1.0
		1.0
		])
	end
	if temp_rates > 30
		return vec([
		-70.02953996060384
		-70.00364683510847
		-70.0013995913228
		1.0
		3.3989773357494646
		0.11430174346528181
		0.0
		25.14694054430742
		0.048232611580821316
		0.004340152415081933
		6.8123069384392556e-6
		19.865575695641734
		0.13442430463415694
		62.32357049746406
		0.8459878113145832
		0.009234347120412881
		0.3607751119415231
		0.003937688645035126
		2.3440067259757305
		1.0593952833042108
		0.02913354160885602
		0.013585064852268125
		1.6677561415335218
		1.3426177862569586
		0.0
		0.0
		0.0
		0.0
		0.012314741437834478
		0.999999826986425
		0.00010831999489540837
		0.03393795777123425
		1.0
		1.0
		])
	end
end

function initial_conditions_continuous_steady_dye(param_synapse)
	@unpack_SynapseParams param_synapse
  return vec([
  -70.79429490923671
  -70.79411507925752
  -70.79373452702642
    1.0
    1.8937985515700873
    0.06544565489559395
   -9.353645289760924e-24
   27.668652955744996
    0.018289183656159797
    0.0012380953924366839
    6.019820446572491e-7
   19.97204506487
    0.027955066618653265
   63.581926032844144
    1.153523449658415
    0.003960046076090567
    0.16522683320931195
    0.0005534077298222614
    0.8323722641818724
    0.12411271177445717
    0.003138175023134278
    0.00047491361819079867
    1.6179504762996062
    0.24615855264482245
    0.0062839007098402155
    0.000982621969867967
    1.953017871747061
    0.30134128704713425
    0.0077134827679826694
    0.0012279303271366144
    0.0
    0.0
    0.0
    1.3854633395928156e-26
    3.243117778864661e-6
    0.9999826149121529
    0.0008260666028628115
    0.0002458563962436573
	])
end

function initial_conditions_discrete(param_synapse)
	@unpack_SynapseParams param_synapse

  return vec([
		  N_ampa,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,   #AMPA		 1-16
		  N_N2A,0,0,0,0,0,0,					  #NMDA		 17-23
		  0,									  #print	 24
		  0, 0, N_caR, 0,			   		      #R-type	 25-28
		  0, 0, N_caT, 0,		   			      #T-type	 29-32
		  N_caL, 0, 0,							  #L-type	 33-35
		  rest_plstcty, 0, 0,
		  N_N2B, 0, 0, 0, 0, 0, 0,				  #GLUN2B
		  N_GABA, 0, 0, 0, 0 ])					  #GABA
end

function initial_conditions_ds(param_synapse)
	@unpack_SynapseParams param_synapse

  return vec([
		  N_ampa,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,   #AMPA		 1-16
		  N_N2A,0,0,0,0,0,0,					  #NMDA		 17-23
		  0,									  #print	 24
		  0, 0, N_caR, 0,			   		      #R-type	 25-28
		  0, 0, N_caT, 0,		   			      #T-type	 29-32
		  N_caL, 0, 0,							  #L-type	 33-35
		  rest_plstcty, 0, 0,					  #PMC
		  N_N2B, 0, 0, 0, 0, 0, 0,				  #GLUN2B
		  N_GABA, 0, 0, 0, 0,					  #GABA
		  N_SK, 0 ])							  #SK
end
