"""
$(TYPEDEF)

# Presynaptic parameters

* Firing events are processed separately from the main simulation (at `src/OnlyStp.jl`) it takes the firing structure as input from the function [`firingPattern`](@ref).
* Using the presynaptic stimulation times, the vesicle depletion and AP induced by EPSP are estimated, however one can use a tag to deactivate it (covering sub-threshold EPSP cases) as in  [`dataProtocol`](@ref).
* The presynaptic part of our model is phenomenological, for instance, the variable `Soma` in `src/OnlyStp.jl:38` was made to represent the voltage and can accumulate for faster frequencies but has an abstract unit.

# Equations

based on D. Sterratt et al book; [`Principles of Computational Modelling in Neuroscience`](https://www.compneuroprinciples.org/ ) page 183 \n
raterec  = (`R_0` - `R`) ⋅ τ_R ⋅ `rrp` \n
raterrp  = (`D_0` - `D`) ⋅ τ_D ⋅ `rec` \n
rateref  = (`R_0` - `R`) ⋅ ref_dt \n

# Fields

$(FIELDS)
"""
@with_kw struct PreSynapseParams
	"recovery constant of pre calcium decay function"
	τ_rec 		 ::Float64		= 20000
	"fraction of decay constant of pre calcium decay f"
	δ_ca		 ::Float64		= .0004
	"decay time constant of pre calcium"
	τ_pre 		 ::Float64		= 20
	"decay time constant for AP induced by EPSP"
	τ_V		 	::Float64		= 40
	"delay to EPSPs onset and evoked AP"
	δ_delay_AP	::Float64		= 15.
	"initial conditions ready releaseble pool"
	D_0		 	::Int64		= 25
	"initial conditions recovery pool"
	R_0		 	::Int64		= 30
	"rate for `D -> R`"
	τ_R 		::Float64		= 5000
	"rate for `R -> D`"
	τ_D 		::Float64		= 45000
	"rate for `infinite reservoir -> R`"
	τ_R_ref 	::Float64		= 40000
	"sigmoid parameter for release probability"
	s 		    ::Float64		= 2.
	"sigmoid parameter for release probability"
	h 		    ::Float64		= .7 #This value changes given Ca external concentration
	"sampling rate for plotting / printing"
	sampling_rate ::Float64 	= 1.
end

"""
$(TYPEDEF)

Postsynaptic parameters.

---

# Units
- time: ms
- length: µm (area µm^2, volume µm^3)
- voltage: mV
- current: pA
- conductance: nS
- resistance: GOhm
- capacitance: pF ( ms.pA/mV = ms.nS = ms/GOhm)
- concentration: µM

---

# Fields

$(FIELDS)
"""
@with_kw struct SynapseParams{Tp}
	"**POLYGONAL THRESHOLDS:** Region LTP"
	"Region LTP"
	LTP_region ::Tp 	=  VPolygon([[6.35,1.4],[10,1.4],[6.35,29.5],[10,29.5]])
	# LTP_region ::Tp 	= VPolygon([SVector(6.35,1.4), SVector(10,1.4),SVector(6.35,29.5), SVector(10,29.5)])
	"Region LTD"
	LTD_region ::Tp	    = VPolygon([[6.35,1.4],[6.35,23.25],[6.35,29.5],[1.85,11.327205882352938],[1.85,23.25],[3.7650354609929075,1.4],[5.650675675675676,29.5]])
	# LTD_region ::Tp		 = VPolygon([SVector(6.35,1.4),SVector(6.35,23.25),SVector(6.35,29.5),SVector(1.85,11.327205882352938),SVector(1.85,23.25),SVector(3.7650354609929075,1.4),SVector(5.650675675675676,29.5)])
	############################################################################
	"**ACTIVATION RATES FOR THRESHOLDS**"
	a_D             ::Float64				         = .1
	b_D             ::Float64				         = .00002
	a_P             ::Float64				         = .2
	b_P             ::Float64				         = .0001
	t_D             ::Float64				         = 18000
	t_P             ::Float64				         = 13000

	"**Sigmoids controlling the rate of plasticity change**"
	K_D  			::Float64 	                     = 80000.
	K_P  			::Float64 	                     = 13000.

	"**Plasticity states**"
	rest_plstcty     ::Int64			             = 100

	############################################################################
	"**SIMULATION**"
	t_end            ::Float64 		                 = 100
	sampling_rate    ::Float64 		                 = 10

	############################################################################
	"**BIOPHYSICAL AND GHK PARAMETERS**"

	temp_rates       ::Float64			             = 35.
	age              ::Float64				         = 60.
	faraday          ::Float64			             = 96485e-6*1e-3
	Ca_ext           ::Float64				         = 2.5e3
	Ca_infty         ::Float64 	 	                 = 50e-3
	tau_ca           ::Float64 	 	                 = 10.
	D_Ca             ::Float64 	 	                 = .3338
	f_Ca             ::Float64 	 	                 = .1
	perm             ::Float64				         = -0.04583333333333333
	z                ::Float64					     = 2.
	gas              ::Float64			             = 8.314e-6
	p_release 		 ::NTuple{4,Float64} 			 = (0.004225803293622208 ,1708.4124496514878 , 1.3499793762587964, 0.6540248201173222)

	############################################################################
	"**BACKPROPAGATION ATTENUATION**"
	trec             ::Float64				         = 2000
	trec_soma        ::Float64					     = 500
	delta_decay      ::Float64		                 = 1.7279e-5
	p_age_decay_bap  ::NTuple{3,Float64} 			 = (0.13525468256031167 , 16.482800452454164 ,5.564691354645679)
	delta_soma ::Float64		 		        	 = 2.5e-5 * (p_age_decay_bap[3]/(1+exp(p_age_decay_bap[1]*(age-p_age_decay_bap[2]))))
	delta_aux   ::Float64		 		        	 = 2.304e-5
	injbap           ::Float64 			             = 2.
	soma_dist        ::Float64			             = 200.
	p_dist			 ::NTuple{4,Float64}  			 = (0.019719018173341547 , 230.3206470553394 ,1.4313810030893268 , 0.10406540965358434)
	ϕ_dist             ::Float64				         = ( p_dist[4] + p_dist[3]/(1+exp(p_dist[1]*(soma_dist-p_dist[2]))))
	I_clamp          ::Float64 	 	 	             = 0.0

	############################################################################
	"**Na and K**"
	gamma_Na         ::Float64			             = 8e2
	Erev_Na          ::Float64			             = 50.
	gamma_K          ::Float64			             = 4e1
	Erev_K           ::Float64				         = -90.

	############################################################################
	"**NMDAr temperature modification**"
	p_nmda_frwd		 ::NTuple{4,Float64}    		 = (-0.09991802053299291, -37.63132907014948, 1239.0673283348326, -1230.6805720050966  )
	frwd_T_chng_NMDA ::Float64					     = ( p_nmda_frwd[4] + p_nmda_frwd[3]/(1+exp(p_nmda_frwd[1]*(temp_rates-p_nmda_frwd[2]))))
	p_nmda_bcwd      ::NTuple{4,Float64} 			 = (-0.10605060141396823, 98.99939433046647, 1621.6168608608068, 3.0368551011554143 )
	bcwd_T_chng_NMDA ::Float64		 			     = ( p_nmda_bcwd[4] + p_nmda_bcwd[3]/(1+exp(p_nmda_bcwd[1]*(temp_rates-p_nmda_bcwd[2]))))
	#0.16031*temp_rates - 0.80775#
	"**NMDAr kinetics (GluN2A type)**"
	NMDA_N2A_ka      ::Float64		 	             = frwd_T_chng_NMDA * 34. * 1e-3			# (uM-1ms-1)
	NMDA_N2A_kb      ::Float64		 	             = frwd_T_chng_NMDA * 17.	 * 1e-3			# (uM-1ms-1)
	NMDA_N2A_kc      ::Float64		 	             = frwd_T_chng_NMDA * 127. * 1e-3			# (uM-1ms-1)
	NMDA_N2A_kd      ::Float64		 	             = frwd_T_chng_NMDA * 580. * 1e-3			# (uM-1ms-1)
	NMDA_N2A_ke      ::Float64		 	             = frwd_T_chng_NMDA * 2508. * 1e-3			# (uM-1ms-1)
	NMDA_N2A_kf      ::Float64		 	             = frwd_T_chng_NMDA * 3449. * 1e-3			# (uM-1ms-1)
	NMDA_N2A_k_f     ::Float64		 	             = bcwd_T_chng_NMDA * 662. * 1e-3			# (ms-1)
	NMDA_N2A_k_e     ::Float64		 	             = bcwd_T_chng_NMDA * 2167. * 1e-3			# (ms-1)
	NMDA_N2A_k_d     ::Float64		 	             = bcwd_T_chng_NMDA * 2610. * 1e-3			# (ms-1)
	NMDA_N2A_k_c     ::Float64		 	             = bcwd_T_chng_NMDA * 161. * 1e-3			# (ms-1)
	NMDA_N2A_k_b     ::Float64		 	             = bcwd_T_chng_NMDA * 120. * 1e-3			# (ms-1)
	NMDA_N2A_k_a     ::Float64		 	             = bcwd_T_chng_NMDA * 60. * 1e-3			# (ms-1)

	"**NMDAr kinetics (GluN2B type)**"
	NMDA_N2B_sa      ::Float64		 	             = frwd_T_chng_NMDA * .25 * 34. * 1e-3			# (uM-1ms-1)
	NMDA_N2B_sb      ::Float64		 	             = frwd_T_chng_NMDA * .25 * 17. * 1e-3			# (uM-1ms-1)
	NMDA_N2B_sc      ::Float64		 	             = frwd_T_chng_NMDA * .25 * 127. * 1e-3			# (uM-1ms-1)
	NMDA_N2B_sd      ::Float64		 	             = frwd_T_chng_NMDA * .25 * 580. * 1e-3			# (uM-1ms-1)
	NMDA_N2B_se      ::Float64		 	             = frwd_T_chng_NMDA * .25 * 2508. * 1e-3		# (uM-1ms-1)
	NMDA_N2B_sf      ::Float64		 	             = frwd_T_chng_NMDA * .25 * 3449. * 1e-3		# (uM-1ms-1)
	NMDA_N2B_s_f     ::Float64		 	             = bcwd_T_chng_NMDA * .23 * 662. * 1e-3			# (ms-1)
	NMDA_N2B_s_e     ::Float64		 	             = bcwd_T_chng_NMDA * .23 * 2167. * 1e-3		# (ms-1)
	NMDA_N2B_s_d     ::Float64		 	             = bcwd_T_chng_NMDA * .23 * 2610. * 1e-3		# (ms-1)
	NMDA_N2B_s_c     ::Float64		 	             = bcwd_T_chng_NMDA * .23 * 161. * 1e-3			# (ms-1)
	NMDA_N2B_s_b     ::Float64		 	             = bcwd_T_chng_NMDA * .23 * 120. * 1e-3			# (ms-1)
	NMDA_N2B_s_a     ::Float64		 	             = bcwd_T_chng_NMDA * .23 * 60. * 1e-3			# (ms-1)

	# NMDA details
	p_nmda			 ::NTuple{4,Float64} 			 = (0.004477162852447629 , 2701.3929349701334 , 58.38819453272428 , 33.949463268365555)
	gamma_nmda       ::Float64			             = (p_nmda[4] + p_nmda[3]/(1+exp(p_nmda[1]* (Ca_ext-p_nmda[2])))) * 1e-3
	p_age 			 ::NTuple{4,Float64} 			 = (0.09993657672916968, 25.102347872464193, 0.9642137892004939, 0.5075183905839776)
	"**Ratio N2B/N2A**"
	r_NMDA_age       ::Float64 	                     = rand(Normal(0,.05)) + p_age[4] + p_age[3]/(1+exp(p_age[1]* (age-p_age[2])))#0.5+1.6*exp(-age/16.66) + rand(Normal(0,.05))
	N_NMDA           ::Int64				         = 15
	N_N2B            ::Int64 	                     = round(N_NMDA*r_NMDA_age/(r_NMDA_age+1))
	N_N2A            ::Int64 	                     = round(N_NMDA/(r_NMDA_age+1)) 					# Using Sinclair ratio
	"**Other NMDAr parameters**"
	Erev_nmda        ::Float64 		                 = 0.
	Mg               ::Float64 		                 = 1.0

	############################################################################
	"**AMPAr temperature modification**"
	p_ampa_frwd 	 ::NTuple{3,Float64} 		 	= (-0.4737773089201679, 31.7248285571622, 10.273135485873242)
	frwd_T_chng_AMPA ::Float64	                     = ( p_ampa_frwd[3]/(1+exp(p_ampa_frwd[1]*(temp_rates-p_ampa_frwd[2]))))#temp_rates*0.78-18.7
	p_ampa_bcwd 	 ::NTuple{3,Float64}  			 = (-0.36705555170278986, 28.976662403966674, 5.134547217640794)
	bcwd_T_chng_AMPA ::Float64	                     = ( p_ampa_bcwd[3]/(1+exp(p_ampa_bcwd[1]*(temp_rates-p_ampa_bcwd[2]))))#temp_rates*0.37-8.25

	"** AMPAr kinetics\n"
	AMPA_k1          ::Float64			             = frwd_T_chng_AMPA*1.6*1e7*1e-6*1e-3 				# (uM-1ms-1)
	AMPA_k_1         ::Float64			             = bcwd_T_chng_AMPA*7400*1e-3						# (ms-1)
	AMPA_k_2         ::Float64			             = bcwd_T_chng_AMPA*0.41*1e-3						# (ms-1)
	AMPA_alpha       ::Float64			             = 2600*1e-3										# (ms-1)
	AMPA_beta        ::Float64			             = 9600*1e-3 										# (ms-1)
	AMPA_delta_1     ::Float64		                 = 1500*1e-3 										# (ms-1)
	AMPA_gamma_1     ::Float64		                 = 9.1*1e-3 										# (ms-1)
	AMPA_delta_2     ::Float64		                 = 170*1e-3 										# (ms-1)
	AMPA_gamma_2     ::Float64		                 = 42*1e-3 											# (ms-1)
	AMPA_delta_0     ::Float64		                 = 0.003*1e-3 										# (ms-1)
	AMPA_gamma_0     ::Float64		                 = 0.83*1e-3 										# (ms-1)

	"**AMPAr conductances**"
	gamma_ampa1      ::Float64		                 = .5*31e-3											# (nS)
	gamma_ampa2      ::Float64		                 = .5*52e-3											# (nS)
	gamma_ampa3      ::Float64		                 = .5*73e-3											# (nS)
	N_ampa           ::Int64 	 	                 = 120 												# AMPAr number [REF: ]
	Erev_ampa        ::Float64 		                 = 0.0 												# mV, AMPAR reversal potential

	############################################################################
	"**GABAr**"
	N_GABA		 	 ::Int64						 = 34
	p_Cl		 	 ::NTuple{4,Float64}  			 = (0.09151696057098718, 0.6919298240788684 ,243.5159017060495,-92.6496083089155)
	"**GABAr, Cl reversal potential**"
	Erev_Cl	 	 	 ::Float64					     = (p_Cl[4]+p_Cl[3]/(1+exp(p_Cl[1]*(age-p_Cl[2]))))
	gamma_GABA	 	 ::Float64					     = 35e-3
	GABA_r_b1 	     ::Float64			             = 1e6 * 1e-6 * 1e-3 * 20
	GABA_r_u1 	     ::Float64			             = 1e3 * 4.6e-3
	GABA_r_b2 	     ::Float64		                 = 1e6 * 1e-6 * 1e-3 * 10
	GABA_r_u2 	     ::Float64		                 = 1e3 * 9.2e-3
	GABA_r_ro1 	     ::Float64		                 = 1e3 * 3.3e-3
	GABA_r_ro2 	     ::Float64		                 = 1e3 * 10.6e-3
	p_GABA	 		 ::NTuple{4,Float64}   			= (0.19127068198185954 , 32.16771140618756,  -1.2798050197287802 , 1.470692263981145)
	GABA_r_c1 	     ::Float64		                 = (p_GABA[4] + p_GABA[3]/(1+exp(p_GABA[1]*(temp_rates-p_GABA[2])))) * 1e3 * 9.8e-3
	GABA_r_c2 	     ::Float64		                 = (p_GABA[4] + p_GABA[3]/(1+exp(p_GABA[1]*(temp_rates-p_GABA[2])))) * 400e-3


	############################################################################
	"**Passive electrical properties**"
	E_leak           ::Float64		                 = -70.0
	g_leak           ::Float64 		                 = 4e-6
	Cm               ::Float64 		                 = 0.6e-2
	R_a              ::Float64 		                 = 1e-2

	############################################################################
	"** MORPHOLOGY: Dendritic properties**"
	D_dend           ::Float64			             = 2. 				 								# (um) dendrite diameter
	L_dend           ::Float64 		                 = 1400				 								# (um) dendrite length - choosen to tune attenuation and however not modified in BaP adaptation for simplicity sake
	A_dend           ::Float64 		                 = 2 * pi * (D_dend / 2) * L_dend 	 				# (um^2) dendrite surface area [5000 gives dendrite input resistance of 200 MOhm]
	Vol_dend         ::Float64	 		             = pi * ((D_dend / 2)^2) * L_dend 	 				# (um^3) dendrite volume
	Cdend            ::Float64                       = Cm * A_dend 					 					# dendritic membrane capacitance
	CS_dend          ::Float64 	                     = pi*(D_dend/2).^2 								# (um2), Dendrite cross sectional area
	g_leakdend       ::Float64                       = g_leak * A_dend									# (nS)

	############################################################################
	"**MORPHOLOGY: Soma properties**"
	D_soma           ::Float64			             = 30 				 								# (um) soma diameter
	A_soma           ::Float64 		                 = pi*D_soma^2						 				# (um^2) soma surface area [5000 gives dendrite input resistance of 200 MOhm]
	Csoma            ::Float64                       = Cm * A_soma 					 					# soma membrane capacitance
	CS_soma          ::Float64 	                     = pi*(D_soma/2).^2 								# (um2), soma cross sectional area
	g_leaksoma       ::Float64                       = 15.#												# (nS)
	g_diff           ::Float64			             = D_dend / (4R_a)									# value subjected to modifications due to BaP adaptation implementation

	############################################################################
	"**Spine properties**"
	Vol_sp           ::Float64 		                 = 0.03 			 								# um^3, Spine head volume [REF: Bartol et al., 2015]
	A_sp             ::Float64 		                 = 4*pi*((3*Vol_sp)/(4*pi))^(2. /3.) 	 			# Spine head surface area
	Csp              ::Float64 		                 = Cm * A_sp 		 								# spine membrane capacitance
	g_leaksp         ::Float64                       = g_leak * A_sp 	 								# spine head leak conductance

	############################################################################
	"**Neck properties**"
	D_neck           ::Float64 		                 = 0.1 					 							# (um), Spine neck diameter [REF: Bartol et al., 2015]
	L_neck           ::Float64 		                 = 0.2 					 							# (um), Neck length [REF: Bartol et al., 2015]
	CS_neck          ::Float64 		                 = pi*(D_neck/2).^2 								# (um^2), Neck cross sectional area
	g_neck           ::Float64 		                 = CS_neck / (L_neck * R_a)
	tau_diff		 ::Float64						 = ((Vol_sp/(2*D_Ca*D_neck))+(L_neck^2/(2*D_Ca)))
	############################################################################
	"**SYNAPTIC GLUTAMATE TRANSIENT PARAMETERS**"
	glu              ::Float64						 = 1.0 												# where Glutamate is ON (1) or OFF (0)
	glu_width        ::Float64 	                     = 1.0 												# ms, 0.1 ms for synapse [REF: arbitrary]
	glu_amp          ::Float64 			             = 1e+3 											# 1 mM [REF: arbitrary]
	glu_cv			 ::Float64					     = .5												# 0.5 from Liu et al, Neuron, 1999
	############################################################################
	"**SK channels**"
	N_SK             ::Int64				 	     = 15												# Number of SK channels
	SK_gamma         ::Float64	                     = 10e-3 											# (nS) Ref: James Maylie 10.1113/jphysiol.2003.049072
	SK_Erev          ::Float64 		                 = -90 												# (mv) Ref: Supp files Mellor and Griffith 2016
	SK_gating_half   ::Float64 	                     = 0.33 											# (uM) Ref: Supp files Mellor and Griffith 2016
	SK_time          ::Float64			             = 6.3 												# (ms) Ref: Supp files Mellor and Griffith 2016
	SK_hill          ::Float64			             = 6												# (ms) Ref: Supp files Mellor and Griffith 2016

	p_SK_bcwd 		 ::NTuple{4,Float64}  			 = (0.09391588258147192 ,98.85165844770867 , -147.61669527876904 ,  149.37767054612135 )
	bcwd_SK 	     ::Float64	 	                 = (p_SK_bcwd[4] + p_SK_bcwd[3]/(1+exp(p_SK_bcwd[1]*(temp_rates-p_SK_bcwd[2]))))
	p_SK_frwd        ::NTuple{4,Float64} 			 = (-0.334167923607112 ,25.590920461511878 , 2.2052151559841193 , 0.005904170174699533)
	frwd_SK 		 ::Float64		                 = (p_SK_frwd[4] + p_SK_frwd[3]/(1+exp(p_SK_frwd[1]*(temp_rates-p_SK_frwd[2]))))

	############################################################################
	"**CaM, CaMKII and CaN**"
	"**Concentrations"
	CaM_con          ::Float64                       = 30.
	mKCaM_con        ::Float64                       = 70. 												# (uM) (100), [Feng et al, Brain Res 2011] (rename)
	mCaN_con         ::Float64 	                     = 20. 												# (uM) Lisman

	############################################################################
	"**Chang Pepke model - CaM reactions I**"
	kon_1C           ::Float64                       = 5e-3
	koff_1C          ::Float64                       = 50e-3
	kon_2C           ::Float64                       = 10e-3
	koff_2C          ::Float64                       = 10e-3
	kon_1N           ::Float64                       = 100e-3
	koff_1N          ::Float64                       = 2000e-3
	kon_2N           ::Float64                       = 200e-3
	koff_2N          ::Float64                       = 500e-3

	############################################################################
	"**Chang Pepke model - CaM reactions II**"
	kf_CaM0          ::Float64		                 = 3.8e-6
	kb_CaM0          ::Float64                       = 5.5e-3
	kf_CaM2C         ::Float64                       = 0.92e-3
	kb_CaM2C         ::Float64                       = 6.8e-3
	kf_CaM2N         ::Float64                       = 0.12e-3
	kb_CaM2N         ::Float64                       = 1.7e-3
	kf_CaM4          ::Float64                       = 30e-3
	kb_CaM4          ::Float64                       = 1.5e-3


	############################################################################
	"**Chang Pepke model - CaMKII reactions**"
	kon_K1C          ::Float64		                 = 44e-3
	koff_K1C         ::Float64		                 = 33e-3
	kon_K2C          ::Float64		                 = 44e-3
	koff_K2C         ::Float64		                 = 0.8e-3
	kon_K1N          ::Float64		                 = 76e-3
	koff_K1N         ::Float64		                 = 300e-3
	kon_K2N          ::Float64		                 = 76e-3
	koff_K2N         ::Float64		                 = 20e-3

	############################################################################
	"**Chang Pepke model - autophosphorilation**"
	p_camkii_q10     ::NTuple{4,Float64}   			 = (0.5118207068695309 ,  45.47503600542303 ,-161.42634157226917 , 162.1718925882677)
	q10              ::Float64			             = p_camkii_q10[4] + p_camkii_q10[3]/(1+exp(p_camkii_q10[1]*(temp_rates-p_camkii_q10[2])))					# change of temp to fit chang 35C
	k1               ::Float64					     = 12.6e-3
	k2               ::Float64					     = q10 * 0.33e-3 #q10 * 0.33e-3
	k3               ::Float64					     = 4 * q10 * 0.17e-3 #q10 * 0.17e-3
	k4               ::Float64					     = 4 * 0.041e-3
	k5               ::Float64					     = 4* q10 * 2 * 0.017e-3 #q10 * 2* 0.017e-3

	#############################################################################
	"**CaM-CaN reactions**"
	p_CaN_frwd       ::NTuple{4,Float64}  			 = (-0.29481489145354556 , 29.999999999999968 ,  0.15940019940354327,	  0.870299900298228)
	kcanf            ::Float64 		                 = (p_CaN_frwd[4] + p_CaN_frwd[3]/(1+exp(p_CaN_frwd[1]*(temp_rates-p_CaN_frwd[2])))) * 1.75e-2								# Ref: Quintana in 22C - 4.6e-2
	p_CaN_bcwd 		 ::NTuple{4,Float64}  			 = (-0.6833299932488973 , 26.277500129849113 , 0.7114524682690591 , 0.29037766196937326)
	kcanb            ::Float64 		                 = (p_CaN_bcwd[4] + p_CaN_bcwd[3]/(1+exp(p_CaN_bcwd[1]*(temp_rates-p_CaN_bcwd[2])))) * 2e-5									# Ref: Quintana in 22C - 1.2e-6

	############################################################################
	"**VGCCs**"

	p_frwd_VGCC 	 ::NTuple{4,Float64}  			 = (1.0485098341579628 , 30.66869198447378 , -0.3040010721391852 ,  2.5032059559264357)
	frwd_VGCC        ::Float64	                     = (p_frwd_VGCC[4] +p_frwd_VGCC[3]/(1+exp(p_frwd_VGCC[1]*(temp_rates-p_frwd_VGCC[2]))))
	p_bcwd_VGCC 	 ::NTuple{4,Float64}  			 = (-0.3302682317933842 ,36.279019647221226,3.2259761593440155,0.7298285671937866)
	bcwd_VGCC        ::Float64	                     = (p_bcwd_VGCC[4] +p_bcwd_VGCC[3]/(1+exp(p_bcwd_VGCC[1]*(temp_rates-p_bcwd_VGCC[2]))))

	Erev_CaT         ::Float64                       = 10.0 											# (mV), Calcium channel reversal potential
	Erev_CaR         ::Float64                       = 10.0 											# (mV), Calcium channel reversal potential
	Erev_CaL         ::Float64                       = 10.0 											# (mV), Calcium channel reversal potential
	gamma_CaT        ::Float64                       = 12e-3 											# (nS) Magee and Jonhston 95
	gamma_CaR        ::Float64                       = 17e-3 											# (nS) Magee and Jonhston 95
	gamma_CaL        ::Float64                       = 27e-3 											# (nS) Magee and Jonhston 95
	N_caT            ::Int64		                 = 3
	N_caR            ::Int64		                 = 3
	N_caL            ::Int64		                 = 3

	############################################################################
	"**Calcium dye and buffers**"
	EGTA_kf          ::Float64                       = 2.7e-3 											# (uM/ms), Zenisek et al 2003, Naraghi, 1997
	EGTA_kb          ::Float64                       = 0.18*EGTA_kf 									# (/ms), assuming Kd of 0.18uM, Naraghi, 1997
	EGTA_con         ::Float64                       = 0.0 												# (uM), 0.2 used by Tigaret et al for imaging, 200 for elecrophysiology
	BAPTA_kf         ::Float64                       = 0.45 											# (/uM/ms), Zenisek et al 2003, Naraghi, 1997
	BAPTA_kb         ::Float64                       = 0.176*BAPTA_kf 									# (/ms), assuming Kd of 0.176uM, Naraghi, 1997
	BAPTA_con        ::Float64                       = 0.0 												# (uM)
	Imbuf_k_on       ::Float64				         = 0.247 											# (/uM/ms), Bartol et al, 2015
	Imbuf_k_off      ::Float64			             = 0.524 											# (/ms), Bartol et al, 2015
	K_buff_diss      ::Float64			             = Imbuf_k_off/Imbuf_k_on
	Imbuf_con        ::Float64				         = 62 												# (uM), 76.7 Bartol et al, 2015
	Imbuf_con_dend   ::Float64			             = Imbuf_con*4

	############################################################################
	"**Calcium fluorescent dyes**"
	ogb1_kf          ::Float64                       = 0.8 												# (/ms), from Bartol et al 2015, assuming a [Ca] = 1uM.
	ogb1_kb          ::Float64                       = 0.16 											# (/ms), from Bartol et al 2015
	fluo4_kf         ::Float64                       = 0.8 												# (/ms), from Bartol et al 2015, assuming a [Ca] = 1uM.
	fluo4_kb         ::Float64                       = 0.24 											# (/ms), from Bartol et al 2015
	dye              ::Float64	                     = 0.
	fluo5f_kf        ::Float64                       = dye*0.01 										# (uM/ms), Zenisek et al J Neurosci 2003, Naraghi, 1997
	fluo5f_kb        ::Float64                       = dye*26*fluo5f_kf 								# from Yasuda et al 2004. Kd = 1.3uM.
	fluo5f_con       ::Float64                       = dye*200.0 										# (uM), concentration used by Tigaret et al
end
