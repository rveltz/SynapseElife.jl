# New synapse model (cian, started 23/03/2018) that replaces NMDA model with fully state-based one from Jahr and Stevens, plus three types of VGCCs (R-type, T-type and L-type), from Magee and Johhston (1995).
function F_synapse(xdot, pop_c, discrete_var, p_synapse::SynapseParams, t, events_bap, bap_by_epsp)
	# vector field used for the continuous variable

	@unpack_SynapseParams p_synapse
	Vsp, Vdend,	Vsoma,
	λ, ImbufCa,
	Ca,	Dye,
	CaM0, CaM2C, CaM2N, CaM4,
	mCaN, CaN4,
	mKCaM, KCaM0, KCaM2N, KCaM2C, KCaM4,
	PCaM0, PCaM2C, PCaM2N, PCaM4,
	P, P2,
	LTD, LTP, act_D, act_P,
	m, h, n, SK, λ_age, λ_aux = pop_c

	##### Stochastic channels/receptors
	n1_ampa	     = discrete_var[14]									#ampa subconductance 1
	n2_ampa	     = discrete_var[15]									#ampa subconductance 2
	n3_ampa	     = discrete_var[16]									#ampa subconductance 3
	n1_nmda_A	 = discrete_var[22]									#nmda subconductance 1
	n2_nmda_A	 = discrete_var[23]									#nmda subconductance 2
	n1_nmda_B	 = discrete_var[44]									#nmda subconductance 1
	n2_nmda_B	 = discrete_var[45]									#nmda subconductance 2
	n_car	     = discrete_var[28]									#vgcc-R opened state
	n_cat	     = discrete_var[32]									#vgcc-T opened state
	n_cal	     = discrete_var[34] + discrete_var[35]				#vgcc-L opened states
	n_gaba1	 	 = discrete_var[49]									#GABA opened state
	n_gaba2	   	 = discrete_var[50]									#GABA opened state


	##### plasticity prediction regions
	CaMKII	     = KCaM0 + KCaM2C + KCaM2N + KCaM4 + PCaM0 + PCaM2C + PCaM2N + PCaM4 + P + P2
	CaN	 	     = CaN4

	#### activation when it is inside the region
	# this following line allocates 74.779 ns (3 allocations: 144 bytes). The 2 lines count for 25% of the performance
	∂LTD	     = SVector(CaN, CaMKII) ∈ LTD_region
	∂LTP	     = SVector(CaN, CaMKII) ∈ LTP_region

	∂act_D     = a_D * ∂LTD - b_D * act_D * (1 - ∂LTD)
	∂act_P     = a_P * ∂LTP - b_P * act_P * (1 - ∂LTP)

	##### Na channel
	m_inf		 = alpha_m(Vsoma) / (alpha_m(Vsoma) + beta_m(Vsoma))
	m_tau		 = 1 / (alpha_m(Vsoma) + beta_m(Vsoma))
	∂m		     = (m_inf - m) / m_tau
	∂h		     = alpha_h(Vsoma) * (1 - h) - beta_h(Vsoma) * h
	I_Na		 = gamma_Na * (m^3) * h * (Erev_Na - Vsoma)


	##### K channel
	n_inf = 1 / (1 + alpha_n(Vsoma) )
	n_tau = max(50 * beta_n(Vsoma) / (1 + alpha_n(Vsoma)), 2.)
	∂n	         = (n_inf - n) / n_tau
	I_K		     = gamma_K * n * (Erev_K - Vsoma)

	##### NMDA
	NMDA         = (n1_nmda_A + n2_nmda_A + n1_nmda_B + n2_nmda_B) * B(Vsp, Mg) * gamma_nmda
	Inmda	     = (Erev_nmda - Vsp) * NMDA # current nmda

	##### AMPA
	Iampa	     = (Erev_ampa - Vsp) * (gamma_ampa1 * n1_ampa + gamma_ampa2 * n2_ampa + gamma_ampa3 * n3_ampa) # current ampa

	##### GABA
	Igaba        = (n_gaba1 + n_gaba2) *  (Erev_Cl - Vdend)  * gamma_GABA

	##### Calcium sources (VGCCs currents, and NMDA calcium contribution)
	ΦCa   = perm * ghk(Vsp, Ca, Ca_ext, p_synapse) #GHK factor
	Ica_nmda     = f_Ca * ΦCa * NMDA
	Icar	     = gamma_CaR * n_car * ΦCa
	Icat	     = gamma_CaT * n_cat * ΦCa
	Ical	     = gamma_CaL * n_cal * ΦCa

	##### SK channel (not stochastic)
	∂SK		     = (SK_chnnl(Ca) * frwd_SK - SK) / (SK_time * bcwd_SK) #SK spine
	Isk		     = SK_gamma * (SK_Erev - Vsp) * SK * N_SK

	##### Backpropgation
	# Post input - for experimentally induced BaPs and those induced by EPSPs
	I_BaP  = inputBaP(t, bap_by_epsp, injbap, I_clamp) + inputBaP(t, events_bap,  injbap, I_clamp)
	# Bap decay/attenuation - two component for adaptation in the Bap
	∂λ	     = (1-λ)/trec - delta_decay * (1/λ_aux) * λ * I_BaP
	∂λ_aux	 = (1-λ_aux)/trec - delta_aux * λ_aux * I_BaP
	gadapt	 = λ * g_diff * ϕ_dist

	# Bap decay/attenuation - age dependent modification factor
	∂λ_age	  = (1-λ_age)/trec_soma - delta_soma * λ_age * I_BaP

	##### Voltage
	# Spine
	∂Vsp	     = (Isk + Inmda + Iampa +Icat + Icar + Ical + g_neck * (Vdend - Vsp) + g_leak * (E_leak - Vsp)) / (Csp)
	# Dendrite
	∂Vdend       = (g_neck * (Vsp - Vdend) + Igaba + g_leakdend * (E_leak - Vdend) + gadapt * (Vsoma - Vdend)) / Cdend
	# Soma
	∂Vsoma       = ((I_BaP + I_Na)*λ_age + I_K + g_leaksoma * (E_leak - Vsoma) +  gadapt * (Vdend - Vsoma)) / Csoma


	##### Buffer and dye (spine only - no neck diffusion)
	∂ImbufCa     = Imbuf_k_on * (Imbuf_con - ImbufCa) * Ca - Imbuf_k_off * ImbufCa
	∂Dye	     = 4*fluo5f_kf * (fluo5f_con - Dye) * Ca - 8*fluo5f_kb * Dye

	##### Ca Downstream
	### CaM-KCaM-rates (coarsed model) from Pepke adapted by
	kf_2C	     = rates_adapt(kon_1C,   kon_2C,   koff_1C,  kon_2C, Ca)
	kb_2C	     = rates_adapt(koff_1C,  koff_2C,  koff_1C,  kon_2C, Ca)
	kf_2N	     = rates_adapt(kon_1N,   kon_2N,   koff_1N,  kon_2N, Ca)
	kb_2N	     = rates_adapt(koff_1N,  koff_2N,  koff_1N,  kon_2N, Ca)
	kf_K2C	     = rates_adapt(kon_K1C,  kon_K2C,  koff_K1C, kon_K2C, Ca)
	kb_K2C	     = rates_adapt(koff_K1C, koff_K2C, koff_K1C, kon_K2C, Ca)
	kf_K2N	     = rates_adapt(kon_K1N,  kon_K2N,  koff_K1N, kon_K2N, Ca)
	kb_K2N	     = rates_adapt(koff_K1N, koff_K2N, koff_K1N, kon_K2N, Ca)
	F            = CaMKII / mKCaM_con

	∂CaM0 = k2*PCaM0 + kb_2C*CaM2C + kb_2N*CaM2N + kb_CaM0*KCaM0 - (1//2)*kf_2C*(Ca^2)*CaM0 - (1//2)*kf_2N*(Ca^2)*CaM0 - kf_CaM0*CaM0*mKCaM
	∂CaM2C = kb_2N*CaM4 + kb_CaM2C*KCaM2C + k2*PCaM2C + (1//2)*kf_2C*(Ca^2)*CaM0 - kb_2C*CaM2C - (1//2)*kf_2N*(Ca^2)*CaM2C - kf_CaM2C*CaM2C*mKCaM
	∂CaM2N = kb_2C*CaM4 + kb_CaM2N*KCaM2N + k2*PCaM2N + (1//2)*kf_2N*(Ca^2)*CaM0 - kb_2N*CaM2N - (1//2)*kf_2C*(Ca^2)*CaM2N - kf_CaM2N*CaM2N*mKCaM
	∂CaM4 = k2*PCaM4 + kcanb*CaN4 + kb_CaM4*KCaM4 + (1//2)*kf_2N*(Ca^2)*CaM2C + (1//2)*kf_2C*(Ca^2)*CaM2N - kb_2C*CaM4 - kb_2N*CaM4 - kcanf*CaM4*mCaN - kf_CaM4*CaM4*mKCaM
	∂mCaN = kcanb*CaN4 - kcanf*CaM4*mCaN
	∂CaN4 = kcanf*CaM4*mCaN - kcanb*CaN4
	∂mKCaM = kb_CaM0*KCaM0 + k3*P + kb_CaM2C*KCaM2C + kb_CaM2N*KCaM2N + kb_CaM4*KCaM4 - kf_CaM0*CaM0*mKCaM - kf_CaM2C*CaM2C*mKCaM - kf_CaM2N*CaM2N*mKCaM - kf_CaM4*CaM4*mKCaM
	∂KCaM0 = kb_K2C*KCaM2C + kb_K2N*KCaM2N + kf_CaM0*CaM0*mKCaM - kb_CaM0*KCaM0 - (1//2)*kf_K2C*(Ca^2)*KCaM0 - F*k1*KCaM0 - (1//2)*kf_K2N*(Ca^2)*KCaM0
	∂KCaM2N = kb_K2C*KCaM4 + kf_CaM2N*CaM2N*mKCaM + (1//2)*kf_K2N*(Ca^2)*KCaM0 - kb_CaM2N*KCaM2N - kb_K2N*KCaM2N - (1//2)*kf_K2C*(Ca^2)*KCaM2N - F*k1*KCaM2N
	∂KCaM2C = kb_K2N*KCaM4 + kf_CaM2C*CaM2C*mKCaM + (1//2)*kf_K2C*(Ca^2)*KCaM0 - kb_CaM2C*KCaM2C - kb_K2C*KCaM2C - F*k1*KCaM2C - (1//2)*kf_K2N*(Ca^2)*KCaM2C
	∂KCaM4 = kf_CaM4*CaM4*mKCaM + (1//2)*kf_K2C*(Ca^2)*KCaM2N + (1//2)*kf_K2N*(Ca^2)*KCaM2C - kb_CaM4*KCaM4 - kb_K2C*KCaM4 - kb_K2N*KCaM4 - F*k1*KCaM4
	∂PCaM0 = F*k1*KCaM0 - k2*PCaM0
	∂PCaM2N = F*k1*KCaM2N - k2*PCaM2N
	∂PCaM2C = F*k1*KCaM2C - k2*PCaM2C
	∂PCaM4 = F*k1*KCaM4 - k2*PCaM4
	∂P = k2*PCaM0 + k5*P2 + k2*PCaM2C + k2*PCaM2N + k2*PCaM4 - k3*P - k4*P
	∂P2 = k4*P - k5*P2


	### Postsynaptic Ca
	∂Ca	         = (Ca_infty - Ca) / tau_ca +
					(Ica_nmda + Icar + Ical + Icat) / (2 * faraday * A_sp) +
					(max(Ca_infty, Ca/3) - Ca) / tau_diff -
					∂ImbufCa -
					∂Dye +
					2kb_2C*CaM2C + 2kb_2C*CaM4 + 2kb_2N*CaM2N + 2kb_2N*CaM4 + 2kb_K2C*KCaM2C + 2kb_K2N*KCaM2N + 2kb_K2C*KCaM4 + 2kb_K2N*KCaM4 - kf_2C*(Ca^2)*CaM0 - kf_2N*(Ca^2)*CaM0 - kf_2N*(Ca^2)*CaM2C - kf_2C*(Ca^2)*CaM2N - kf_K2C*(Ca^2)*KCaM0 - kf_K2C*(Ca^2)*KCaM2N - kf_K2N*(Ca^2)*KCaM0 - kf_K2N*(Ca^2)*KCaM2C

	### Xdot update
	xdot[1] =  ∂Vsp
	xdot[2] =  ∂Vdend
	xdot[3] =  ∂Vsoma
	xdot[4] =  ∂λ
	xdot[5] =  ∂ImbufCa
	xdot[6] =  ∂Ca
	xdot[7] =  ∂Dye
	xdot[8] =  ∂CaM0
	xdot[9] =  ∂CaM2C
	xdot[10] = ∂CaM2N
	xdot[11] = ∂CaM4
	xdot[12] = ∂mCaN
	xdot[13] = ∂CaN4
	xdot[14] = ∂mKCaM
	xdot[15] = ∂KCaM0
	xdot[16] = ∂KCaM2N
	xdot[17] = ∂KCaM2C
	xdot[18] = ∂KCaM4
	xdot[19] = ∂PCaM0
	xdot[20] = ∂PCaM2C
	xdot[21] = ∂PCaM2N
	xdot[22] = ∂PCaM4
	xdot[23] = ∂P
	xdot[24] = ∂P2
	xdot[25] = ∂LTD
	xdot[26] = ∂LTP
	xdot[27] = ∂act_D
	xdot[28] = ∂act_P
	xdot[29] = ∂m
	xdot[30] = ∂h
	xdot[31] = ∂n
	xdot[32] = ∂SK
	xdot[33] = ∂λ_age
	xdot[34] = ∂λ_aux
end

function R_synapse(rate, xc, xd, p_synapse::SynapseParams, t, sum_rate, glu = 0)

	@unpack_SynapseParams p_synapse

	############### Voltage ###################
	Vsp = xc[1]

	############### Glutamate & GABA ###################
	Glu = glu_amp * glu

	############### AMPA ###################
	#2line-GO
	rate[1]  = 4 * AMPA_k1 * Glu * xd[1]
	rate[2]  = 3 * AMPA_k1 * Glu * xd[2]
	rate[3]  = 2 * AMPA_k1 * Glu * xd[3]
	rate[4]  = 1 * AMPA_k1 * Glu * xd[4]
	#2line-BACK
	rate[5]  = 4 * AMPA_k_1 * xd[5]
	rate[6]  = 3 * AMPA_k_1 * xd[4]
	rate[7]  = 2 * AMPA_k_1 * xd[3]
	rate[8]  = 1 * AMPA_k_1 * xd[2]
	#3line-GO
	rate[9]   = 3 * AMPA_k1 * Glu * xd[6]
	rate[10]  = 3 * AMPA_k1 * Glu * xd[7]
	rate[11]  = 2 * AMPA_k1 * Glu * xd[8]
	rate[12]  = 1 * AMPA_k1 * Glu * xd[9]
	#3line-BACK
	rate[13]  = 3 * AMPA_k_1 * xd[10]
	rate[14]  = 2 * AMPA_k_1 * xd[9]
	rate[15]  = 1 * AMPA_k_1 * xd[8]
	rate[16]  = 1 * AMPA_k_2 * xd[7]
	#4line-GO
	rate[17]  = 2 * AMPA_k1 * Glu * xd[11]
	rate[18]  = 1 * AMPA_k1 * Glu * xd[12]
	#4line-BACK
	rate[19]  = 2 * AMPA_k_1 * xd[13]
	rate[20]  = 1 * AMPA_k_1 * xd[12]
	#1column-GO-BACK
	rate[21]  = 4 * AMPA_delta_0 * xd[1]
	rate[22]  = 1 * AMPA_gamma_0 * xd[6]
	#2column-GO-BACK
	rate[23]  = 1 * AMPA_delta_1 * xd[2]
	rate[24]  = 1 * AMPA_gamma_1 * xd[7]
	#3column-GO
	rate[25]  = 1 * AMPA_alpha * xd[14]
	rate[26]  = 2 * AMPA_delta_1 * xd[3]
	rate[27]  = 1 * AMPA_delta_2 * xd[8]
	#3column-BACK
	rate[28]  = 1 * AMPA_gamma_2 * xd[11]
	rate[29]  = 1 * AMPA_gamma_1 * xd[8]
	rate[30]  = 2 * AMPA_beta * xd[3]
	#4column-GO
	rate[31]  = 1 * AMPA_alpha * xd[15]
	rate[32]  = 3 * AMPA_delta_1 * xd[4]
	rate[33]  = 2 * AMPA_delta_2 * xd[9]
	#4column-BACK
	rate[34]  = 1 * AMPA_gamma_2 * xd[12]
	rate[35]  = 1 * AMPA_gamma_1 * xd[9]
	rate[36]  = 2 * AMPA_beta * xd[4]
	#5column-GO
	rate[37]  = 1 * AMPA_alpha * xd[16]
	rate[38]  = 4 * AMPA_delta_1 * xd[5]
	rate[39]  = 3 * AMPA_delta_2 * xd[10]
	#5column-BACK
	rate[40]  = 1 * AMPA_gamma_2 * xd[13]
	rate[41]  = 1 * AMPA_gamma_1 * xd[10]
	rate[42]  = 4 * AMPA_beta * xd[5]

	############### NMDA ###################
	#1line-GO
	rate[43]  =   NMDA_N2A_ka * xd[17] * Glu
	rate[44]  =   NMDA_N2A_kb * xd[18] * Glu
	rate[45]  =   NMDA_N2A_kc * xd[19]
	rate[46]  =   NMDA_N2A_kd * xd[20]
	rate[47]  =   NMDA_N2A_ke * xd[21]
	rate[48]  =   NMDA_N2A_kf * xd[22]
	#1line-BACK

	rate[49]  =   NMDA_N2A_k_f * xd[23]
	rate[50]  =   NMDA_N2A_k_e * xd[22]
	rate[51]  =   NMDA_N2A_k_d * xd[21]
	rate[52]  =   NMDA_N2A_k_c * xd[20]
	rate[53]  =   NMDA_N2A_k_b * xd[19]
	rate[54]  =   NMDA_N2A_k_a * xd[18]

	################### Sampling ###################
	rate[55]  = sampling_rate

	################### R-type VGCC ###################
	alpha_m_r, beta_m_r = rates_m_r(Vsp)
	alpha_h_r, beta_h_r = rates_h_r(Vsp)

	rate[56] = xd[25] * alpha_m_r * frwd_VGCC
	rate[57] = xd[26] * beta_m_r  * bcwd_VGCC
	rate[58] = xd[25] * alpha_h_r * frwd_VGCC
	rate[59] = xd[27] * beta_h_r  * bcwd_VGCC
	rate[60] = xd[26] * alpha_h_r * frwd_VGCC
	rate[61] = xd[28] * beta_h_r  * bcwd_VGCC
	rate[62] = xd[27] * alpha_m_r * frwd_VGCC
	rate[63] = xd[28] * beta_m_r  * bcwd_VGCC

	################### T-type VGCC  ###################
	alpha_m_t, beta_m_t = rates_m_t(Vsp)
	alpha_h_t, beta_h_t = rates_h_t(Vsp)

	rate[64] = xd[29] * alpha_m_t * frwd_VGCC
	rate[65] = xd[30] * beta_m_t  * bcwd_VGCC # this one can have a high rate
	rate[66] = xd[29] * alpha_h_t * frwd_VGCC
	rate[67] = xd[31] * beta_h_t  * bcwd_VGCC
	rate[68] = xd[30] * alpha_h_t * frwd_VGCC
	rate[69] = xd[32] * beta_h_t  * bcwd_VGCC
	rate[70] = xd[31] * alpha_m_t * frwd_VGCC
	rate[71] = xd[32] * beta_m_t  * bcwd_VGCC # this one can have a high rate

	################### L-type VGCC  ###################
	alpha_l, beta_1_l, beta_2_l = rates_l(Vsp)
	rate[72] = xd[33] * alpha_l  * frwd_VGCC
	rate[73] = xd[34] * beta_1_l * bcwd_VGCC
	rate[74] = xd[33] * alpha_l  * frwd_VGCC
	rate[75] = xd[35] * beta_2_l * bcwd_VGCC

	################### LTD/LTP  ###################
	#the 6 lines take 50ns on 200ns, 1/4 of computations are here!!
	D_rate = plasticityRate(xc[27], 2, K_D) / t_D
	P_rate = plasticityRate(xc[28], 2, K_P) / t_P
	rate[76] = xd[36] * D_rate
	rate[77] = xd[37] * P_rate
	rate[78] = xd[36] * P_rate
	rate[79] = xd[38] * D_rate

	############### NMDA GLUN2B ###################
	#1line-GO
	rate[80]  =   NMDA_N2B_sa * xd[39] * Glu
	rate[81]  =   NMDA_N2B_sb * xd[40] * Glu
	rate[82]  =   NMDA_N2B_sc * xd[41]
	rate[83]  =   NMDA_N2B_sd * xd[42]
	rate[84]  =   NMDA_N2B_se * xd[43]
	rate[85]  =   NMDA_N2B_sf * xd[44]
	#1line-BACK

	rate[86]  =   NMDA_N2B_s_f * xd[45]
	rate[87]  =   NMDA_N2B_s_e * xd[44]
	rate[88]  =   NMDA_N2B_s_d * xd[43]
	rate[89]  =   NMDA_N2B_s_c * xd[42]
	rate[90]  =   NMDA_N2B_s_b * xd[41]
	rate[91]  =   NMDA_N2B_s_a * xd[40]

	############### GABA ###################

	rate[92]  =  GABA_r_b1  * xd[46] * Glu #to simplify, we use the same ammount at the same time
	rate[93]  =  GABA_r_u1  * xd[47]
	rate[94]  =  GABA_r_b2  * xd[47] * Glu
	rate[95]  =  GABA_r_u2  * xd[48]
	rate[96]  =  GABA_r_ro1 * xd[47]
	rate[97]  =  GABA_r_c1  * xd[49]
	rate[98]  =  GABA_r_ro2 * xd[48]
	rate[99]  =  GABA_r_c2  * xd[50]

	bound = 0.
	if sum_rate == false
		return 0., bound
	else
		return sum(rate), bound
	end
end

function buildTransitionMatrix()
	matrix_list = [AMPA_matrix()]
	push!(matrix_list, NMDA_matrix()) #for GluN2A
	push!(matrix_list, Matrix{Int64}(I, 1, 1)) #Print from Poisson Rate
	push!(matrix_list, R_channel_matrix())
	push!(matrix_list, T_channel_matrix())
	push!(matrix_list, L_channel_matrix())
	push!(matrix_list, LTP_LTD_matrix())
	push!(matrix_list, NMDA_matrix())  #for GluN2B
	push!(matrix_list, GABA_matrix())
	return sparse(jump_matrix(matrix_list))
end

"""
$(SIGNATURES)


"""
function pdmpsynapse(xc, xd, t1, t2, events_bap, bap_by_epsp, glu, p_synapse, nu; algo = CHV(:lsoda), kwargs...)
	problem = PDMP.PDMPProblem(
		(xdot, xc, xd, p, t) -> F_synapse(xdot, xc, xd, p, t, events_bap, bap_by_epsp),
		(rate, xc, xd, p, t, sum_rate) -> R_synapse(rate, xc, xd, p, t, sum_rate, glu),
		nu, xc, xd, p_synapse, (t1, t2);
		Ncache = 12) # this option is for AD in PreallocationTools
	return solve(problem, algo; kwargs...)
end

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
		events_sorted_times,
		is_pre_or_post_event,
		bap_by_epsp,
		is_glu_released,
		algos;
		verbose = false, progress = false,
		abstol = 1e-8, reltol = 1e-7,
		save_positions = (false, true),
		nu = buildTransitionMatrix(), kwargs...) where T

		tt, XC, XD = evolveSynapse_noformat(xc0, xd0, p_synapse,
						events_sorted_times,
						is_pre_or_post_event,
						bap_by_epsp,
						is_glu_released,
						algos;
				verbose = verbose, progress = progress,
				abstol = abstol, reltol = reltol, save_positions = save_positions,
				nu = nu, kwargs...)

		# format the output to make it convenient to parse
		# this is wasting a lot of ressources but is convenient for plotting
		verbose && @printf("=> done! parsing results")
		out = formatSynapseResult(tt, XC, XD)
end


"""
Same as `evolveSynapse` but do not format the output because it takes time.
"""
function evolveSynapse_noformat(xc0::Vector{T}, xd0, p_synapse::SynapseParams,
		events_sorted_times,
		is_pre_or_post_event,
		bap_by_epsp,
		is_glu_released,
		algos;
		verbose = false, progress = false,
		abstol = 1e-8, reltol = 1e-7, save_positions = (false, true),
		nu = buildTransitionMatrix(), kwargs...) where T

	if save_positions isa Tuple{Bool, Bool}
		save_positionsON = save_positions
		save_positionsOFF = save_positions
		# save_positionsLAST = save_positions
	else
		save_positionsON = save_positions[1]
		save_positionsOFF = save_positions[2]
		# save_positionsLAST = save_positions
	end

	@assert eltype(is_pre_or_post_event) == Bool "Provide booleans for glutamate releases."
	@assert eltype(is_glu_released) == Bool "Provide booleans for glutamate indices."

	verbose && printstyled(color=:red,"\n+++++++++++++++++++++++++++++\n")
	verbose && printstyled(color=:red,"Synapse simulation")

	XC = VectorOfArray([xc0]) # vector to hold continuous variables
	XD = VectorOfArray([xd0]) # vector to hold discrete variables
	tt = [0.0]				  # vector of times

	# results from the simulation of an external event
	res = PDMP.PDMPResult([0., 0.], copy(XC), copy(XD))

	# we collect which external events correspond to BaPs
	events_bap = events_sorted_times[is_pre_or_post_event .== false]

	# function to simulate the synapse when Glutamate is ON
	SimGluON = (xc, xd, t1, t2, glu) -> pdmpsynapse(xc, xd, t1, t2, events_bap, bap_by_epsp, glu, p_synapse, nu; algo = algos[1], save_positions = save_positionsON, reltol = reltol, abstol = abstol, kwargs...)

	# function to simulate the synapse when Glutamate is OFF
	SimGluOFF = (xc, xd, t1, t2)	 -> pdmpsynapse(xc, xd, t1, t2, events_bap, bap_by_epsp, zero(T), p_synapse, nu;  algo = algos[2], save_positions = save_positionsOFF, reltol = reltol, abstol = abstol, kwargs...)

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
			res = SimGluOFF(res.xc[:,end], res.xd[:,end], tt[end], eve)
			append!(XC, res.xc);  append!(XD, res.xd);  append!(tt, res.time)
			gluamp = rand(gluDist)
			verbose && @printf("=> Glu on, %4d, t ∈ [%9.4e, %9.4e]\n", eveindex, eve, eve+ p_synapse.glu_width )

			# simulate the event with Glutamate ON
			# variability here
			res = SimGluON(res.xc[:,end], res.xd[:,end], eve, eve + p_synapse.glu_width,  ifelse(is_glu_released[eveindex], gluamp, zero(T)))
			append!(XC,res.xc);  append!(XD,res.xd);  append!(tt,res.time)
		end
		# update the progress bar
		progress && next!(pbar; showvalues = [(:steps, eveindex), (:t, tt[end])])
	end

	# reaching tend: we simulate the synapse with Glutamate OFF until simulation end time required
	# by the user. In  most protocol, this is taking most of the time.
	verbose && @printf("=> Reaching the end, t ∈ [%9.4e, %9.4e]\n",tt[end], p_synapse.t_end)
	res = @time SimGluOFF(res.xc[:,end], res.xd[:,end], tt[end], p_synapse.t_end)
	@info "last bit" length(res.time) tt[end] p_synapse.t_end
	append!(XC, res.xc);  append!(XD, res.xd);  append!(tt, res.time)

	# update the progress bar
	progress && next!(pbar; showvalues = [(:steps, length(events_sorted_times) + 1), (:t, p_synapse.t_end)])

	@assert res.time[end] == p_synapse.t_end "Error in PDMP. Did not reach requested simulated time","complete simulated time"

	return (t = tt, XC = XC, XD = XD)
end

function formatSynapseResult(tt, XC, XD)
	namesC = (:Vsp, :Vdend, :Vsoma, :λ, :ImbufCa, :Ca, :Dye,
				 :CaM0, :CaM2C, :CaM2N, :CaM4,
				 :mCaN, :CaN4, :mKCaM,
				 :KCaM0, :KCaM2N, :KCaM2C, :KCaM4,
				 :PCaM0, :PCaM2C, :PCaM2N, :PCaM4,
				 :P, :P2, :LTD, :LTP, :act_D, :act_P,
				 :m, :h, :n, :SK ,:λ_age, :λ_aux)
	values = (XC[i, :] for i in 1:length(namesC))
	return (t = tt, XD = XD, XC = XC, zip(namesC, values)...)
end

function indexOfVariable(name::Symbol)
	names = (:Vsp, :Vdend, :Vsoma, :λ, :ImbufCa, :Ca, :Dye,
				 :CaM0, :CaM2C, :CaM2N, :CaM4,
				 :mCaN, :CaN4, :mKCaM,
				 :KCaM0, :KCaM2N, :KCaM2C, :KCaM4,
				 :PCaM0, :PCaM2C, :PCaM2N, :PCaM4,
				 :P, :P2, :LTD, :LTP, :act_D, :act_P,
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
