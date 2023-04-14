
function construct_ma_p(p_synapse::SynapseParams, nu)
	@unpack_SynapseParams p_synapse
	p = (
		#2line-GO
		4 * AMPA_k1 * Glu, # 1
		3 * AMPA_k1 * Glu, # 2
		2 * AMPA_k1 * Glu, # 3
		1 * AMPA_k1 * Glu, # 4
		#2line-BACK
		4 * AMPA_k_1, # 5
		3 * AMPA_k_1, # 6
		2 * AMPA_k_1, # 7
		1 * AMPA_k_1, # 8
		#3line-GO
		3 * AMPA_k1 * Glu, # 9
		3 * AMPA_k1 * Glu, # 10
		2 * AMPA_k1 * Glu, # 11
		1 * AMPA_k1 * Glu, # 12
		#3line-BACK
		3 * AMPA_k_1, # 13
		2 * AMPA_k_1, # 14
		1 * AMPA_k_1, # 15
		1 * AMPA_k_2, # 16
		#4line-GO
		2 * AMPA_k1 * Glu, # 17
		1 * AMPA_k1 * Glu, # 18
		#4line-BACK
		2 * AMPA_k_1, # 19
		1 * AMPA_k_1, # 20
		#1column-GO-BACK
		4 * AMPA_delta_0, # 21
		1 * AMPA_gamma_0, # 22
		#2column-GO-BACK
		1 * AMPA_delta_1, # 23
		1 * AMPA_gamma_1, # 24
		#3column-GO
		1 * AMPA_alpha, # 25
		2 * AMPA_delta_1, # 26
		1 * AMPA_delta_2, # 27
		#3column-BACK
		1 * AMPA_gamma_2, # 28
		1 * AMPA_gamma_1, # 29
		2 * AMPA_beta, # 30
		#4column-GO
		1 * AMPA_alpha, # 31
		3 * AMPA_delta_1, # 32
		2 * AMPA_delta_2, # 33
		#4column-BACK
		1 * AMPA_gamma_2, # 34
		1 * AMPA_gamma_1, # 35
		2 * AMPA_beta, # 36
		#5column-GO
		1 * AMPA_alpha, # 37
		4 * AMPA_delta_1, # 38
		3 * AMPA_delta_2, # 39
		#5column-BACK
		1 * AMPA_gamma_2, # 40
		1 * AMPA_gamma_1, # 41
		4 * AMPA_beta, # 42

		############### NMDA ###################
		#1line-GO
		NMDA_N2A_ka * Glu, # 43
		NMDA_N2A_kb * Glu, # 44
		NMDA_N2A_kc, # 45
		NMDA_N2A_kd, # 46
		NMDA_N2A_ke, # 47
		NMDA_N2A_kf, # 48
		#1line-BACK

		NMDA_N2A_k_f, # 49
		NMDA_N2A_k_e, # 50
		NMDA_N2A_k_d, # 51
		NMDA_N2A_k_c, # 52
		NMDA_N2A_k_b, # 53
		NMDA_N2A_k_a, # 54

		############# Original pararameters ############
		p_synapse,

	)
end

function ma_jumps(p_synapse::SynapseParams, nu)

	p = construct_ma_p(p_synapse, nu)

	reactions, species = size(nu)

	reactant_stoich = Vector{Vector{Pair{Int,Int}}}(undef, reactions)
	reactants = (nu .< 0)
	for i in 1:reactions
		js = findall(reactants[i, :])
		stoich = [j => -nu[i, j] for j in js] 
		reactant_stoich[i] = stoich
	end

	nets = (nu .> 0)
	net_stoich = Vector{Vector{Pair{Int,Int}}}(undef, reactions)
	for i in 1:reactions
		js = findall(nets[i, :])
		stoich = [j => nu[i, j] for j in js] 
		net_stoich[i] = stoich
	end

	param_idxs = 1:reactions

	jump = MassActionJump(reactant_stoich, net_stoich; scale_rates = false, param_idxs)

	return p, jump
end

function buildMATransitionMatrix()
	matrix_list = [AMPA_matrix()]
	push!(matrix_list, NMDA_matrix()) #for GluN2A
	return sparse(jump_matrix(matrix_list))
end

function R_jumps()
	nu = R_channel_matrix()

	stoich = nu[1, :] .!= 0 
	function rate56(u, p, t)
		@unpack_SynapseParams p[end]
		alpha_m_r, beta_m_r = rates_m_r(Vsp)
		Vsp = u[49]
		return u[25] * alpha_m_r * frwd_VGCC
	end
	function urate56(u, p, t)
	end
	function rateinterval56(u, p, t)
	end
	function affect56!(integrator)
		for i in stoich
			integrator.u[24+i] += nu[1, i]
		end
	end
	jump56 = VariableRateJump(rate56, affect56!; urate=urate56, rateinterval=rateinterval56)

end


# function F_synapse()
# end

# function R_synapse()
# 	jumps = [
# 		@synapse_rate 4 * AMPA_k1 * Glu * u[1]
# 	]
# end

# function rate1(u, p::SynapseParams, t)
# 	@unpack_SynapseParams p
# 	Glu = glu_amp * glu
# 	return 4 * AMPA_k1 * Glu * u[1]
# end

# function affect1!(u, p::SynapseParams, t)
# end

# function rate2(u, p::SynapseParams, t)
# 	@unpack_SynapseParams p
# 	Glu = glu_amp * glu
# 	return 3 * AMPA_k1 * Glu * u[2]
# end

# function affect2!(u, p::SynapseParams, t)
# end

# rate_factory((u, p, t) -> 4 * AMPA_k1 * Glu * u[1])

# macro simple_rate(i, nu, rate)
# 	return quote
# 		_, js, nu_ij = findnz(nu[i, :]) 
# 		function rate(u, p, t)
# 			@unpack_SynapseParams p
# 			Glu = glu_amp * glu
# 			return $ex
# 		end
# 		function affect(u, p, t)
# 			for j, k

# 		end

# 	end
# end
#


# function pdmpsynapse()

# end
