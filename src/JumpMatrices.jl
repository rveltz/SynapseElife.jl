######### jumps matrix
function stp_build_transition_matrix()
	matrix_list = [Matrix{Int64}(I,1,1), STP_matrix()]
	return jump_matrix(matrix_list)
end

# Catalyst uses transposed transition matrix since version > 8.3
get_stoichmatrix(model) = transpose(Catalyst.netstoichmat(model)) |> Matrix

jump_matrix(matrix_list) = cat(matrix_list..., dims=(1,2))

ampa_model = @reaction_network begin
	#1line
	#2line-GO
	1,	C0 → C1
	1,	C1 → C2
	1,	C2 → C3
	1,	C3 → C4
	#2line-BAC
	1,	C4 → C3
	1,	C3 → C2
	1,	C2 → C1
	1,	C1 → C0
	#3line-GO
	1,	D0 → D1
	1,	D1 → D2
	1,	D2 → D3
	1,	D3 → D4
	#3line-BACK
	1,	D4 → D3
	1,	D3 → D2
	1,	D2 → D1
	1,	D1 → D0
	#4line-GO
	1,	D22 → D23
	1,	D23 → D24
	#4line-BACK
	1,	D24 → D23
	1,	D23 → D22
	#1column-GO-BACK
	1,	C0 → D0
	1,	D0 → C0
	#2column-GO-BACK
	1,	C1 → D1
	1,	D1 → C1
	#3column-GO
	1,	O2 → C2
	1,	C2 → D2
	1,	D2 → D22
	#3column-BACK
	1,	D22 → D2
	1,	D2 → C2
	1,	C2 → O2
	#4column-GO
	1,	O3 → C3
	1,	C3 → D3
	1,	D3 → D23
	#4column-BACK
	1,	D23 → D3
	1,	D3 → C3
	1,	C3 → O3
	#5column-GO
	1,	O4 → C4
	1,	C4 → D4
	1,	D4 → D24
	#5column-BACK
	1,	D24 → D4
	1,	D4 → C4
	1,	C4 → O4
end

plasticity = @reaction_network begin
	1,	NC → LTD
	1,	LTD → NC
	1,	NC → LTP
	1,	LTP → NC
end

nmda_model_v2 = @reaction_network begin #same structure for N2A and N2B
	#1line-GO
	1,	A0 → A1
	1,	A1 → A2
	1,	A2 → A3
	1,	A3 → A4
	1,	A4 → AO1
	1,	AO1 → AO2
	#2line-BACK
	1,	AO2 → AO1
	1,	AO1 → A4
	1,	A4 → A3
	1,	A3 → A2
	1,	A2 → A1
	1,	A1 → A0
end

gaba_destexhe = @reaction_network begin
	#1line-GO
	(1,1),	CO ↔ C1
	(1,1),	C1 ↔ C2
	(1,1),	C1 ↔ O1
	(1,1),	C2 ↔ O2
end

#SK channel not used due to tonic flickering
sk_model = @reaction_network begin
	#1line-GO
	1,	C → O
	1,	O → C
end

AMPA_matrix()       = get_stoichmatrix(ampa_model)
LTP_LTD_matrix()    = get_stoichmatrix(plasticity)
NMDA_matrix()       = get_stoichmatrix(nmda_model_v2)
SK_channel_matrix() = get_stoichmatrix(sk_model)
GABA_matrix()       = get_stoichmatrix(gaba_destexhe)

STP_matrix() = [
	[ -1  1]; # add vesicle to RRP
	[  1 -1]; # move vesicle from recovery pool to RRP
	[  0  1]; # add veiscle to recovery pool
]

L_channel_matrix() = [
	[-1   1  0]; # CaL C -> O1
	[ 1  -1  0]; # CaL O1 -> C
	[-1   0  1]; # CaL C -> O2
	[ 1   0 -1]; # CaL O2 -> C
]

R_channel_matrix() = [  
	[ -1  1  0  0]; # CaR m0h0 -> m1h0
	[  1 -1  0  0]; # CaR m1h0 -> m0h0
	[ -1  0  1  0]; # CaR m0h0 -> m0h1
	[  1  0 -1  0]; # CaR m0h1 -> m0h0
	[  0 -1  0  1]; # CaR m1h0 -> O
	[  0  1  0 -1]; # CaR    O -> m1h0
	[  0  0 -1  1]; # CaR m0h1 -> O
	[  0  0  1 -1]; # CaR    O -> m0h1
]

T_channel_matrix() = [ 
	[ -1  1  0  0]; # CaT m0h0 -> m1h0
	[  1 -1  0  0]; # CaT m1h0 -> m0h0
	[ -1  0  1  0]; # CaT m0h0 -> m0h1
	[  1  0 -1  0]; # CaT m0h1 -> m0h0
	[  0 -1  0  1]; # CaT m1h0 -> O
	[  0  1  0 -1]; # CaT    O -> m1h0
	[  0  0 -1  1]; # CaT m0h1 -> O
	[  0  0  1 -1]; # CaT    O -> m0h
]

function buildTransitionMatrix()
	matrix_list = [AMPA_matrix()]
	push!(matrix_list, NMDA_matrix())          # for GluN2A
	push!(matrix_list, Matrix{Int64}(I, 1, 1)) # print from Poisson Rate
	push!(matrix_list, R_channel_matrix())
	push!(matrix_list, T_channel_matrix())
	push!(matrix_list, L_channel_matrix())
	push!(matrix_list, LTP_LTD_matrix())
	push!(matrix_list, NMDA_matrix())          # for GluN2B
	push!(matrix_list, GABA_matrix())
	return sparse(jump_matrix(matrix_list))
end
