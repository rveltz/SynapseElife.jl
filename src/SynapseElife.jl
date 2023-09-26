module SynapseElife
	using Parameters, Distributions, SparseArrays, DataFrames, ProgressMeter

	using DocStringExtensions

	using Catalyst
	using Printf
	using LazySets # for using ∈ plasticity region
	using StaticArrays # for using LazySets in efficient way
	using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, RecursiveArrayTools, Plots
	const PDMP = PiecewiseDeterministicMarkovProcesses

	include("ParamsSynapse.jl")
	include("UtilsData.jl")
	include("UtilsDynamics.jl")
	include("JumpMatrices.jl")
	include("SynapseModel.jl")
	include("OnlyStp.jl")
	include("CaM-KCaM-reactions.jl")

	export PreSynapseParams, SynapseParams, SynapseParamsDet,  firingPattern, initial_conditions_continuous_temp, initial_conditions_discrete, initial_conditions_deterministic

	export F_synapse, R_synapse, F_synapse_ds, R_synapse_ds

	export writeEquations

	export dataProtocol, buildTransitionMatrix, buildTransitionMatrix_ds

	export stp, evolveSynapse_ds, evolveSynapse_noformat_ds, evolveSynapse, evolveSynapse_noformat
	export getCaM, getCamKII, getCaN
end
