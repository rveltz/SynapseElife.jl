module Synapse
	using Parameters, Distributions, Documenter, SparseArrays, DataFrames, ProgressMeter

	using DocStringExtensions

	using Catalyst
	using Printf
	using LazySets # for using ∈ plasticity region
	using StaticArrays # for using LazySets in efficient way
	using Plots
	using LinearAlgebra, RecursiveArrayTools
	using DifferentialEquations, LSODA, Sundials
	using JumpProcesses
	using PiecewiseDeterministicMarkovProcesses
	const PDMP = PiecewiseDeterministicMarkovProcesses

	include("ParamsSynapse.jl")
	include("UtilsData.jl")
	include("UtilsDynamics.jl")
	include("JumpMatrices.jl")
	include("ModelPDMP.jl")
	include("ModelODE.jl")
	include("Solve.jl")
	include("OnlyStp.jl")
	include("CaM-KCaM-reactions.jl")

	export PreSynapseParams, SynapseParams, SynapseParamsDet,  firingPattern, initial_conditions_continuous_temp, initial_conditions_discrete, initial_conditions_deterministic

	export F_synapse, R_synapse, F_synapse_ds, R_synapse_ds

	export G_synapse, J_synapse, model_jump, buildRxDependencyGraph

	export writeEquations

	export dataProtocol, buildTransitionMatrix, buildTransitionMatrix_ds

	export stp, evolveSynapse_ds, evolveSynapse_noformat_ds, evolveSynapse, evolveSynapse_noformat
	export getCaM, getCamKII, getCaN
end
