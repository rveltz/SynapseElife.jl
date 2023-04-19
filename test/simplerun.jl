# using Revise, Plots
using Test
using Synapse
using PiecewiseDeterministicMarkovProcesses, JumpProcesses, OrdinaryDiffEq, Sundials, LSODA

# this holds the spine parameters
p_synapse = SynapseParams(t_end = 1000.)

# initial conditions for the spine
xc0 = initial_conditions_continuous_temp(p_synapse) # initial conditions deterministic vars
xd0 = initial_conditions_discrete(p_synapse) # initial conditions stochastic channels

# we put a presynaptic pulse at 500ms and simulate for 1s
events_times = [500.]
is_pre_or_post_event = [true]

# transition matrix
nu = buildTransitionMatrix()

# ODE time stepper
# algos = (CHV(:lsoda), CHV(:lsoda))
# algos = (CHV(CVODE_BDF()), CHV(CVODE_BDF()))
# algos = (Tsit5(), Tsit5())
algos = (TRBDF2(), TRBDF2())
# algos = (lsoda(), lsoda())
# algos = (AutoTsit5(Rosenbrock23()), AutoTsit5(Rosenbrock23()))
# algos = (CVODE_BDF(), CVODE_BDF())

# Discrete aggregator
# agg = nothing
agg = Direct()
# agg = CoevolveSynced()

# simulate the spine dynamics
result = @time evolveSynapse(
	xc0,
	xd0,
	p_synapse,
	events_times,           # external events
	is_pre_or_post_event,   # pre or post?
	Float64[],
	[true],
	nu,
	algos,
	agg;
	save_positions = (false, true),
)

@test ~isnothing(result)

# # plot the discrete variables
# Synapse.plot_discrete(result.t, result.XC, result.XD)
#
# # plot specific variable
Synapse.plot_variable(result.t, result.XC, result.XD, :Vsp; xlim = (480., 520))
