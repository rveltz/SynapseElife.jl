# using Revise, Plots
using Test
using Synapse
using PiecewiseDeterministicMarkovProcesses, JumpProcesses, OrdinaryDiffEq, Sundials

# this holds the spine parameters
param_synapse = SynapseParams(t_end = 1000.)

# initial conditions for the spine
xc0 = initial_conditions_continuous_temp(param_synapse) # initial conditions deterministic vars
xd0 = initial_conditions_discrete(param_synapse) # initial conditions stochastic channels

# we put a presynaptic pulse at 500ms and simulate for 1s
events_times = [500.]
is_pre_or_post_event = [true]

# ODE time stepper
ode = :lsoda
ode = CVODE_BDF()
ode = Tsit5()

# Discrete aggregator
agg = Direct()

# simulate the spine dynamics
result = @time evolveSynapse(
		xc0,
		xd0,
		param_synapse,
		events_times,			# external events
		is_pre_or_post_event,	# pre or post?
		Float64[],				# auxiliary BaP
		[true],
		# (CHV(ode), CHV(ode));
		(ode, ode);
		save_positions = (false, true),
		agg = agg)

@test ~isnothing(result)

# # plot the discrete variables
# Synapse.plot_discrete(result.t, result.XC, result.XD)
#
# # plot specific variable
# Synapse.plot_variable(result.t, result.XC, result.XD, :Vsp; xlim = (480., 520))
