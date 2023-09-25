# Tutorial

We show here a very basic simulation where we apply a glutamate pulse to the postsynaptic side.

```@example TUT
using Plots, SynapseElife
using Sundials, LSODA # ODE solvers
# we load this package to simulate this particular class of Markov processes
using PiecewiseDeterministicMarkovProcesses

# this holds the spine parameters
# we change the simulation time parameter
param_synapse = SynapseParams(t_end = 1000.)

# initial conditions for the spine
# initial conditions deterministic vars
xc0 = initial_conditions_continuous_temp(param_synapse)
# initial conditions stochastic channels
xd0 = initial_conditions_discrete(param_synapse)

# we put a presynaptic pulse at 500ms, actually we put a glutamate pulse.
events_times = [500.]
# we say that this is a pre-synaptic simulation. 
# If we put `false`, it would be a post-synaptic stimulus
is_pre_or_post_event = [true]

# ODE time stepper
ode = :lsoda
ode = CVODE_BDF()

# simulate the spine dynamics
result = @time evolveSynapse(
		xc0,
		xd0,
		param_synapse,
		events_times,			# external events
		is_pre_or_post_event,	# pre or post?
		Float64[],			# auxiliary BaP, empty list
		[true],				# whether the pre-synaptic event leads to Glu release
		(CHV(ode), CHV(ode));
		save_positions = (false, true) # we want to record only the post-jump values
		)
nothing #hide
```

We can for example plot all the discrete variables

```@example TUT
# plot the discrete variables
SynapseElife.plot_discrete(result.t, result.XC, result.XD)
```

We can also select a specific variable to be printed:

```@example TUT
# plot specific variable
SynapseElife.plot_variable(result.t, result.XC, result.XD, :Vsp; xlim = (480., 520))
```
