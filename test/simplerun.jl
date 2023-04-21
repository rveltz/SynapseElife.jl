# using Revise, Plots
using Test
using Synapse
using PiecewiseDeterministicMarkovProcesses, JumpProcesses, OrdinaryDiffEq, Sundials, LSODA
const PDMP = PiecewiseDeterministicMarkovProcesses

# this holds the spine parameters
p_synapse = SynapseParams(t_end = 1000.)

# initial conditions for the spine
xc0 = initial_conditions_continuous_temp(p_synapse) # initial conditions deterministic vars
xd0 = initial_conditions_discrete(p_synapse) # initial conditions stochastic channels

# we put a presynaptic pulse at 500ms and simulate for 1s
events_sorted_times = [500.]
is_pre_or_post_event = [true]

# transition matrix
nu = buildTransitionMatrix()

# ODE time stepper
algos = (CHV(:lsoda), CHV(:lsoda))
# algos = (CHV(CVODE_BDF()), CHV(CVODE_BDF()))
# algos = (Tsit5(), Tsit5())
# algos = (TRBDF2(), TRBDF2())
# algos = (lsoda(), lsoda())
algos = (AutoTsit5(Rosenbrock23()), AutoTsit5(Rosenbrock23()))
# algos = (CVODE_BDF(), CVODE_BDF())

# Discrete aggregator
# agg = nothing
# agg = Direct()
agg = CoevolveSynced()

# simulate a single problem
glu = 0.0
events_bap = events_sorted_times[is_pre_or_post_event .== false]
bap_by_epsp = Float64[]
# u = ExtendedJumpArray(xc0, xd0)
u = vcat(xc0, xd0)
p, jumps = J_synapse(p_synapse, nu, glu)
t1 = 0.
t2 = 500.
save_positions = (false, false)
oprob = ODEProblem((du, u, p, t) -> G_synapse(du, u, p_synapse, t, events_bap, bap_by_epsp), u, 
	(t1, t2), p)
dep_graph = buildRxDependencyGraph(nu)
jprob = JumpProblem(oprob, agg, jumps...; dep_graph = dep_graph, save_positions = save_positions)
jsol = @time solve(jprob, algos[1]);

pdmpprob = PDMP.PDMPProblem(
	(xdot, xc, xd, p, t) -> F_synapse(xdot, xc, xd, p, t, events_bap, bap_by_epsp),
	(rate, xc, xd, p, t, sum_rate) -> R_synapse(rate, xc, xd, p, t, sum_rate, glu),
	nu, xc0, xd0, p_synapse, (t1, t2);
	Ncache = 12) # this option is for AD in PreallocationTools
pdmpsol = @time solve(pdmpprob, algos[1])

# simulate the spine dynamics
result = @time evolveSynapse(
	xc0,
	xd0,
	p_synapse,
	events_sorted_times,           # external events
	is_pre_or_post_event,   # pre or post?
	bap_by_epsp,
	[true],
	nu,
	algos,
	agg;
	save_positions = save_positions,
)

@test ~isnothing(result)

# # plot the discrete variables
# Synapse.plot_discrete(result.t, result.XC, result.XD)
#
# # plot specific variable
Synapse.plot_variable(result.t, result.XC, result.XD, :Vsp; xlim = (480., 520))
