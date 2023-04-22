# using Revise, Plots
using Test
using Synapse
using PiecewiseDeterministicMarkovProcesses, JumpProcesses, OrdinaryDiffEq, Sundials, LSODA
const PDMP = PiecewiseDeterministicMarkovProcesses

##### Parameters
p_synapse = SynapseParams(t_end = 1000.);
glu = 0.0;
events_sorted_times = [500.];
is_pre_or_post_event = [true];
events_bap = events_sorted_times[is_pre_or_post_event .== false];
bap_by_epsp = Float64[];
t1 = 0.;
t2 = 500.;
nu = buildTransitionMatrix();

##### Initial conditions
xc0 = initial_conditions_continuous_temp(p_synapse);
xd0 = initial_conditions_discrete(p_synapse);
u = vcat(xc0, xd0)

##### Jump problem

jsave_positions = (false, true);

# jalgos = (Tsit5(), Tsit5());
# jalgos = (TRBDF2(), TRBDF2());
# jalgos = (lsoda(), lsoda());
# jalgos = (CVODE_BDF(), CVODE_BDF());
jalgos = (AutoTsit5(Rosenbrock23()), AutoTsit5(Rosenbrock23()));

p, jumps = J_synapse(p_synapse, nu, glu);
oprob = ODEProblem((du, u, p, t) -> G_synapse(du, u, p[end], t, events_bap, bap_by_epsp), u, (t1, t2), p);
dep_graph = buildRxDependencyGraph(nu);

# Coevolve

coagg = CoevolveSynced();

coprob = JumpProblem(oprob, coagg, jumps...; dep_graph = dep_graph, save_positions = jsave_positions);
cosol = solve(coprob, jalgos[1]);
cosol = @time solve(coprob, jalgos[1]);

coresult = @time evolveSynapse(
	xc0,
	xd0,
	p_synapse,
	events_sorted_times,    # external events
	is_pre_or_post_event,   # pre or post?
	bap_by_epsp,
	[true],
	nu,
	jalgos,
	coagg;
	save_positions = jsave_positions,
);

# Direct

diagg = Direct();

diprob = JumpProblem(oprob, diagg, jumps...; save_positions = jsave_positions);
disol = solve(diprob, jalgos[1]);
disol = @time solve(diprob, jalgos[1]);

diresult = @time evolveSynapse(
	xc0,
	xd0,
	p_synapse,
	events_sorted_times,    # external events
	is_pre_or_post_event,   # pre or post?
	bap_by_epsp,
	[true],
	nu,
	jalgos,
	diagg;
	save_positions = jsave_positions,
);

##### PDMP problem

pdmpsave_positions = (false, true);

pdmpagg = nothing;

# pdmpalgos = (CHV(:lsoda), CHV(:lsoda));
# pdmpalgos = (CHV(CVODE_BDF()), CHV(CVODE_BDF()));
pdmpalgos = (CHV(AutoTsit5(Rosenbrock23())), CHV(AutoTsit5(Rosenbrock23())));

pdmpprob = PDMP.PDMPProblem(
	(xdot, xc, xd, p, t) -> F_synapse(xdot, xc, xd, p, t, events_bap, bap_by_epsp),
	(rate, xc, xd, p, t, sum_rate) -> R_synapse(rate, xc, xd, p, t, sum_rate, glu),
	nu, xc0, xd0, p_synapse, (t1, t2);
	Ncache = 12) # this option is for AD in PreallocationTools
pdmpsol = solve(pdmpprob, pdmpalgos[1]);
pdmpsol = @time solve(pdmpprob, pdmpalgos[1]);

pdmpresult = @time evolveSynapse(
	xc0,
	xd0,
	p_synapse,
	events_sorted_times,    # external events
	is_pre_or_post_event,   # pre or post?
	bap_by_epsp,
	[true],
	nu,
	pdmpalgos,
	pdmpagg;
	save_positions = pdmpsave_positions,
);

@test ~isnothing(result);

##### Plots
using Plots

plot(coresult.t, coresult.XD[1, :], label="CoevolveSynced");
plot!(diresult.t, diresult.XD[1, :], label="Direct");
plot!(pdmpresult.t, pdmpresult.XD[1, :], label="PDMP");
title!("N_ampa")

plot(coresult.t, coresult.XC[1, :], label="CoevolveSynced");
plot!(diresult.t, diresult.XC[1, :], label="Direct");
plot!(pdmpresult.t, pdmpresult.XC[1, :], label="PDMP");
title!("Vsp")

# # plot the discrete variables
# Synapse.plot_discrete(result.t, result.XC, result.XD)
#
# # plot specific variable
# Synapse.plot_variable(result.t, result.XC, result.XD, :Vsp; xlim = (480., 520))
