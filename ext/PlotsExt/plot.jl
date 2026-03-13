"""
$(SIGNATURES)

Plotting function

# Arguments
- `tt` times
- `xc` continuous variable
- `xd` discrete variable
- `s::Symbol = :ampa` which variable to plot. Must be `:ampa, :nmda, :vgcc_t, :vgcc_r, :vgcc_l` or `:Vsp,:Vdend,:Vsoma,:λ,:ImbufCa,:Ca,:Dye,:CaM0,:CaM2C,:CaM2N,:CaM4,:mCaN,:CaN4,:mKCaM,:KCaM0,:KCaM2N,:KCaM2C,:KCaM4,:PCaM0,:PCaM2C,:PCaM2N,:PCaM4,:P,:P2,:LTD,:LTP,:LTD_act,:LTP_act,:m,:h,:n,:SK,:λ_age,:λ_aux`.

# Optional arguments
- all arguments from Plots.jl. For example `xlims=(0,10), legend=false`
"""
function plot_variable(tt, xc, xd, s = :ampa; tspan = (0., Inf64), kwargs...)
	out = get_names(xc, xd)
	if s in [:ampa, :nmda, :vgcc_t, :vgcc_r, :vgcc_l]
		st = statistics_jumps(tt, out[s]'; tspan = tspan)
	else
		st = [0, 0]
	end
	plot(tt, out[s]; title = "#jumps ($s) = $(sum(st[2]))", label = "$s", xlims = (tspan[1], min(tt[end], tspan[2])), kwargs...)
end

"""
$(SIGNATURES)

Plot all discrete variables.

# Arguments
- `tt` times
- `xc` continuous variables (result from PDMP)
- `xd` discrete variables (result from PDMP)

# Optional arguments
- all arguments from Plots.jl. For example `xlims = (0, 10), legend = false`
"""
function plot_discrete(tt, xc, xd; tspan = (0., Inf64), kwargs...)
	outd = get_names(xc, xd)
	Plots.plot(layout=(3,2))
	for (ind, s) in enumerate([:ampa, :nmda, :vgcc_t, :vgcc_r, :vgcc_l])
		st = statistics_jumps(tt, outd[s]'; tspan = tspan)
		plot!(tt, outd[s]; title = "#jumps ($s) = $(sum(st[2]))", subplot=ind, label="", titlefontsize = 10, xlims = (tspan[1], min(tt[end], tspan[2])), line = :step, kwargs...) |> display
	end
	if :Vsp in keys(outd)
		plot!(tt, outd[:Vsp]; title = "#jumps (Vsp) = 0", subplot=6, label="", titlefontsize = 10, xlims = (tspan[1], min(tt[end], tspan[2])),  kwargs...) |> display
	end
	Ind = (tt .> tspan[1]) .* (tt .< tspan[2])
	annotate!(3000,1,"total #jumps = $(sum(Ind))",subplot=5)
end