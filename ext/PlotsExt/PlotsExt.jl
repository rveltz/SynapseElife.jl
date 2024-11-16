module PlotsExt
    using Plots, SynapseElife, DocStringExtensions
    import SynapseElife: plot_discrete,
                         plot_variable,
                         get_names,
                         statistics_jumps
    include("plot.jl")
end
