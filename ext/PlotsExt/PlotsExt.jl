module PlotsExt
    using Plots, SynapseElife, DocStringExtensions
    import SynapseElife: plot_discrete,
                         plot_variable
    include("plot.jl")
end
