using Documenter, Synapse, Setfield

makedocs(doctest = false,
	sitename = "Model of excitatory synapse in Julia",
	format = false,
	pages = Any[
		"Home" => "index.md",
		# "Simple example" => "tutorials.md",
		# "Figure 1" => "figure1.md",
		"Library" => "library.md"
	]
	)
	
deploydocs(
	repo = "github.com/rveltz/SynapseElife.git",
	devbranch = "main"
)
