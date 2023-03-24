# Synapse.jl

This is the main repository for the Synapse simulation package written in Julia language. The associated publication is here on [BiorXiv](https://www.biorxiv.org/content/10.1101/2021.03.30.437703v1).

## Installation

To run the package `Synapse`, we need some bleeding edge versions of other packages.

### Installation from website

To install the package, please run

```julia
pkg> add PiecewiseDeterministicMarkovProcesses#master
pkg> add https://gitlab.inria.fr/rveltz/synapsepaper.git
```

You can then use it like

```julia
using Synapse
```


### Packages for Figure1 and Figure2

The following packages must be installed (only once) to run `examples/Figure1.jl`:

```julia
pkg> add Revise Random Plots ColorSchemes Parameters
```

The following packages must be installed (only once) to run `examples/Figure2.jl`:

```julia
pkg> add Revise Random Plots ColorSchemes Parameters
```
