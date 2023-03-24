# Synapse.jl

| **Documentation** | **Build Status** | 
|:-----------------:|:----------------:|
|  |  [![CI](https://github.com/rveltz/SynapseElife/actions/workflows/ci.yml/badge.svg)](https://github.com/rveltz/SynapseElife/actions/workflows/ci.yml) [![codecov](https://codecov.io/gh/rveltz/SynapseElife/branch/main/graph/badge.svg?token=SQPQFOGJWT)](https://codecov.io/gh/rveltz/SynapseElife)


This is the main repository for the Synapse simulation package written in Julia language. The associated publication is on [BiorXiv](https://www.biorxiv.org/content/10.1101/2021.03.30.437703v1).

### Authors

Yuri Rodrigues, Cian O'Donnell and Romain Veltz wrote the code.

## Installation

> You need at least Julia v1.5.0 to run this code

To run the package `Synapse`, we need some bleeding edge versions of other packages.

### Installation from website

To install the package, please run

```julia
pkg> add https://github.com/rveltz/SynapseElife.git
```

> Note that this has to be done only once! You don't need to do it for every session


You can then use it like

```julia
using Synapse
```

## Website

You can build the website associated with this package by running `make.jl` in the folder `path/Synapse/docs`


## Documentation

If you want to have information about a function use `?`. For example, you can type `?SynapseParams` to get

```julia
help?> SynapseParams
search: SynapseParams

  struct SynapseParams


  Postsynaptic parameters

  Units
  ≡≡≡≡≡≡≡

    •  time: ms

    •  length: µm (area µm^2, volume µm^3)

    •  voltage: mV
...
```

## Parameter definition / list

In order to know the parameters used in the model, you can do for example (for the presynaptic ones):

```julia
using Synapse
?PreSynapseParams
```

which gives

```julia
help?> PreSynapseParams
search: PreSynapseParams

   struct PreSynapseParams


  Presynaptic parameters

    •    Firing events are processed separately from the main simulation (at /Users/rveltz/only_stp.jl) it takes the firing
        structure as input from the function [compilation][/Users/rveltz/utils_data.jl:672].

    •    Using the EPSP times vesicle depletion and AP induced by EPSP is estimated, however one can use a tag in order to
        deactivate it (covering subthreshold EPSP cases).

    •    The presynaptic part of our model is phenomelogical, for instance, soma variable was made to have intended to
        represent the voltage and is able to accumulate for faster frequencies but has an abstract unit.

  Equations
  ≡≡≡≡≡≡≡≡≡≡≡

  based on D. Sterratt et al book; Principles of Computational Modelling in Neuroscience (https://www.compneuroprinciples.org/
  ) page 183

  raterec = (n_rec_t0 - rec) * rec_dt * rrp

  raterrp = (n_rrp_t0 - rrp) * rrp_dt * rec

  rateref = (n_rec_t0 - rec) * ref_dt

  Fields
  ≡≡≡≡≡≡≡≡

    •    trec_ca::Float64

        recovery constant of pre calcium decay function Default: 20000

    •    delta_ca::Float64

        fraction of decay constant of pre calcium decay f Default: 0.0004

    •    tau_pre::Float64

        decay time constant of pre calcium Default: 20

    •    tau_soma::Float64

        decay time constant for AP induced by EPSP Default: 40

    •    n_rrp_t0::Int64

        initial conditions ready releaseble pool Default: 25

    •    n_rec_t0::Int64

        initial conditions recovery pool Default: 30

    •    rec_dt::Float64

        rate rrp -> rec Default: 1 / 5000

    •    rrp_dt::Float64

        rate rec -> rrp Default: 1 / 45000

    •    ref_dt::Float64

        rate infinite reservoir -> rec Default: 1 / 40000

    •    slope::Float64

        sigmoid parameter for release probability Default: 2.0

    •    half::Float64

        sigmoid parameter for release probability Default: 0.7

    •    ap_slope::Float64

        sigmoid parameter for AP-by-EPSP probability Default: 2.0

    •    ap_half::Float64

        sigmoid parameter for AP-by-EPSP probability Default: 0.7
```

## Specifying synapse parameters

If you want to alter the default value of a parameter, you can create a variable as follows

```julia
param_synapse = SynapseParams(
				t_end           = 1000.,
				soma_dist 		= 200.,)
```

