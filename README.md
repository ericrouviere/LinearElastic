# LinearElastic

A Julia package for a two-dimensional elastic network model of a protein that is evolved in silico. Each residue is a node of a spring network, and the sequence sets the stiffness, rest-length mismatch and contact energy of every spring. Ligands bind as a spring across an active site. Populations evolve by Wright–Fisher dynamics under selection for binding that switches between two ligands. The question is whether this selection for evolvability produces latent allostery in a model that has no built-in bistable switch.

[`methods.md`](methods.md) describes the model and the simulations in detail.

Written by Eric Rouviere, started October 2024. Sister project of [SpinModel](https://github.com/ericrouviere/SpinModel).

## Installation

LinearElastic needs Julia 1.12 (install it with [juliaup](https://github.com/JuliaLang/juliaup)). Plotting uses matplotlib through [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl). PyPlot usually installs matplotlib itself; if PyCall uses your system Python, run `pip install matplotlib`.

The package is not in the Julia General registry, so it installs from GitHub. To use it as a package, run this from the Julia REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/ericrouviere/LinearElastic.git")
using LinearElastic
```

To read or change the source code, clone the repository instead:

```bash
git clone https://github.com/ericrouviere/LinearElastic.git
cd LinearElastic
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then load it with `julia --project=.` and `using LinearElastic`.
