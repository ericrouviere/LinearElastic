# LinearElastic

A Julia package for a two-dimensional elastic network model of a protein that is evolved in silico. It asks whether selection in fluctuating environments, which selects for evolvability, produces latent allostery in a model that has no built-in bistable switch.

Written by Eric Rouviere, started October 2024. It succeeds an earlier "v15" model and is a sister project of [SpinModel](https://github.com/ericrouviere/SpinModel).

## The model

A protein is an elastic network of $N$ nodes (residues) on a randomly perturbed triangular lattice (115 nodes by default), with springs between neighboring nodes. Each node carries one of $q$ amino acid types, and the sequence $\mathbf{s}$ sets three properties of each spring $\langle ij \rangle$ through fixed random tables: its stiffness $k_{ij}$, its rest-length mismatch $d_{ij}$ and a contact energy $\varepsilon_{ij}$. For node displacements $\mathbf{u}$, the energy with ligand $\lambda$ bound is

$$U(\mathbf{u}; \mathbf{s}, \lambda) = \sum_{\langle ij \rangle} \left[ \tfrac{1}{2} k_{ij} \left(e_{ij}(\mathbf{u}) - d_{ij}\right)^2 + \varepsilon_{ij} \right] + \tfrac{1}{2} \kappa \left(\mathbf{b} \cdot \mathbf{u} - \lambda\right)^2,$$

where $e_{ij}$ is the extension of spring $\langle ij \rangle$ and the last term is a ligand spring across a two-node active site. The energy is quadratic, so the free energy $F(\mathbf{s}, \lambda)$ has a closed form: the ground-state energy plus a fluctuation term $\frac{1}{2\beta} \ln \mathrm{Det}'(\beta H / 2\pi)$.

- **Binding.** The solvent ($\lambda_0 = 0$) and two ligands that open ($\lambda_1 > 0$) or close ($\lambda_2 < 0$) the active site define the binding free energies $\Delta F_k = F(\mathbf{s}, \lambda_k) - F(\mathbf{s}, \lambda_0)$. Binding one ligand tightly costs binding of the other.
- **Evolution.** A population evolves by Wright–Fisher dynamics with selection weights $e^{-\alpha \Delta F_k}$, where the selected ligand $k$ switches every $\tau$ generations. An optional stability threshold gives zero weight to sequences whose unbound free energy is above it.
- **Speed.** Mutations and ligands change only a few springs, so the code updates free energies with low-rank (Woodbury) corrections instead of re-inverting the Hessian. A free energy takes microseconds, a mutation about 40 μs, and a generation at population size 100 a few milliseconds.

[`methods.md`](methods.md) describes the model and the simulation protocol in full.

## Requirements

- **Julia 1.12.** The easiest install is [juliaup](https://github.com/JuliaLang/juliaup):
  ```bash
  curl -fsSL https://install.julialang.org | sh    # macOS / Linux
  winget install --name Julia --id 9NJNWW8PVKMN -e # Windows
  ```
  `Manifest.toml` is resolved for Julia 1.12. Other versions should work after `Pkg.resolve()`, but are untested.
- **matplotlib**, used for plotting through [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl). You usually don't need to do anything: PyPlot installs matplotlib into a private Conda Python the first time it loads. If PyCall is set up to use your system Python instead (the default on Linux), install it there with `pip install matplotlib`.

## Installation

LinearElastic is not in the Julia General registry, so it installs from its GitHub URL. Pick one of the two options below.

### Option 1: add it as a package

Use this if you only want to use LinearElastic from your own scripts. From the Julia REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/ericrouviere/LinearElastic.git")
```

Or press `]` to open the package prompt and run:

```
pkg> add https://github.com/ericrouviere/LinearElastic.git
```

This installs LinearElastic and its dependencies into your active environment. The first `using LinearElastic` precompiles the package, which takes a few minutes. Update to the latest version with `Pkg.update("LinearElastic")`.

It is good practice to give each project its own environment instead of the global one. Activate it before adding:

```julia
using Pkg
Pkg.activate("path/to/my_project")   # creates Project.toml there if needed
Pkg.add(url="https://github.com/ericrouviere/LinearElastic.git")
```

### Option 2: clone the repository

Use this if you want to read or change the source code, or run the tests.

```bash
git clone https://github.com/ericrouviere/LinearElastic.git
cd LinearElastic
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Pkg.instantiate()` installs the exact dependency versions recorded in `Manifest.toml`.

To use your clone from a different environment, add it with `Pkg.develop` instead of `Pkg.add`. The environment then tracks your local files, so edits to the source take effect:

```julia
using Pkg
Pkg.develop(path="path/to/LinearElastic")
```

## Loading LinearElastic

### In the REPL or a notebook

```julia
using LinearElastic
```

When editing the source, load [Revise](https://github.com/timholy/Revise.jl) first so changes are picked up without restarting Julia:

```julia
using Revise
using LinearElastic
```

### In a script

If you installed with Option 1, or used `Pkg.develop`, activate that environment and load the package:

```julia
using Pkg
Pkg.activate("path/to/my_project")
using LinearElastic
```

If you cloned the repository (Option 2), activate the clone itself. Building the path from `@__DIR__`, the folder of the script, keeps the script working wherever it is run from:

```julia
using Pkg
Pkg.activate(joinpath(@__DIR__, "path/to/LinearElastic"))   # path relative to this script
using LinearElastic
```

Alternatively, leave `Pkg.activate` out of the script and choose the environment on the command line:

```bash
julia --project=path/to/LinearElastic my_script.jl
```

Setting the `JULIA_PROJECT` environment variable to the same path has the same effect.

### With parallel workers

The ensemble functions (`evolvePopEnsemble`, `evolvePopStaticEnsemble`, `computeBindingDMSEnsemble`) run on worker processes through `pmap`. Each worker must start in the same environment and load LinearElastic:

```julia
using Distributed
addprocs(4; exeflags="--project=$(Base.active_project())")
@everywhere using LinearElastic
```

Each replicate seeds its worker's random number generator, so ensemble results do not depend on how replicates are scheduled across workers.

## Quick start

This script builds a random network, computes the binding free energies of a random sequence, and evolves a population under selection that switches between the two ligands. It runs in about half a minute, most of it compilation.

```julia
using Random
using LinearElastic

Random.seed!(1)

# Model: an 11 × 11 triangular lattice of 115 nodes and random interaction tables
W, L, offset = 11, 11, 0.15          # lattice size and random displacement of the nodes
q = 5                                # number of amino acid types
k_min, k_max = 1e-2, 1e1             # range of spring stiffnesses
energyVar, mismatchVar = 0.0, 0.1    # spread of contact energies and rest-length mismatches
β = 1.0                              # inverse temperature
tables = makeTables(W, L, offset, q, k_min, k_max, energyVar, mismatchVar)

# Solvent and two ligands: a spring of stiffness 1 across the active site (nodes 47 and 59)
# that holds it at its reference length (solvent), opens it (ligand 1) or closes it (ligand 2)
site = (47, 59)
ligands = Ligands([Perturbation([site], [1.0], [0.0]),    # 1: solvent (unbound)
                   Perturbation([site], [1.0], [+0.2]),   # 2: ligand 1
                   Perturbation([site], [1.0], [-0.2])])  # 3: ligand 2

# Free energies of a random sequence
net = buildRandNetwork(tables)
F = computeFreeEnergies(net, ligands, tables, β)   # [F_solvent, F_1, F_2]
ΔF = F[2:3] .- F[1]                                # binding free energies; negative = binds

# Evolve a population under selection that switches between the ligands every τ = 16 generations
_, ΔF_0, _ = computeBindingRandSeqs(1000, ligands, tables, β)   # scale of binding energies
evo = EvoParams(100, 512, 3e-3, 16, 1 / ΔF_0)                  # P, generations, μ, τ, α
population = [buildRandNetwork(tables) for _ in 1:evo.P]
fitBucket, sampled, F0Bucket = evolvePop!(population, ligands, tables, β, evo; N_seqs2Keep=50)

# fitBucket[n, k, t]: fitness -ΔF_k of individual n at generation t (P × 2 × generations)
meanFitness(t) = round.(vec(sum(fitBucket[:, :, t], dims=1)) ./ (evo.P * ΔF_0), digits=2)
println("ΔF of a random sequence: ", ΔF)
println("mean fitness (ligand 1, ligand 2) in units of ΔF_0: ", meanFitness(1), " → ", meanFitness(evo.N))
```

The last generation falls in a ligand-2 phase of the cycle, so the population ends up binding ligand 2 better than ligand 1.

From here:

- **Other selective pressures.** `evolvePopStatic!(nets, ligands, tables, β, evo, assay)` selects on a fixed `Assay`: `Binding1()`, `Binding2()`, `DoubleBinding()` (the generalist score $\min(-\Delta F_1, -\Delta F_2)$), `TargetBinding(ΔF1_target, ΔF2_target)` or `Allostery()`.
- **Stability selection.** `EvoParams(P, N, μ, τ, α, F_unfolded)` adds a stability threshold: a sequence whose solvent-state free energy is above `F_unfolded` cannot reproduce.
- **Replicates.** `evolvePopEnsemble(ligands, tables, β, evo, seeds)` evolves one seeded population per seed in parallel.
- **Mutational scans.** `computeBindingDMS(net, ligands, tables, β)` returns the effect of every single mutation on both binding energies. `computeConfChange` returns the conformational change on binding.
- **Plots.** `plotNetwork!`, `plotBondStiffness!`, `plotDisplacement!` and `plotMutSensVsDist!` draw on a matplotlib axis.

## Repository layout

| Path | Contents |
|------|----------|
| `src/types.jl`, `src/buildnetwork.jl` | Types (`Tables`, `Network`, `Perturbation`, `Ligands`, `EvoParams`) and network construction |
| `src/physics.jl` | Energy parameters, free energies and low-rank updates |
| `src/fitness.jl` | Selective pressures (`Assay` types) and the stability threshold |
| `src/mutate.jl`, `src/evolve.jl` | Mutations and Wright–Fisher evolution |
| `src/ensemble.jl` | Parallel replicate populations and mutational scans |
| `src/analysis.jl`, `src/plotfuncs.jl` | Mutational scans, sensitivity, stress and normal-mode analysis; plotting |
| `test/` | Tests |
| `methods.md` | Description of the model and simulation protocol |

## Tests

From a clone of the repository:

```bash
julia --project=. test/test_physics.jl   # low-rank updates against building from scratch, mutational scans
julia --project=. test/test_evolve.jl    # assays, selection, seeded and parallel ensembles
```

`test/test_rigidpy.jl` checks the network Hessian against the Python package [Rigidpy](https://github.com/VardaHagh/Rigidpy). It is skipped unless Rigidpy is cloned into `Rigidpy/` at the repository root.
