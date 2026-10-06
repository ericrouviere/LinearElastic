# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

LinearElastic is a Julia package (author Eric Rouviere, started Oct 2024) implementing a 2D linear elastic network
model of a protein that is evolved in silico. The question: does selection in fluctuating environments (selection for
evolvability) produce latent allostery in a model with no bistable switch? It is the successor of an earlier "v15"
model, and the sister project of `../SpinModel`, whose organization, naming and experiment workflow it follows.

## Julia

```julia
# Load the package (from within Julia REPL, at repo root)
using Pkg; Pkg.activate("."); using LinearElastic

# During development, use Revise so edits are picked up without restart
using Revise
using LinearElastic
```

`tutorial.ipynb` (repo root, Julia 1.12 kernel, run from the repo root) walks through the model for a beginner: tables,
ligands, a random sequence, its energy landscape and free energy, binding energies (checked against the closed form),
mutations, a fluctuating-selection run, the conformational change and the mutational sensitivity. It is committed with
its outputs; re-execute it with `jupyter nbconvert --to notebook --execute --inplace tutorial.ipynb` after API changes.

The environment is resolved for Julia 1.12 (installed via juliaup). Tests are standalone scripts in `test/`:

```bash
julia --project=. test/test_physics.jl   # low-rank updates vs. building from scratch, DMS
julia --project=. test/test_evolve.jl    # assays, fitness map, seeded/parallel ensembles
```

### Code Style & Idioms

**Structure**
- Divide programs into functions early (reusable, testable, faster); avoid long top-level blocks. Functions take arguments, not globals (except constants like `pi`).
- Organize code into modules with explicit `export`; treat exported methods as the interface and fields as implementation details.
- Use multiple dispatch (structs + methods on types). Avoid Python-style OOP and Python syntax (`def`, `list.append()`).
- 4 spaces per indent. No parens around `if`/`while` conditions.

**Arrays & loops**
- Julia is strictly 1-indexed; never use 0-indexed arrays.
- Prefer broadcasting (`.+`, `f.(x)`) over explicit loops where natural.
- Prefer `eachindex(x)` over `1:length(x)`; use `@views` for array slices.
- Don't overuse splatting: `[a; b]` not `[a..., b...]`.

**Types & dispatch**
- Write type-stable code.
- Keep arguments generic; don't over-restrict. `f(A::Matrix{Float64})` silently excludes sparse, GPU, Float32, dual numbers. Use `AbstractMatrix` when 2-D structure matters, `AbstractArray` when it doesn't, unannotated when any input works. Prefer abstract args (`Integer`) over concrete (`Int`); Julia specializes automatically (no perf penalty). Annotate to control dispatch and resolve ambiguities, not to document intent.
- Don't use unnecessary static params: `foo(x::Real)`, not `foo(x::T) where {T<:Real}`, when `T` is unused.
- Avoid strange `Union`s (`Union{Function,AbstractString}`) and elaborate container types (prefer `Vector{Any}` over packing alternatives; annotate uses like `a[i]::Int`).
- Default to instance dispatch over `Type{T}`; model enumerations as instances of one immutable type.
- Test types with `isa`/`<:`, not `==`. Constructors `T(x)` must return a `T`.
- Handle excess argument diversity in the caller (convert at the call site), not inside the function.

**Mutation & naming**
- `!` suffix for mutating functions (`sort!(x)`); they typically also return the mutated arg. IO/RNG exception: `!` marks mutation *beyond* the inherent IO/RNG mutation (`read!(io, x)` mutates both).
- Types/modules `CamelCase`, constants `UPPER_SNAKE`, functions `snake_case`. `_` prefix = internal. Avoid abbreviation (`indexin`, not `indxin`).

**Functions & calls**
- Don't wrap a named function: `map(f, a)`, not `map(x->f(x), a)`.
- Avoid redundant keyword syntax: when a variable name matches the keyword name, use `f(; max_iter)` not `f(; max_iter=max_iter)` (applies at call sites, NamedTuple construction, etc.).
- Argument order (Base convention): function, IO, mutated input, type, input, key, value, varargs, kwargs.

**Numeric**
- Avoid float literals in generic numeric code (they promote and disrupt input types); prefer `Int`, then `Rational{Int}`.

**Docs & safety**
- Write docstrings (`"""..."""`) above definitions for public API; use `#`/`#= =#` for reader-only notes.
- Prefer avoiding errors over relying on `try`/`catch`.
- Don't overuse macros (prefer functions); `eval` inside a macro is a red flag.
- Don't overload Base methods on container types you don't own; avoid type piracy (extending Base/others' methods on foreign types).
- Don't expose unsafe operations at the interface; put `unsafe` in the name or check safety.

## Core data types (`src/types.jl`, `src/fitness.jl`)

| Type | Role |
|------|------|
| `Tables` | Fixed for an experiment: bead positions `R = [x1,y1,x2,y2,...]`, adjacency `A`, stiffness table `K`, mismatch table `D`, contact-energy table `E` (each bond adds `E[s_i,s_j]` to `c`; cancels in binding energies). All three must be symmetric (`checkTables`); to use chosen values, build `Tables(t.R, t.A, t.K, t.D, E)` |
| `Network` | One evolving sequence and its cached energy parameters (of `U = ½ r'Hr − h'r + c`) `G = pinv(H)`, `lpdH = log pdet(H)`, `h`, `c`. The `δ*` fields are scratch buffers |
| `Perturbation` | `bonds`, `Δk`, `Δd`: changes in stiffness and rest-length mismatch of a set of bonds. Used for ligands and mutations |
| `Ligands` | `list::Vector{Perturbation}`; `list[1]` is the solvent, `list[2]`, `list[3]` the two environments |
| `EvoParams` | Evolutionary parameters: `P` (population), `N` (generations), `μ` (per-site mutation rate), `τ` (env timescale), `α` (selection strength). Built from settings with `convertEvoParams(S)` |
| `Assay` | Selective pressure, passed as a value: `Binding()`, `DoubleBinding()`, `Allostery()`, `TargetBinding(ΔF1_target, ΔF2_target)` |

## Architecture

### Physics (`src/buildnetwork.jl`, `src/physics.jl`)
- Beads sit on a perturbed triangular lattice of size `W`×`L` (`L` odd). The default 11×11 gives 115 beads. 2D only
  (`N_sdims=2` is asserted in places).
- Each bead has a type in `1:N_types`, and a sequence is a `Vector{Int}` of types. The bond between beads i and j has
  stiffness `K[seq[i],seq[j]]` (log-uniform in `[k_min,k_max]`) and mismatch `D[seq[i],seq[j]]` (Gaussian, sd `mismatchVar`).
- Energy: `U(r) = ½ r'Hr − h'r + c` (the sum over bonds of `½ k (e − d)²`, `e` the bond elongation) with `H = T'ΓT`, `h = T'ΓM·1`, `c = ½·1'ΓM²1`. Here `T` is the beads-to-bonds matrix
  (`makeBeads2BondsMatrix`), `Γ = diag(k)` and `M = diag(d)`. The free energy is the ground-state energy (`-½h'Gh + c`)
  plus the fluctuation term from `log pdet(H)` (3 zero modes in 2D). Unlike SpinModel, `β` is a real parameter here:
  it weighs the fluctuation term against the ground-state energy.
- **Low-rank perturbations (the hot path).** A `Network` caches `G = H⁺`, `lpdH`, `h`, `c`, the ground state `u = G h` and
  `hGh = h'G h`. `perturbedFreeEnergy(net, bonds, Δk, Δd, tables, β)` gives the free energy after changing m bonds from
  the m×m block `T G T'` and `T u` (Woodbury + determinant lemma, written as `log det(I + ΔS)`), in O(m³) without
  touching `G`; `computeFreeEnergy(net, per, ...)` and `computeFreeEnergies` use it. `applyPerturbation!` updates all
  cached parameters in place with one rank-m BLAS update of `G`; `mutateBead!` = `mutationPerturbation` +
  `applyPerturbation!` (+ contact energy). `computeBindingDMS` evaluates each mutant×ligand state as one merged
  perturbation of the wild type (`mergePerturbations`), never updating `G`. `computePerturbedNetParams!` (full `δG`, writes
  the `δ*` scratch fields) is kept only for `computeConfChange`/`computePerturbedGroundState`.
- Typical costs (11×11 network, 1 BLAS thread): `buildNetwork` 6.4 ms (SVD), 3-state free energies 2 μs, mutation 41 μs,
  DMS 10 ms per sequence, evolution ~3.6 ms per generation at P = 100.
- Binding energy is `ΔF = F_lig − F_solvent` (negative = good binder). Population fitness in environment k is
  `F_solvent − F_lig_k` (positive = good binder).

### Fitness assays (`src/fitness.jl`)
Selective pressures are `Assay` types passed as values, not strings. `computeFitness(net, ligands, tables, assay, β)`
dispatches on the type; add a pressure with a new struct plus one `computeFitness(::NewAssay, energies)` method.
`Binding` needs 2 ligand states, `Binding1`/`Binding2` (binding of ligand 1 or 2), `DoubleBinding` (the generalist score
`min(-ΔF1, -ΔF2)`) and `TargetBinding` need 3, `Allostery` needs 4.

**Folding.** `foldedFraction(F, F_unfolded, β) = 1/(1 + exp(β(F − F_unfolded)))`, with `F = F(s, λ0)` the solvent-state free
energy and `F_unfolded` the unfolded-state free energy, a physical parameter passed as the keyword `F_unfolded` to every
population-evolution function (default `Inf`: never unfolds, identical to no unfolding). Selection weights are
`p_fold · exp(α f)` (`fits2weights(f, α, logFoldedFraction)`). `F_unfolded` must be chosen relative to the free energies the
code computes (e.g. a quantile of random sequences' `F(s, λ0)`); the code's constant counts all 2N coordinates in
`n log(2π/β)` (methods Eq. 4 uses 2N − 3), a sequence-independent shift of −(3/2β) log(2π/β).

### Evolution (`src/evolve.jl`)
- `evolveMC(net, ligands, tables, assay, N_iter, ζ, β; Δt, ϕ_terminate)`: single-trajectory Metropolis walk under an `Assay`.
- `evolvePop!(nets, ligands, tables, β, evo; k_init, sampleLastHalf, N_seqs2Keep, F_unfolded)`: population of `evo.P`
  networks under fluctuating selection; the environment switches between ligand 1 and 2 (`Binding1`, `Binding2`) every
  `τ` generations. Returns `(fitBucket, netBucket, F0Bucket)`: fitness `P × 2 × N`, `N_seqs2Keep` networks sampled from
  the last half of the run (aligned to the env cycle), and the solvent-state free energy `P × N`.
- `evolvePopStatic!(nets, ligands, tables, β, evo, assay; ...)`: static selection on one `Assay` (τ unused);
  `fitBucket` is `P × 1 × N`. Both share `_evolvePopulation!`.
- `evolvePop(ligands, tables, β, evo, seed; ...)` / `evolvePopStatic(..., assay, seed; ...)`: seeded runners that build a
  random population, evolve it and return only the sampled sequences.
- Selection is copy-on-write (`applySelection!` + `mutateCopyOnWrite!`): extra draws share their parent's `Network` object
  and get a private copy (recycling an undrawn network from `spare`) only when they mutate; `_unshare!` gives every slot
  its own network at the end. The slot layout matches the old copy-based selection (first draws keep their slot, extra
  draws fill undrawn slots in order), so the same seed gives the same trajectory as before the speedups. Never mutate
  a population member in place inside the loop without going through `mutateCopyOnWrite!`.

### Ensemble (`src/ensemble.jl`)
Naming convention: append `Ensemble` to the single-run function name (`evolvePopEnsemble`, `computeBindingDMSEnsemble`).
- `evolvePopEnsemble(ligands, tables, β, evo, seeds; k_init_vec, N_seqs2Keep, F_unfolded, pool)` evolves one replicate per seed over
  a Distributed worker pool. Each replicate seeds its worker's RNG, so results do not depend on scheduling. `k_init_vec`
  defaults to a random initial environment per replicate. Output is `Vector{Vector{Vector{Int}}}` (rep → seq → beads);
  `convertPopulations2PopArray` turns it into an `Int8` array `(N_beads, N_seqs2Keep, N_reps)` for saving.
- `evolvePopStaticEnsemble(ligands, tables, β, evo, assay, seeds; N_seqs2Keep, F_unfolded, pool)`: static counterpart.
- `computeBindingDMSEnsemble(populations, ligands, tables, β; pool)` runs one flat `pmap` over all sequences (batched,
  since a scan takes ~10 ms) and returns `DMS[type, bead, lig, seq, rep]` and `ΔF[lig, seq, rep]`.
- All `pmap`s go through `withPool(f, pool)`: with `pool = nothing` (default) it uses a temporary `CachingPool(workers())`
  and clears it, so the closure and the `tables`/`ligands` it captures ship to each worker once; pass your own pool to
  reuse it across calls (as `run.jl` does). Tasks return only small results (sequences, scans). Workers use 1 BLAS
  thread, so parallel DMS values can differ from a serial run by ~1e-13; evolved populations are bit-identical.

### Analysis and plotting (`src/analysis.jl`, `src/plotfuncs.jl`)
DMS and double-mutant scans, mutational sensitivity vs. distance from the active site, stress and normal-mode overlaps.
Plotting uses PyPlot (matplotlib via PyCall); in-place variants take an `ax` and end in `!`.

## Running simulations (`exps/`)

A simulation lives in a numbered experiment directory `exps/NNN/` with a `scripts/` subdir holding `make_settings.jl`,
`run.jl`, `run.apply` and `run.sbatch`. New experiments are made by copying the latest one (011 is the current template).
Run from inside `exps/NNN/scripts/`:

```bash
julia make_settings.jl   # build the settings.jld2 files
bash run.apply           # submit one job per ens directory (sbatch); swap in `julia run.jl` to run locally
```

1. **`make_settings.jl`** activates the project via `@__DIR__`, builds seeded `Tables` (`table_seed`) and the ligands,
   sets the evolutionary params and run controls (`seed` = master RNG seed, `N_workers`, `N_reps`, `N_seqs2Keep`), and
   writes one `ens/<NNNN>/scripts/settings.jld2` per (τ, μ) pair. Selection strength is `α = a / ΔF_0`, with `ΔF_0`
   the std of random-sequence binding energies.
2. **`run.apply`** exports `JULIA_PROJECT` (repo root, derived from the script's location), copies `run.jl` and
   `run.sbatch` into each `ens/<NNNN>/scripts/`, and submits. `--cpus-per-task` in `run.sbatch` should be `N_workers + 1`.
3. **`run.jl`** loads `settings.jld2`, calls `addprocs(S["N_workers"])`, derives per-replicate seeds and initial
   environments from the master seed, evolves `N_reps` populations in parallel, then runs the DMS. Each analysis phase
   is timed and isolated.

Outputs land in `ens/<NNNN>/data/`: `populations`, `dms`, `binding_energies`, plus `seeds` (master seed, per-replicate
seeds, `k_init_vec`) and `run_meta` (Julia version, worker count, git commit, timestamp). To reproduce a run, load
`seeds.jld2` and call `evolvePopEnsemble(ligands, tables, β, evo, seeds; k_init_vec, N_seqs2Keep)`.

**Analysis**: `exps/NNN/analysis/make_figs.jl` is a headless script version of the analysis notebook (011 has one, ported
from 010's `analysis_before_after.ipynb`). Run it from that directory with `julia make_figs.jl > make_figs.log 2>&1`
(`N_WORKERS` env var, default 8). It activates the project via `@__DIR__`, reads τ/μ and the tables from the ens
settings instead of hard-coding them, fails fast if any run output is missing, seeds its random-sequence baseline,
writes every figure as `.svg` and `.pdf`, and prints a summary table. Gotchas: paths given to `pmap` must be absolute
(`cd` only moves the main process, not the workers), and with `text.usetex` on in matplotlibrc, labels must not contain
literal unicode or a bare `|` (use math mode, e.g. `L"Stress $|h|$"`).

The cluster is UChicago **midway** (account `pi-ranganathanr`, partition `caslake`, remote path
`/home/rouviere/LinearElastic`); `exps/download_data_midway.sh` scp's `ens/` back. `exps/README` is the experiment log:
- All use W=L=11, β=1, mismatchVar=0.1.
- 001/002: N_types=10, N_gen=1024, μ=1e-3, τ=4 fixed.
- 003: N_types=5. τ sweep from here on.
- 004: N_gen=2048, two-bond site `(4,7),(4,16)`.
- 005: single bond site `(4,7)`.
- 006: k_min=1e-2.
- 007: `selParams=(0.25/s,1)`.
- 008: k_min=1e-2, μ=3e-3, N_reps=100, N_seqs2Keep=50.
- 009: same as 008, but the solvent reference has Δk=0.2 instead of 1.
- 010: same as 008, but with a single interior bond site `(47,59)`.
- 011: new template, same physics as 010. Not yet run.

**001–010 and the notebooks use the old API** (`selParams` sigmoid selection, string assays, `*Ens` names, positional
evolution parameters). Run them against the `old-api` git tag (`git worktree add ../LinearElastic-old old-api`).
Their saved data is unaffected.

## Methods write-up
`methods.md` describes the model and the exps/011 protocol as a journal methods section (energy, free energy,
ligand binding, phenotype, evolutionary dynamics, parameter table). With the current ligands every state inserts the
same spring stiffness, so the fluctuation term cancels and binding depends on the sequence only through the site's
pre-strain `e₀ = b'G'h` and compliance `χ = b'G'b`; keep the note in sync when the model or protocol changes.

## Bugs fixed (Oct 2026) and their effect on existing data
- `updateNetParams!` now sets `net.lpdH`. It used to be a local no-op. Binding energies were never affected because
  `lpdH` cancels in them.
- `computeBindingRandSeqs` now returns `(ΔF, std, mean)`, the order every caller unpacks (`ΔF, s, m`). Experiments
  001–010 scaled selection by the **mean** of the random-sequence binding cloud, not its std. Their stored
  `S["selParams"]` and `S["stdBindRandSeqs"]` (which actually hold the mean) record what was really used.
- 001–010 started every replicate in environment 1 (`k_init_vec` was two `ones(...)` blocks).
- The sigmoid selection of 001–010 equals the current exponential map with `α = a` whenever `a (max f − x0) < −100`
  (e.g. 010). Experiment 009 (`a ≈ 0.94`) was in the sigmoid's flat regime, so selection there was weak (fitness
  weights varied by only ~20% across the population).

## Key conventions
- `exps/`, `notebooks/` and `Rigidpy/` (a clone of github.com/VardaHagh/Rigidpy, used by the unit-test notebook) are
  git-ignored large/data directories. Only `src/`, `test/`, the package files, `README`, `CLAUDE.md`, `methods.md` and `tutorial.ipynb` are tracked.
- `src/old/` contains deprecated code kept verbatim for reference (`<file>_old_funcs.jl`, and removed exports in
  `LinearElastic_old_exports.jl`); do not use or modify it.
- `computeHessian`, `computePerturbedNetParams_old`, `pdetdiagm` and the matrix method of `logpdet` are kept only for
  `notebooks/unit_tests/unit_tests.ipynb`.
- `settings.jld2` stores `Tables`/`Ligands`/`Perturbation` structs. Changing those type definitions breaks loading old
  experiment settings with JLD2.

