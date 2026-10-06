# CLAUDE.md

## Project overview

Julia package `LinearElastic` (author Eric Rouviere, started Oct 2024). It is a 2D linear elastic network
model of a protein that is evolved in silico. The question: does selection in fluctuating environments
(selection for evolvability) produce latent allostery in a model with no bistable switch? It is the successor of an
earlier "v15" model. Plotting uses PyPlot (matplotlib) and data uses JLD2. Julia 1.12 is installed via juliaup.

### Physics (src/physics.jl, src/buildnetwork.jl)
- Beads sit on a perturbed triangular lattice of size `W`×`L` (`L` must be odd). The default is 11×11, which gives 115 beads.
  Positions are stored as a flat vector `R = [x1,y1,x2,y2,...]` (2D only, so `N_sdims=2` is asserted in places).
- Each bead has a type in `1:N_types`, and a sequence is a `Vector{Int}` of types. The bond between beads i and j has
  stiffness `K[seq[i],seq[j]]` (log-uniform in `[k_min,k_max]`) and rest-length mismatch `D[seq[i],seq[j]]`
  (Gaussian with sd `mismatchVar`). The `E` contact-energy table exists but is unused (`energyVar=0`).
- Energy: `U(r) = ½ r'Hr + h'r + c` with `H = T'ΓT`, `h = T'ΓM·1` and `c = ½·1'ΓM²1`. Here `T` is the beads-to-bonds matrix
  (`makeBeads2BondsMatrix`), `Γ = diag(k)` and `M = diag(d)`. The free energy is the ground-state energy (`-½h'Gh + c`)
  plus the fluctuation term from `log pdet(H)`. There are 3 zero modes in 2D.
- The `Network` struct caches `G = pinv(H)`, `lpdH = log pdet(H)`, `h` and `c`. The `δ*` fields are **scratch buffers**.
  `computeFreeEnergy(net, per, ...)` overwrites `net.δG`/`net.δh` with `G+δG`/`h+δh`.
- A `Perturbation(bonds, Δk, Δd)` changes the stiffness and mismatch of a set of bonds. `computePerturbedNetParams!` applies it
  as a low-rank Woodbury update to `G` plus the matrix determinant lemma for `lpdH`, so no re-SVD is needed. Mutations
  (`mutateBead!`) are also Perturbations on all bonds of the mutated bead. They are applied in place and accumulate
  numerical error over many mutations.
- `Ligands.list[1]` is the **solvent/reference** state, and `list[2]`, `list[3]` are ligands 1 and 2. Binding energy is
  `ΔF = F_lig − F_solvent` (negative means a good binder). Population fitness is `F_solvent − F_lig` (positive means a good binder), and
  column `k` of `fits` is environment k, i.e. ligand k.

### Evolution (src/evolve.jl, src/ensemble.jl)
- `evolveMC`: a single-trajectory Metropolis walk under an assay string (`"Binding"`, `"Generalist"`, `"Cooperativity"`,
  `"TargetBinding"`; see `fitness.jl`). It is used in the notebooks.
- `evolvePop!`: Wright-Fisher-like population of size `P`. Per-site mutation probability is `μ`. The environment alternates between
  ligand 1 and ligand 2 every `τ` generations. Selection weights are a sigmoid `1/(1+exp(-a(f-x0)))` with
  `selParams=(a,x0)`. It samples `N_seqs2Keep` sequences from the last half of the run.
- `evolvePopEns` maps over seeds (serially, not with pmap). Its output is `Vector{Vector{Vector{Int}}}` (rep → seq → beads).
  `convertPopulations2PopArray` turns that into an `Int8` array `(N_beads, N_seqs2Keep, N_reps)` for saving.
- `computeBindingDMSEns` returns `DMS[type, bead, lig, seq, rep]` (ΔΔF of every single mutant) and `ΔF[lig, seq, rep]`.

### Directory layout
- `src/` holds the module. `LinearElastic.jl` includes the files and lists the exports. `old.jl` and `potts.jl` are **not**
  included (dead code).
- `exps/NNN/`: one directory per experiment. Each one has three parts:
  - `scripts/make_settings.jl` builds a settings `Dict S` (including the `Tables` and `Ligands` objects) and **wipes and
    recreates** `../ens/`. It writes `ens/MMM/scripts/settings.jld2` for each value of a swept parameter (`the_list_1`,
    usually `τ ∈ [1,2,4,...,256]`).
  - `scripts/run.apply` copies `run.jl` and `run.sbatch` into each `ens/*/scripts/` and `sbatch`es them. `run.jl` then writes
    `ens/MMM/data/{seeds,populations,dms,binding_energies}.jld2`.
  - `analysis/analysis_before_after.ipynb` loads the `ens` data and makes the result figures.
  - `ens/` is about 4.4 GB, is git-ignored, and was produced on the UChicago **midway** cluster (account `pi-ranganathanr`,
    partition `caslake`, remote path `/home/rouviere/LinearElastic`). `exps/download_data_midway.sh` scp's it back.
- `exps/README` is the experiment log. Summary (all use W=L=11, β=1, mismatchVar=0.1):
  - 001/002: N_types=10, N_gen=1024, μ=1e-3, τ=4 fixed. They set `S["a"]` per ens member, but the current `run.jl` doesn't read it.
  - 003: N_types=5. τ sweep from here on.
  - 004: N_gen=2048, two-bond site `(4,7),(4,16)`.
  - 005: single bond site `(4,7)`.
  - 006: k_min=1e-2.
  - 007: `selParams=(0.25/s,1)`.
  - 008: k_min=1e-2, μ=3e-3, N_reps=100, N_seqs2Keep=50.
  - 009: same as 008, but the solvent reference has Δk=0.2 instead of 1.
  - 010: same as 008, but with a single interior bond site `(47,59)`.
  - 001–008 scripts use `push!(LOAD_PATH, ...)` and the old name `evolvePopEnsemble` (now `evolvePopEns`), so they will not run
    as-is against current `src/`. 009 and 010 use `using LinearElastic` directly, so the package must be in the active
    environment (e.g. `julia --project=<repo root>` or `]dev`).
- `notebooks/`: exploratory Jupyter (IJulia) notebooks plus figures. `unit_tests/unit_tests.ipynb` is the only test suite
  (no `test/` dir). It checks the Hessian against the Python `Rigidpy` package via PyCall.
  `notes.ipynb` has the cross-experiment analysis (it hard-codes absolute Dropbox paths).
- `Rigidpy/`: a clone of github.com/VardaHagh/Rigidpy with its own `.git`. It is git-ignored and used only by the unit-test notebook.

### Gotchas / known issues (not yet fixed; confirm with user before changing)
- `updateNetParams!` writes `lpdH = lpdH` (a local no-op), so `net.lpdH` is **not** refreshed. Binding energies are
  differences in which `lpdH` cancels, so they are unaffected. Absolute free energies after `updateNetParams!` are stale.
- In `make_settings.jl`, `ΔF, s, m = computeBindingRandSeqs(...)` but the function returns `(ΔF, mean, std)`. So `s` is
  actually the mean of the centroid, and `selParams = (1/s, 1)` uses that rather than the std of the random-sequence cloud.
- `evolvePopEns` builds `k_init_vec` from two `ones(...)` blocks, so every replicate starts in environment 1. Alternating
  starts may have been intended.
- The kwarg default `N_seqs2Keep=length(nets)` in `evolvePop`/`evolvePopEns` references an undefined variable. Always pass it.
- `settings.jld2` stores `Tables`/`Ligands`/`Perturbation` structs. Changing those type definitions breaks loading old
  experiment settings with JLD2.
- `make_settings.jl` draws random `Tables` without a seed. The saved `settings.jld2` is the only record of an experiment's tables.



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
