
"""
    evolveMC(net, ligands, tables, assay, N_iter, ζ, β; ϕ_terminate=Inf, Δt=0)

Evolve a single network by Metropolis Monte Carlo under the selective pressure
`assay`. Each of the `N_iter` steps mutates one random bead to a new type and accepts
the mutant with probability `min(1, exp(ζ (ϕ_mut - ϕ)))`, where `ϕ` is the fitness and
`ζ` the evolutionary inverse temperature. Stops early once `ϕ > ϕ_terminate`.

Returns `(times, fitnesses, networks)` recorded every `Δt` steps, or only the final
network if `Δt = 0`. `net` itself is not modified.
"""
function evolveMC(net::Network,
                  ligands::Ligands,
                  tables::Tables,
                  assay::Assay,
                  N_iter::Int, # Number of iterations
                  ζ::Real,
                  β::Real;
                  ϕ_terminate::Number=Inf, # terminate when ϕ>ϕ_terminate
                  Δt=0)
    net = deepcopy(net)
    net_wt = deepcopy(net)
    updateNetParams!(net, tables) # make sure params are correct
    ϕ = computeFitness(net, ligands, tables, assay, β)
    times, fitnesses, networks = Int[], Float64[], Network[]
    Δt == 0 && (Δt = N_iter+1) 

    for t in 1:N_iter
        networkCopy!(net_wt, net)
        mutateBead!(net, tables)
        ϕ_mut = computeFitness(net, ligands, tables, assay, β)
        # selection
        if rand() < exp(ζ*(ϕ_mut-ϕ))
            ϕ = ϕ_mut
        else
            networkCopy!(net,net_wt)
        end
        # recording:
        if mod(t, Δt) == 0
            push!(times, t)
            push!(fitnesses, ϕ)
            push!(networks, deepcopy(net))
        end 
        ϕ > ϕ_terminate && break # terminate if fitness thresh is acheived.
    end

    # save last network
    if Δt == N_iter+1
        push!(times, N_iter)
        push!(fitnesses, ϕ)
        push!(networks, deepcopy(net))
    end 
    return times, fitnesses, networks
end

"""
    evolvePop!(nets, ligands, tables, β, evo; k_init=1, sampleLastHalf=true, N_seqs2Keep=length(nets))

Evolve the population `nets` in place for `evo.N` generations under fluctuating
selection. `ligands` holds the solvent and the two ligand states, one per environment.
Each generation:
1. scores every network in both environments, `f_k = F_solvent - F_k` (`Binding1`, `Binding2`),
2. resamples the population with weights `exp(α f_k)` in the current environment `k`,
   and weight 0 for unfolded networks, those with `F_solvent > evo.F_unfolded`
   (`isFolded`; none for `F_unfolded = Inf`); throws an error if every network is unfolded,
3. mutates every bead of every network with probability `evo.μ`,
4. switches environment every `evo.τ` generations, starting in `k_init`.

Returns `(fitBucket, netBucket, F0Bucket)`: the fitnesses of every network in both
environments at every generation (`P × 2 × N`); if `sampleLastHalf`, `N_seqs2Keep`
networks sampled at random generations of the last half of the run (aligned to whole
environmental periods), otherwise the final population; and the solvent-state free
energy of every network at every generation (`P × N`).
"""
function evolvePop!(nets::Vector{Network},
                    ligands::Ligands,
                    tables::Tables,
                    β::Number, # Physical Temp
                    evo::EvoParams;
                    k_init::Int=1, # Initial enviroment
                    sampleLastHalf::Bool=true, # if true sample last half, if false return final
                    N_seqs2Keep::Int=length(nets))
    @assert length(ligands.list) == 3 "fluctuating selection needs the solvent and two ligand states"
    return _evolvePopulation!(nets, ligands, tables, β, evo, (Binding1(), Binding2()), evo.τ;
                              k_init, sampleLastHalf, N_seqs2Keep)
end

"""
    evolvePopStatic!(nets, ligands, tables, β, evo, assay; sampleLastHalf=true, N_seqs2Keep=length(nets))

Evolve the population `nets` in place for `evo.N` generations under static selection:
every generation is scored by `computeFitness` with the same `assay` (e.g. `Binding1()`
or `DoubleBinding()`), so `evo.τ` is not used. Otherwise as `evolvePop!`; returns
`(fitBucket, netBucket, F0Bucket)` with `fitBucket` of size `P × 1 × N`, and samples
from the whole last half of the run.
"""
function evolvePopStatic!(nets::Vector{Network},
                          ligands::Ligands,
                          tables::Tables,
                          β::Number,
                          evo::EvoParams,
                          assay::Assay;
                          sampleLastHalf::Bool=true,
                          N_seqs2Keep::Int=length(nets))
    return _evolvePopulation!(nets, ligands, tables, β, evo, (assay,), 1;
                              k_init=1, sampleLastHalf, N_seqs2Keep)
end

"""
Wright-Fisher evolution of `nets` with selection on `assays[k]`, where `k` cycles through
the assays every `τ` generations starting at `k_init`; shared by `evolvePop!` (two
assays) and `evolvePopStatic!` (one assay). See `evolvePop!` for the return values.

Selection is copy-on-write: after resampling, an individual drawn several times shares
its parent's `Network` object, and gets a private copy (reusing the network of an
individual that was not drawn) only when it mutates.
"""
function _evolvePopulation!(nets::Vector{Network}, ligands::Ligands, tables::Tables,
                            β::Number, evo::EvoParams, assays::Tuple, τ::Integer;
                            k_init::Int, sampleLastHalf::Bool, N_seqs2Keep::Int)
    N_gen, μ, α, F_unfolded = evo.N, evo.μ, evo.α, evo.F_unfolded
    N_assays = length(assays)
    @assert iseven(N_gen)
    @assert 0 <= μ <= 1 # constrain mutation rate to be a probability
    @assert 0 < k_init <= N_assays # Env index has to be positive
    if sampleLastHalf
        period = N_assays * τ # one cycle through the environments
        @assert 2period <= N_gen # Env periods can't be longer that 1/4 the total evo traj.
        @assert 2N_seqs2Keep < N_gen # cant keep more seqs than half the number of generations
        gen2StartSampling = Int(N_gen - (floor((N_gen/2)/period)*period)) + 1
        gens2Sample = sort((gen2StartSampling:N_gen)[randperm(length(gen2StartSampling:N_gen))[1:N_seqs2Keep]])
    end
    P = length(nets)
    k = k_init # environment index
    fits = zeros(P, N_assays)
    F0 = zeros(P) # solvent-state free energy
    folded = trues(P) # F0 <= F_unfolded
    mutated = trues(P) # whose fitness must be recomputed
    spare = Network[] # networks of individuals that were not drawn, reused for copies
    netBucket = Network[] # list to save networks
    fitBucket = zeros(P, N_assays, N_gen)
    F0Bucket = zeros(P, N_gen)

    for t in 1:N_gen # generations
        updateAllFitnesses!(fits, F0, folded, nets, ligands, tables, assays, mutated, β, F_unfolded)
        fitBucket[:, :, t] .= fits
        F0Bucket[:, t] .= F0

        applySelection!(nets, spare, fits, F0, folded, α, k)
        mutateCopyOnWrite!(nets, spare, mutated, μ, tables)

        # switch environments
        t % τ == 0 && (k = k % N_assays + 1) # increment the env cyclically

        # Sample sequence
        if sampleLastHalf && (t in gens2Sample)
            push!(netBucket, deepcopy(nets[rand(1:P)]))
        end
    end
    _unshare!(nets, spare)
    if !sampleLastHalf
        append!(netBucket, nets)
    end
    return fitBucket, netBucket, F0Bucket
end


"""
    evolvePop(ligands, tables, β, evo, seed; k_init, sampleLastHalf, N_seqs2Keep)

Seeded variant of `evolvePop!`: build a random population of `evo.P` networks
and evolve it, returning only the sampled sequences. Seeds the global RNG with
`seed` first, so the initial population, mutations, selection and sampling are
exactly reproducible.
"""
function evolvePop(ligands::Ligands,
                   tables::Tables,
                   β::Number, # Physical Temp
                   evo::EvoParams,
                   seed::Integer;
                   k_init::Int=1, # Initial enviroment
                   sampleLastHalf::Bool=true, # if true sample last half, if false return final
                   N_seqs2Keep::Int=evo.P)
    Random.seed!(seed)
    nets = [buildRandNetwork(tables) for i in 1:evo.P]
    _, netBucket, _ = evolvePop!(nets, ligands, tables, β, evo; k_init, sampleLastHalf, N_seqs2Keep)
    return map(x -> x.seq, netBucket) # return only the sequences.
end

"""
    evolvePopStatic(ligands, tables, β, evo, assay, seed; sampleLastHalf, N_seqs2Keep)

Seeded variant of `evolvePopStatic!`: build a random population of `evo.P` networks,
evolve it under static selection on `assay`, and return only the sampled sequences.
"""
function evolvePopStatic(ligands::Ligands,
                         tables::Tables,
                         β::Number,
                         evo::EvoParams,
                         assay::Assay,
                         seed::Integer;
                         sampleLastHalf::Bool=true,
                         N_seqs2Keep::Int=evo.P)
    Random.seed!(seed)
    nets = [buildRandNetwork(tables) for i in 1:evo.P]
    _, netBucket, _ = evolvePopStatic!(nets, ligands, tables, β, evo, assay; sampleLastHalf, N_seqs2Keep)
    return map(x -> x.seq, netBucket)
end


"""
Recompute, for every network `i` with `mutated[i]`, its fitness under each assay
(`fits[i, a]`), its solvent-state free energy `F0[i]` and whether it is folded,
`folded[i] = isFolded(F0[i], F_unfolded)`. The other networks have not changed since
they were scored.
"""
function updateAllFitnesses!(fits::Matrix, F0::Vector, folded::AbstractVector{Bool},
                             nets::Vector{Network}, ligands::Ligands, tables::Tables,
                             assays::Tuple, mutated::BitVector, β::Number, F_unfolded::Real)
    for i in eachindex(nets)
        mutated[i] || continue
        F = computeFreeEnergies(nets[i], ligands, tables, β)
        for (a, assay) in enumerate(assays)
            fits[i, a] = computeFitness(assay, F)
        end
        F0[i] = F[1]
        folded[i] = isFolded(F[1], F_unfolded)
    end
    return nothing
end


"""
    fits2weights(fits, α)

Exponential map from fitness to selection weight, `w = exp(α (f - max f))`.
Larger `α` means stronger selection. Subtracting `max f` sets the largest
weight to 1, which avoids overflow without changing relative selection.
"""
function fits2weights(fits::AbstractVector, α::Number)
    return exp.(α .* (fits .- maximum(fits)))
end

"""
    fits2weights(fits, α, folded)

Selection weights `w ∝ exp(α f)` for folded individuals (`folded[i]`) and `w = 0` for
unfolded ones, normalized so that the largest weight is 1. Throws an error if no
individual is folded, since then no one can reproduce.
"""
function fits2weights(fits::AbstractVector, α::Number, folded::AbstractVector{Bool})
    any(folded) || error("every individual is unfolded (F_solvent > F_unfolded): no one can reproduce")
    z = α .* fits
    zmax = maximum(z[i] for i in eachindex(z) if folded[i])
    return ifelse.(folded, exp.(z .- zmax), zero(eltype(z)))
end

"""
Resample the population in place: draw `P` survivors with replacement, with weights
`fits2weights(fits[:, k], α, folded)` from the fitnesses in environment `k`. A
survivor drawn for the first time stays in its slot; each extra draw fills the slot of
an individual that was not drawn (in increasing slot order) and refers to the same
`Network` as its parent, so several slots may share one. Networks that are no longer
referenced are added to `spare` for reuse.
"""
function applySelection!(nets::Vector{Network}, spare::Vector{Network},
                         fits::Matrix, F0::Vector, folded::AbstractVector{Bool},
                         α::Number, # selection strength
                         k::Int) # the environment index
    P = length(nets)
    w = fits2weights(view(fits, :, k), α, folded)
    winners = wsample(1:P, w, P)
    drawn = falses(P)
    repeats = Int[] # extra draws of a survivor, in the order drawn
    for i in winners
        drawn[i] ? push!(repeats, i) : (drawn[i] = true)
    end
    holes = findall(!, drawn) # slots of individuals that were not drawn
    parents = copy(nets)
    for (hole, parent) in zip(holes, repeats)
        nets[hole] = parents[parent]
        fits[hole, :] .= view(fits, parent, :)
        F0[hole] = F0[parent]
        folded[hole] = folded[parent]
    end
    survivors = IdDict{Network, Nothing}(n => nothing for n in nets)
    released = IdDict{Network, Nothing}()
    for n in parents # a network shared by several slots is released only once
        if !haskey(survivors, n) && !haskey(released, n)
            released[n] = nothing
            push!(spare, n)
        end
    end
    return nothing
end

"""
Mutate each bead of each individual with probability `μ`, setting `mutated[i]` for the
individuals that changed. An individual whose `Network` is shared with other slots is
first given a private copy (reusing a network from `spare`), so mutations never leak
between individuals.
"""
function mutateCopyOnWrite!(nets::Vector{Network}, spare::Vector{Network},
                            mutated::BitVector, μ::Number, tables::Tables)
    refs = IdDict{Network, Int}()
    for n in nets
        refs[n] = get(refs, n, 0) + 1
    end
    for i in eachindex(nets)
        mutated[i] = false
        for j in eachindex(nets[i].seq)
            rand() < μ || continue
            if !mutated[i] && refs[nets[i]] > 1 # first mutation of a shared network
                refs[nets[i]] -= 1
                nets[i] = _privateCopy(nets[i], spare)
            end
            mutated[i] = true
            mutateBead!(nets[i], tables, j)
        end
    end
    return nothing
end

"Copy of `net` written into a network taken from `spare`, or a new one if `spare` is empty."
function _privateCopy(net::Network, spare::Vector{Network})
    isempty(spare) && return deepcopy(net)
    copy = pop!(spare)
    networkCopy!(copy, net)
    return copy
end

"Give every slot of `nets` its own `Network`, so the population can be used independently."
function _unshare!(nets::Vector{Network}, spare::Vector{Network})
    seen = IdDict{Network, Nothing}()
    for i in eachindex(nets)
        if haskey(seen, nets[i])
            nets[i] = _privateCopy(nets[i], spare)
        end
        seen[nets[i]] = nothing
    end
    return nothing
end
