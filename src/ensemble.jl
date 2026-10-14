
"""
    withPool(f, pool)

Call `f(pool)`. If `pool` is `nothing`, use a temporary `CachingPool` of all workers
instead and clear it afterwards: a `CachingPool` sends the function (and the data it
captures, such as `tables` and `ligands`) to each worker once, rather than with every
task. Without added workers, `pmap` runs everything on the main process.
"""
function withPool(f, pool)
    pool === nothing || return f(pool)
    nprocs() == 1 && return f(default_worker_pool())
    cachingPool = CachingPool(workers())
    try
        return f(cachingPool)
    finally
        clear!(cachingPool)
    end
end

"""
Evolve replicate MC trajectories, in parallel over the worker pool.
"""
function evolveMCEnsemble(N_reps::Int,
                    ligands::Ligands,
                    tables::Tables,
                    assay::Assay,
                    N_iter::Int, # Number of iterations
                    ζ::Real,
                    β::Real;
                    pool=nothing)
    
    function kernel(x)
        net = buildRandNetwork(tables)
        times, fitnesses, networks = evolveMC(net, ligands, tables, assay, N_iter, ζ, β; Δt=0)
        return networks[end]
    end
    return withPool(p -> pmap(kernel, p, 1:N_reps), pool)
end


"""
    evolvePopEnsemble(ligands, tables, β, evo, seeds; k_init_vec, N_seqs2Keep, pool=nothing)

Evolve one replicate population per entry of `seeds` under fluctuating selection
(`evolvePop`), in parallel over the worker pool (serially when no workers have been
added). Replicate `i` seeds its worker's RNG with `seeds[i]` and starts in environment
`k_init_vec[i]`, so the ensemble is exactly reproducible from the saved `seeds` and
`k_init_vec` regardless of how replicates are scheduled across workers. Each task
returns only its sampled sequences. Without a `pool`, a temporary `CachingPool` ships
`ligands` and `tables` to each worker once (`withPool`).
Returns `Vector{Vector{Vector{Int}}}` (replicate → sequence → beads).
"""
function evolvePopEnsemble(ligands::Ligands,
                           tables::Tables,
                           β::Number, # Physical Temp
                           evo::EvoParams,
                           seeds::AbstractVector{<:Integer};
                           k_init_vec=rand(1:(length(ligands.list)-1), length(seeds)),
                           N_seqs2Keep::Int=evo.P,
                           pool=nothing)
    @assert length(seeds) == length(k_init_vec)
    return withPool(pool) do p
        pmap((s, k) -> evolvePop(ligands, tables, β, evo, s; k_init=k, N_seqs2Keep),
             p, seeds, k_init_vec)
    end
end

"""
    evolvePopStaticEnsemble(ligands, tables, β, evo, assay, seeds; N_seqs2Keep, pool=nothing)

Static-selection counterpart of `evolvePopEnsemble`: one seeded replicate per entry of
`seeds`, each evolved by `evolvePopStatic` on `assay`, in parallel over the worker pool.
Returns `Vector{Vector{Vector{Int}}}` (replicate → sequence → beads).
"""
function evolvePopStaticEnsemble(ligands::Ligands,
                                 tables::Tables,
                                 β::Number,
                                 evo::EvoParams,
                                 assay::Assay,
                                 seeds::AbstractVector{<:Integer};
                                 N_seqs2Keep::Int=evo.P,
                                 pool=nothing)
    return withPool(pool) do p
        pmap(s -> evolvePopStatic(ligands, tables, β, evo, assay, s; N_seqs2Keep), p, seeds)
    end
end


"""
    computeBindingDMSEnsemble(seqsList, ligands, tables, β; pool=nothing)

Deep mutational scan of every sequence in `seqsList` (replicate → sequence →
beads), run as one flat `pmap` over all sequences. Returns
`DMS[type, bead, lig, seq, rep]` (ΔΔF of every single mutant) and the binding
energies `ΔF[lig, seq, rep]`. Tasks are sent in batches (a few per worker) because
one scan takes only milliseconds; without a `pool`, a temporary `CachingPool` ships
`ligands` and `tables` to each worker once (`withPool`).
"""
function computeBindingDMSEnsemble(seqsList::Vector{Vector{Vector{T}}},
                                   ligands::Ligands,
                                   tables::Tables,
                                   β::Number;
                                   pool=nothing) where T <: Integer

    function kernel(seq)
        net = buildNetwork(seq, tables)
        return computeBindingEnergies(net, ligands, tables, β), computeBindingDMS(net, ligands, tables, β)
    end

    N_beads = length(seqsList[1][1])
    P = length(seqsList[1])
    N_reps = length(seqsList)
    N_ligs = length(ligands.list)-1
    N_types = size(tables.K, 1)

    seqsMat = [seqsList[i][j] for j in 1:P, i in 1:N_reps]
    batch_size = max(1, length(seqsMat) ÷ (8 * nworkers()))
    out = withPool(p -> pmap(kernel, p, seqsMat; batch_size), pool)
    DMS = Array{Float64, 5}(undef, N_types, N_beads, N_ligs, P, N_reps)
    ΔF = Array{Float64, 3}(undef, N_ligs, P, N_reps)
    for i in 1:N_reps, j in 1:P
        ΔF[:,j,i], DMS[:,:,:,j,i] = out[j,i]
    end
    return DMS, ΔF
end


"""
Apply `countMutSensitivity` to every sequence of an ensemble scan
`DMS[type, bead, lig, seq, rep]`. Returns the five statistics as `P × N_reps` matrices.
"""
function countMutSensitivityEnsemble(DMS::Array{T, 5},
                                tables::Tables,
                                siteBeads::Vector{Int},
                                normalizedThresh::Number;
                                scale::Number=0.0) where T <:Number
    
    N_types, N_beads, N_ligs, P, N_reps = size(DMS)
    stepsVec = getStepsVector(tables.A, siteBeads)
    N_sensMuts_array = Matrix{Int64}(undef, P, N_reps)
    N_sensPos_array = similar(N_sensMuts_array)
    evolvability_array = Matrix{Float64}(undef, P, N_reps)
    scale_array = similar(evolvability_array)
    maxDist_array = similar(N_sensMuts_array)
    for i in 1:N_reps, j in 1:P
        @views dms = DMS[:,:,:,j,i]
        N_sensMuts, N_sensPos, evolvability, maxDist, scales =
            countMutSensitivity(dms, stepsVec, normalizedThresh; scale)
        N_sensMuts_array[j,i] = N_sensMuts
        N_sensPos_array[j,i] = N_sensPos
        evolvability_array[j,i] = evolvability
        maxDist_array[j,i] = maxDist
        scale_array[j,i] = scales
    end
    return N_sensMuts_array, N_sensPos_array, evolvability_array, maxDist_array, scale_array
end


"""
Stress `|h|` (`computeStress`) of each sequence in `seqs`.
"""
function computeStressEnsemble(seqs::Vector{Vector{T}}, tables::Tables) where T <: Integer
    return map( x -> computeStress(x, tables), seqs)
end

"""
Stress `|h|` of each sequence of each replicate (replicate → sequence → beads), as a
`P × N_reps` matrix.
"""
function computeStressEnsemble(seqs::Vector{Vector{Vector{T}}}, tables::Tables) where T <: Integer
    stress = zeros(length(seqs[1]), length(seqs))
    for j in axes(stress, 2), i in axes(stress, 1)
        stress[i,j] = computeStress(seqs[j][i], tables)
    end
    return stress
end

"""
Binding stress `|h ⋅ b|` (`computeBindingStress`) of each sequence in `seqs`.
"""
function computeBindingStressEnsemble(seqs::Vector{Vector{T}}, ligands::Ligands, tables::Tables) where T <: Integer
    return map( x -> computeBindingStress(x, ligands, tables), seqs)
end

"""
Binding stress `|h ⋅ b|` of each sequence of each replicate (replicate → sequence →
beads), as a `P × N_reps` matrix.
"""
function computeBindingStressEnsemble(seqs::Vector{Vector{Vector{T}}}, ligands::Ligands, tables::Tables) where T <: Integer
    stress = zeros(length(seqs[1]), length(seqs))
    for j in axes(stress, 2), i in axes(stress, 1)
        stress[i,j] = computeBindingStress(seqs[j][i], ligands, tables)
    end
    return stress
end


"""
Largest normal-mode overlap (`computeMaxOverlap`) of each sequence in `seqs`.
"""
function computeMaxOverlapEnsemble(seqs::Vector{Vector{T}}, ligands::Ligands, tables::Tables) where T <: Integer
    return map( x -> computeMaxOverlap(x, ligands, tables), seqs)
end

"""
Largest normal-mode overlap of each sequence of each replicate (replicate → sequence →
beads), computed in parallel with `pmap`; a `P × N_reps` matrix.
"""
function computeMaxOverlapEnsemble(seqs::Vector{Vector{Vector{T}}}, ligands::Ligands, tables::Tables;
                                   pool=nothing) where T <: Integer
    newseqs = hcat(seqs...)
    return withPool(p -> pmap(x -> computeMaxOverlap(vec(x), ligands, tables), p, newseqs), pool)
end


##############################################################3
##### Ensemble Helper functions ##############################
##############################################################3

"""
Pack populations (replicate → sequence → beads) into an `Int8` array
`(N_beads, P, N_reps)` for saving.
"""
function convertPopulations2PopArray(populations::Vector{Vector{Vector{T}}}) where T<: Integer

    N_reps = length(populations)
    P = length(populations[1])
    N_beads = length(populations[1][1])
    popArray = Array{Int8,3}(undef, N_beads, P, N_reps)
    for l in 1:N_reps
        popArray[:,:,l] = hcat(populations[l]...)
    end
    return popArray
end

"""
Unpack an `Int8` population array `(N_beads, P, N_reps)` into replicate → sequence →
beads vectors of `Int`.
"""
function convertPopArray2Populations(popArray::Array{Int8, 3})
    return [ Vector{Int}[ eachcol(popArray[:,:,i])... ]  for i in 1:size(popArray, 3) ]
    #return [eachcol( popArray[:,:,i])...]
end


