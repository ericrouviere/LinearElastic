
"""
Evolve replicate MC trajectories.
"""
function evolveMCEnsemble(N_reps::Int,
                    ligands::Ligands,
                    tables::Tables,
                    assay::Assay,
                    N_iter::Int, # Number of iterations
                    ζ::Real,
                    β::Real)
    
    function kernel(x)
        net = buildRandNetwork(tables)
        times, fitnesses, networks = evolveMC(net, ligands, tables, assay, N_iter, ζ, β; Δt=0)
        return networks[end]
    end
    return pmap( kernel, 1:N_reps)
end


function evolvePopEnsemble(seeds::Vector{UInt},
                      P::Integer,
                      ligands::Ligands,
                      tables::Tables,
                      β::Number, # Physical Temp
                      N_gen::Int, # number of generations
                      μ::Number, # mutation rate
                      τ::Int, # enviromental timescale
                      selParams::Tuple;
                      N_seqs2Keep::Int=P)
    
    # evolve many replicate populations.
    N_reps = length(seeds)
    k_init_vec = [ones(Int, Int(floor(N_reps / 2))); ones(Int, Int(ceil(N_reps / 2)))]
    seqsList= map( (x,y) -> evolvePop_seqsOnly(P, ligands, tables, β, N_gen, μ, τ, selParams, x;
                             k_init=y, N_seqs2Keep), seeds, k_init_vec)
    return seqsList
end


function computeBindingDMSEnsemble(seqsList::Vector{Vector{Vector{T}}},
                             ligands::Ligands,
                             tables::Tables,
                             β::Number) where T <: Integer
    # pokesListIndex is the index that tell what pokes to use from the
    # pokes list to measure the energy of. 
    # pokesIndex allows you choose to measure the energy of only
    # one strain. pokesIndex=0 used the defalt pokes = pokesList[pokesListIndex].
    # compute the change the energy of each applied strain for every single mutant.

    N_beads = length(seqsList[1][1])
    P = length(seqsList[1])
    N_reps = length(seqsList)
    N_ligs = length(ligands.list)-1
    N_types = size(tables.K, 1)

    DMS = Array{Float64, 5}(undef, N_types, N_beads, N_ligs, P, N_reps)
    ΔF = Array{Float64, 3}(undef, N_ligs, P, N_reps)

    net = buildNetwork(seqsList[1][1], tables)
    for i in 1:N_reps, j in 1:P
        net.seq .= seqsList[i][j] 
        updateNetParams!(net,  tables)
        F = computeFreeEnergies(net, ligands, tables, β)
        ΔF[:,j,i] = [F[2] - F[1], F[3] - F[1]]
        DMS[:,:,:,j,i] = computeBindingDMS(net, ligands, tables, β)
    end
    return DMS, ΔF
end


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


function computeStressEnsemble(seqs::Vector{Vector{T}}, tables::Tables) where T <: Integer
    return map( x -> computeStress(x, tables), seqs)
end

function computeBindingStressEnsemble(seqs::Vector{Vector{T}}, ligands::Ligands, tables::Tables) where T <: Integer
    return map( x -> computeBindingStress(x, ligands, tables), seqs)
end


function computeMaxOverlapEnsemble(seqs::Vector{Vector{T}}, ligands::Ligands, tables::Tables) where T <: Integer
    return map( x -> computeMaxOverlap(x, ligands, tables), seqs)
end


##############################################################3
##### Ensemble Helper functions ##############################
##############################################################3

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

function convertPopArray2Populations(popArray::Array{Int8, 3})
    return [ Vector{Int}[ eachcol(popArray[:,:,i])... ]  for i in 1:size(popArray, 3) ]
    #return [eachcol( popArray[:,:,i])...]
end


