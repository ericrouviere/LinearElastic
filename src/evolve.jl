
function evolveMC(net::Network,
                  ligands::Ligands,
                  tables::Tables,
                  assay::Assay,
                  N_iter::Int, # Number of iterations
                  ζ::Real,
                  β::Real;
                  ϕ_terminate::Number=Inf, # terminate when ϕ>ϕ_terminate
                  Δt=0)
    # evolve network using monte carlo evolution.
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

function evolvePop!(nets::Vector{Network}, # population of networks
                    ligands::Ligands,
                    tables::Tables,
                    β::Number, # Physical Temp
                    evo::EvoParams;
                    k_init::Int=1, # Initial enviroment
                    sampleLastHalf::Bool=true, # if true sample last half, if false return final
                    N_seqs2Keep::Int=length(nets))

    # Evolve a popluation of sequences under a fluctuating selection
    N_gen, μ, τ, α = evo.N, evo.μ, evo.τ, evo.α
    N_env = length(ligands.list) - 1
    @assert N_env == 2 # this is temporary
    @assert iseven(N_gen) 
    @assert 0 <= μ <= 1 # constrain mutation rate to be a probability
    @assert 0 < k_init <= N_env # Env index has to be positive
    if sampleLastHalf
        @assert 4τ <= N_gen # Env periods can't be longer that 1/4 the total evo traj.
        @assert 2N_seqs2Keep < N_gen # cant keep more seqs than half the number of generations
        gen2StartSampling = Int(N_gen - (floor((N_gen/2)/(2τ))*2τ)) + 1
        gens2Sample = sort((gen2StartSampling:N_gen)[randperm(length(gen2StartSampling:N_gen))[1:N_seqs2Keep]])
    end
    P, N_beads = length(nets), length(nets[1].seq)
    k = k_init # environment index
    fits = zeros(length(nets), N_env)
    mutatedRecord = BitVector(undef, P); mutatedRecord .= 1 # keep track of who was mutated
    netBucket = Network[] # list to save networks
    fitBucket = zeros(P, 2, N_gen)
    
    for t in 1:N_gen # generations
        # compute Fitness
        updateAllFitnesses!(fits, nets, ligands, tables, mutatedRecord, β)
        fitBucket[:,:,t] .= fits

        # apply selection 
        applySelection!(nets, fits, α, k) 

        # mutation 
        mutatedRecord = mutateAtRate!.(nets, μ, [tables]) # mutation all sequences every generation.

        # switch environments
        t % τ == 0 && (k = k % N_env + 1) # increment the env cyclically

        # Sample sequence
        if sampleLastHalf && (t in gens2Sample)
            push!(netBucket, deepcopy(nets[rand(1:P)]))
        end
    end
    if !sampleLastHalf
        append!(netBucket, nets)
    end
    return fitBucket, netBucket
end


function evolvePop_seqsOnly(ligands::Ligands,
                   tables::Tables,
                   β::Number, # Physical Temp
                   evo::EvoParams,
                   seed::UInt; # set the seed
                   k_init::Int=1, # Initial enviroment
                   sampleLastHalf::Bool=true, # if true sample last half, if false return final
                   N_seqs2Keep::Int=evo.P) 
    # generate population of networks and evolve it
    # return only a vector or sequences
    Random.seed!(seed)
    nets = [buildRandNetwork(tables) for i in 1:evo.P]
    fitBucket, netBucket = evolvePop!(nets, ligands, tables, β, evo; k_init, sampleLastHalf, N_seqs2Keep);
    return map(x -> x.seq, netBucket) # return only the sequences. 
end


function updateAllFitnesses!(fits::Matrix,
                          nets::Vector{Network},
                          ligands::Ligands, # rename this to make more sense
                          tables::Tables,
                          mutatedRecord::BitVector,
                          β::Number)
    # Compute fitness for all sequences for enviroment all.
    # omit those sequences that have not been mutated, (0 in mutatedRecord)
    P = length(nets)
    for i in 1:P
        if mutatedRecord[i] # Bool value
            F = computeFreeEnergies(nets[i], ligands, tables, β)
            fits[i,1] = F[1] - F[2]
            fits[i,2] = F[1] - F[3] 
        end
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

function applySelection!(nets::Vector{Network},
                         fits::Matrix,
                         α::Number, # selection strength
                         k::Int) # the environment index
    # apply exponential probablity selection.
    # 1. fitnesses in env k pass through fits2weights to get a weight
    # 2. sequences are drawn for survival according to weight
    @views f = fits[:,k]
    P = length(f) 
    w = fits2weights(f, α)
    winners = wsample(1:P, w, P)
    
    uniques, repeats, holes = splitIntoUniqesAndRepeats(winners)
    for i in eachindex(repeats)
        fits[holes[i],:] = fits[repeats[i],:]
        networkCopy!(nets[holes[i]], nets[repeats[i]])
    end
    return nothing
end

function splitIntoUniqesAndRepeats(list::Vector{T}) where T
    # split list in to a list of unique elements and repeated elements
        
    uniques = eltype(list)[]
    repeats = eltype(list)[]
    #holes = Int[]
    for (i,el) in enumerate(list)
        if el in uniques
            push!(repeats, el)
        else
            push!(uniques, el)
        end
    end
    holes = setdiff(1:length(list), uniques)
    return uniques, repeats, holes
end


