
"""
Random bead type different from `oldType`, uniform over the other `numTypes - 1` types.
"""
function randNewType(numTypes, oldType)
    a = rand( 1:(numTypes-1) )
    return  ((oldType-1 + a) % numTypes) +1 # new type
end

"""
Mutate bead `mutatedBead` of `net` to type `newType` and update the cached energy
parameters by a low-rank update over the bead's bonds (`computePerturbedNetParams!`).
"""
function mutateBead!(net::Network, tables::Tables, mutatedBead::Int, newType::Integer)
    A, R, K, D, E= tables.A, tables.R, tables.K, tables.D, tables.E
    N_types = size(tables.K, 1)
    oldType = net.seq[mutatedBead]
    #newType = randNewType(N_types, seq[mutatedBead]) # mutate sequence
    connectedBeads = findall(x -> x!=0, A[:, mutatedBead]) # find beads connected to the mutated bead
    bonds = [(mutatedBead, conBead) for conBead in connectedBeads]
    N_m = length(connectedBeads)
    Δk = zeros(N_m)
    Δd = zeros(N_m)
    for (i, bead_i) in enumerate(connectedBeads)
        beadType_i = net.seq[bead_i]
        Δk[i] = K[newType, beadType_i] - K[oldType, beadType_i]
        Δd[i] = D[newType, beadType_i] - D[oldType, beadType_i]
    end
    mutation = Perturbation(bonds, Δk, Δd)
    computePerturbedNetParams!(net, mutation, tables)
    net.seq[mutatedBead] = newType
    net.G .+= net.δG
    net.h .+= net.δh
    net.lpdH += net.δlpdH
    net.c += net.δc
    return nothing
end

"""
Mutate bead `mutatedBead` of `net` to a random different type.
"""
function mutateBead!(net::Network, tables::Tables, mutatedBead::Int)
    N_types = size(tables.K, 1)
    newType = randNewType(N_types, net.seq[mutatedBead]) # mutate sequence
    mutateBead!(net, tables, mutatedBead, newType)
    return nothing
end

"""
Mutate a random bead of `net` to a random different type.
"""
function mutateBead!(net::Network, tables::Tables)
    N_beads = length(net.seq)
    mutateBead!(net, tables, rand(1:N_beads))
    return nothing
end


"""
Mutate each bead of `net` independently with probability `μ`. Returns `true` if any
bead was mutated.
"""
function mutateAtRate!(net::Network, μ::Number, tables::Tables)
    didMutate = false
    N_types = size(tables.K, 1)
    for i in eachindex(net.seq)
        if rand() < μ
            didMutate = true
            mutateBead!(net, tables, i)
        end
    end
    return didMutate
end


