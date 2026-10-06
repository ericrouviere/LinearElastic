
function randNewType(numTypes, oldType)
    a = rand( 1:(numTypes-1) )
    return  ((oldType-1 + a) % numTypes) +1 # new type
end

function mutateBead!(net::Network, tables::Tables, mutatedBead::Int, newType::Integer)
    # Mutate sequence at position `mutatedBead` with a bead of newType.
    # update the quadratic energy params.
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

function mutateBead!(net::Network, tables::Tables, mutatedBead::Int)
    # Mutate sequence at position `mutatedBead` with random new bead type.
    # update the quadratic energy params.
    N_types = size(tables.K, 1)
    newType = randNewType(N_types, net.seq[mutatedBead]) # mutate sequence
    mutateBead!(net, tables, mutatedBead, newType)
    return nothing
end

function mutateBead!(net::Network, tables::Tables)
    # Mutate sequence at random position , update the quadratic energy params
    N_beads = length(net.seq)
    mutateBead!(net, tables, rand(1:N_beads))
    return nothing
end


function mutateAtRate!(net::Network, μ::Number, tables::Tables)
    # Mutate each position of sequence with probability μ.
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


