
"""
Random bead type different from `oldType`, uniform over the other `numTypes - 1` types.
"""
function randNewType(numTypes, oldType)
    a = rand( 1:(numTypes-1) )
    return  ((oldType-1 + a) % numTypes) +1 # new type
end

"""
    mutationPerturbation(net, tables, mutatedBead, newType)

Bonds of bead `mutatedBead` and the changes of their stiffness `Δk` and mismatch `Δd`
if the bead were mutated to `newType`, plus the change `ΔE` of their contact energy.
Returns `(bonds, Δk, Δd, ΔE)`; `net` is not modified.
"""
function mutationPerturbation(net::Network, tables::Tables, mutatedBead::Integer, newType::Integer)
    A, K, D, E = tables.A, tables.K, tables.D, tables.E
    oldType = net.seq[mutatedBead]
    connectedBeads = findall(!iszero, view(A, :, mutatedBead))
    bonds = [(mutatedBead, b) for b in connectedBeads]
    Δk = [K[newType, net.seq[b]] - K[oldType, net.seq[b]] for b in connectedBeads]
    Δd = [D[newType, net.seq[b]] - D[oldType, net.seq[b]] for b in connectedBeads]
    ΔE = sum((E[newType, net.seq[b]] - E[oldType, net.seq[b]] for b in connectedBeads); init=0.0)
    return bonds, Δk, Δd, ΔE
end

"""
Mutate bead `mutatedBead` of `net` to type `newType` and update the cached energy
parameters by a low-rank update over the bead's bonds (`applyPerturbation!`), plus
the change in the contact energy of those bonds.
"""
function mutateBead!(net::Network, tables::Tables, mutatedBead::Integer, newType::Integer)
    bonds, Δk, Δd, ΔE = mutationPerturbation(net, tables, mutatedBead, newType)
    applyPerturbation!(net, bonds, Δk, Δd, tables)
    net.seq[mutatedBead] = newType
    net.c += ΔE
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


