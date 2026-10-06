# Commented-out code moved verbatim from src/analysis.jl.
# Not included in the module; kept for reference.

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#function scanAllostericSurface(net::Network, ligands::Ligands, tables::Tables, β::Number)

#    @assert length(ligands.list) == 3
#    # Define list of bonds
#    sep = 2
#    N_diffperts = 2
#    W = 11
#    k_allo = 1.0
#    N_beads = length(net.seq)
#    alloMismatch = 1.0
    
#    alloBonds = map( x -> (x, x-sep), N_beads-(W-4):N_beads)
#    N_alloBonds = length(alloBonds)
#    alloPerts = [ Perturbation([alloBonds[i]], [k_allo], [ifelse(isodd(j), alloMismatch, -alloMismatch)]) for i in 1:N_alloBonds, j in 1:N_diffperts]
#    net_wt, net = deepcopy(net), deepcopy(net)
#    N_ligs = length(ligands.list) - 1
#    updateNetParams!(net, tables) # compute net params from scratch.
#    ΔF_wt = computeBindingEnergies(net, ligands, tables, β)
#    alloScan = zeros(N_alloBonds, N_diffperts, N_ligs)

#    for j in 1:N_diffperts, i in 1:N_alloBonds
#        per = alloPerts[i,j]
#        computePerturbedNetParams!(net, per, tables)
#        net.G .+= net.δG
#        net.h .+= net.δh
#        net.lpdH += net.δlpdH
#        net.c += net.δc
#        F = computeBindingEnergies(net, ligands, tables, β)
#        alloScan[i, j, :] .= F .- ΔF_wt
#        networkCopy!(net, net_wt)
#    end
#    return alloScan
#end

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#function computeAllostery(net::Network, ligands::Ligands, tables::Tables, β::Number)
#   return maximum(abs, scanAllostericSurface(net, ligands, tables, β))
#end

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#computeStress(net::Network) = norm(net.h)
