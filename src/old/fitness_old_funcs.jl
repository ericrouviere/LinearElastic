# Commented-out code moved verbatim from src/fitness.jl.
# Not included in the module; kept for reference.

#function computeFitness(net::Network, ligands::Ligands, tables::Tables,
#                        assay::String, β::Number; selParams::Tuple=())
#    
#    if assay == "PottsBinding"
#        nothing
#    end
#
#    energies = computeFreeEnergies(net, ligands, tables, β)
#
#    if assay == "Binding"
#        ϕ = binding(energies...)
#    elseif assay == "Generalist"
#        ϕ = generalist(energies...)
#    elseif assay == "Cooperativity"
#        ϕ = cooperativity(energies...)
#    elseif assay == "TargetBinding"
#        isempty(selParams) && error("`selParams` not specified for TargetBiding")
#        ϕ = targetBinding( energies[2]-energies[1], energies[3]-energies[1], selParams...)
#    else
#        error("Error: Assay "*assay*" not supported.")
#    end
#    return ϕ
#end
