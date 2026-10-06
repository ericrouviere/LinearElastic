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

# String-based assay dispatch, replaced by the Assay types in src/fitness.jl.
# generalist(ΔF::Vector) lives on in src/analysis.jl.
#binding(F0::Number, F1::Number) = F0-F1 # strong binders > 0
#generalist(F0::Number, F1::Number, F2::Number) = min(binding(F0,F1), binding(F0,F2))
#generalist(ΔF::Vector) = min( -ΔF[1], -ΔF[2]) # here ΔF is the binding energy
#cooperativity(E00::Number, E10::Number, E01::Number, E11::Number) = E10 + E01 - E11 - E00

#function targetBinding(ΔF1::Number, ΔF2::Number, ΔF1_target::Number, ΔF2_target::Number)
#    # return a fitness function that increased (to zero) as
#    # the binding energies approach the target binding energies.
#    err1 = abs(ΔF1 - ΔF1_target)
#    err2 = abs(ΔF2 - ΔF2_target)
#    return -max(err1, err2)
#end


#function computeFitness(assay::String, energies::Vector; selParams::Tuple=())
#    # compute fitness from energies
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

#function computeFitness(net::Network, ligands::Ligands, tables::Tables,
#                        assay::String, β::Number; selParams::Tuple=())
#    return computeFitness(assay, computeFreeEnergies(net, ligands, tables, β); selParams)
#end
