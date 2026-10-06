# Commented-out code moved verbatim from src/ensemble.jl.
# Not included in the module; kept for reference.

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#function computeStressEns(seqs::Vector{Vector{Vector{T}}}, tables::Tables) where T <: Integer
#    stress = zeros(length(seqs[1]), length(seqs))
#    for j in axes(stress, 2), i in axes(stress, 1)
#        stress[i,j] = computeStress(seqs[j][i], tables)
#    end
#    return stress
#end

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#function computeBindingStressEns(seqs::Vector{Vector{Vector{T}}}, ligands::Ligands, tables::Tables) where T <: Integer
#    stress = zeros(length(seqs[1]), length(seqs))
#    for j in axes(stress, 2), i in axes(stress, 1)
#        stress[i,j] = computeBindingStress(seqs[j][i], ligands, tables)
#    end
#    return stress
#end

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#function computeMaxOverlapEns(seqs::Vector{Vector{Vector{T}}}, ligands::Ligands, tables::Tables) where T <: Integer
#    newseqs = hcat(seqs...)
#    return pmap( x -> computeMaxOverlap(vec(x), ligands, tables), newseqs)
#end
