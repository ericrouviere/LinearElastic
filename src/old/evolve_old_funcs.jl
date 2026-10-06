# Commented-out code moved verbatim from src/evolve.jl.
# Not included in the module; kept for reference.

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#function evolvePop(P::Integer, # population size
#                   ligands::Ligands,
#                   tables::Tables,
#                   β::Number, # Physical Temp
#                   N_gen::Int, # number of generations
#                   μ::Number, # mutation rate
#                   τ::Int, # enviromental timescale
#                   selParams::Tuple,
#                   seed::UInt; # set the seed
#                   k_init::Int=1, # Initial enviroment
#                   sampleLastHalf::Bool=true, # if true sample last half, if false return final
#                   N_seqs2Keep::Int=P) 
#    # generate population of networks and evolve it.
#    Random.seed!(seed)
#    nets = [buildRandNetwork(tables) for i in 1:P]
#    return evolvePop!(nets, ligands, tables, β, N_gen, μ, τ, selParams; k_init, sampleLastHalf, N_seqs2Keep);
#end
