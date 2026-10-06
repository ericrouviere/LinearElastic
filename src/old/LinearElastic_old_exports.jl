# Names removed from the export list in src/LinearElastic.jl, grouped under the
# section each was listed in.

# Commented out because the function is unused (its code is in src/old/*_old_funcs.jl):

# from physics.jl
#    getNetParams_old,

# from evolve.jl
#    evolvePop,

# from analysis.jl
#    scanAllostericSurface,
#    computeAllostery,

# from plotfuncs.jl
#    plotNetwork,
#    plotBindingEvoTraces,
#    plotGeneralistSpecialistEvoTraces,

# Exported but never defined anywhere in src/:

# from physics.jl
#    formatLigandParams,
#    computeBindingFreeEnergy,

# from tools.jl
#    convertSettings2NetParams,
#    convertSettings2EvoParams,

# from evolve.jl (renamed to evolvePop, the seeded single-replicate runner as in SpinModel)
#    evolvePop_seqsOnly,
