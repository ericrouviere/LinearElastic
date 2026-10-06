module LinearElastic

using LinearAlgebra
using PyPlot
using Statistics
using Random
using JLD2
using Printf
using Distributed
using StatsBase: wsample
import PyPlot.PyObject
import LinearAlgebra: pinv

include("types.jl")
include("buildnetwork.jl")
include("physics.jl")
include("fitness.jl")
include("evolve.jl")
include("mutate.jl")
include("tools.jl")
include("analysis.jl")
include("plotfuncs.jl")
include("ensemble.jl")

export

    #from types.jl
    Tables,
    Network,
    Perturbation,
    Ligands,
    EvoParams,

    # from buildnetwork.jl
    buildStructure,
    buildSpringTable,
    buildMismatchTable,
    buildContactEnergyTable,
    makeTables,
    randSeq,
    buildNetwork,
    buildRandNetwork,

    # from physics.jl
    makeBeads2BondsMatrix,
    seq2NetParams,
    computeFreeEnergy,
    computeFreeEnergies,
    logpdet,
    computeBondVector,
    computeBondElongationMode,
    computePerturbedNetParams!,
    computePerturbedNetParams_old,
    pdetdiagm,
    computePerturbedGroundState,
    getNetParams,
    updateNetParams!,
    pinv,
    computeBindingRandSeqs,
    computeBindingEnergies,
    computeHessian!, # keep for unit test
    computeHessian, # keep for unit test

    # from evolve.jl
    evolveMC,
    evolvePop!,
    evolvePop,
    updateAllFitnesses!,
    fits2weights,
    applySelection!,
    splitIntoUniqesAndRepeats,

    # from fitness.jl
    computeFitness,
    Assay,
    Binding,
    DoubleBinding,
    Allostery,
    TargetBinding,

    # from mutate.jl
    mutateBead!,
    mutateAtRate!,

    # from tools.jl
    networkCopy!,
    convertEvoParams,

    # from analysis.jl
    computeBindingDMS,
    computeBindingDCS,
    countMutSensitivity,
    getStepsVector,
    getStepsMat,
    getSiteBeads,
    getSeqWithTargetBinding,
    compareSeqWithSameBinding,
    computeStress,
    computeBindingStress,
    computeConfChange,
    computeOverlaps,
    computeMaxOverlap,
    generalist,

    # from plotfuncs.jl
    plotNetwork!,
    plotDisplacement!,
    plotBondStiffness!,
    plotNetworkColoredBonds!,
    plotDensityBindingSpaceScatter,
    plotEvoPathBindingSpace,
    makePopulationFrame,
    makePopulationDynamicsMovieFrames,
    plotDMS,
    plotMutSensVsDist!,
    plotMutSensVsDist,

    # from ensemble.jl
    evolveMCEnsemble,
    evolvePopEnsemble,
    computeBindingDMSEnsemble,
    countMutSensitivityEnsemble,
    computeStressEnsemble,
    computeMaxOverlapEnsemble,
    computeBindingStressEnsemble,
    convertPopulations2PopArray,
    convertPopArray2Populations
    
end
