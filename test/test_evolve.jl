# Tests for fitness assays, selection and population evolution.
#
#     julia --project=. test/test_evolve.jl
#
# Pins down the Assay scoring functions, the exponential fitness map, and that
# seeded ensembles are exactly reproducible and independent of worker scheduling.

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")); io=devnull)
using Test, Distributed, LinearElastic, Random

const β = 1.0

function testSetup(; seed=1)
    Random.seed!(seed)
    tables = makeTables(7, 7, 0.15, 5, 1e-2, 1e1, 0.0, 0.1)
    i = size(tables.A, 1) ÷ 2
    j = findfirst(!iszero, tables.A[:, i])
    p(d) = Perturbation([(i, j)], [1.0], [d])
    return tables, Ligands([p(0.0), p(0.2), p(-0.2)])
end

@testset "LinearElastic evolution" begin

    @testset "assays score raw free energies" begin
        @test computeFitness(Binding(), [1.0, 0.25]) == 0.75
        @test computeFitness(DoubleBinding(), [1.0, 0.25, 0.5]) == 0.5
        @test computeFitness(Allostery(), [1.0, 2.0, 3.0, 4.5]) == 2.0 + 3.0 - 4.5 - 1.0
        # binding energies -0.25 and -0.5 against targets -0.2 and -0.4
        @test computeFitness(TargetBinding(-0.2, -0.4), [1.0, 0.75, 0.5]) ≈ -0.1
    end

    @testset "computeFitness on a network matches its free energies" begin
        tables, ligands = testSetup()
        net = buildRandNetwork(tables)
        F = computeFreeEnergies(net, ligands, tables, β)
        @test computeFitness(net, ligands, tables, DoubleBinding(), β) == min(F[1] - F[2], F[1] - F[3])
    end

    @testset "exponential fitness map" begin
        f = [0.1, -0.3, 0.05]
        w = fits2weights(f, 4.0)
        @test maximum(w) == 1
        @test w[2] / w[1] ≈ exp(4.0 * (f[2] - f[1]))
        @test all(isfinite, fits2weights([1e3, -1e3], 10.0)) # no overflow
    end

    @testset "EvoParams from a settings dict" begin
        evo = convertEvoParams(Dict("P" => 10, "N" => 64, "μ" => 1e-2, "τ" => 4, "α" => 50))
        @test (evo.P, evo.N, evo.μ, evo.τ, evo.α) == (10, 64, 0.01, 4, 50.0)
    end

    @testset "seeded ensembles are reproducible" begin
        tables, ligands = testSetup()
        evo = EvoParams(10, 64, 1e-2, 4, 50.0)
        seeds = UInt.([11, 12, 13])
        k_init_vec = [1, 2, 1]
        pops = evolvePopEnsemble(ligands, tables, β, evo, seeds; k_init_vec, N_seqs2Keep=5)
        @test length(pops) == 3
        @test all(p -> length(p) == 5, pops)
        @test pops == evolvePopEnsemble(ligands, tables, β, evo, seeds; k_init_vec, N_seqs2Keep=5)
        @test pops[2] == evolvePop(ligands, tables, β, evo, seeds[2]; k_init=2, N_seqs2Keep=5)
        # the default N_seqs2Keep is the population size
        @test length(evolvePop(ligands, tables, β, evo, UInt(1))) == evo.P

        # same populations when replicates are spread over worker processes
        addprocs(2)
        @everywhere using LinearElastic
        par = evolvePopEnsemble(ligands, tables, β, evo, seeds; k_init_vec, N_seqs2Keep=5,
                                pool=CachingPool(workers()))
        rmprocs(workers())
        @test par == pops
    end

    @testset "sampleLastHalf=false returns the final population" begin
        tables, ligands = testSetup()
        evo = EvoParams(6, 8, 1e-2, 2, 50.0)
        seqs = evolvePop(ligands, tables, β, evo, UInt(1); sampleLastHalf=false)
        @test length(seqs) == evo.P
    end

    @testset "binding DMS ensemble layout" begin
        tables, ligands = testSetup()
        seqsList = [[randSeq(tables.R, tables.K) for _ in 1:2] for _ in 1:3]
        DMS, ΔF = computeBindingDMSEnsemble(seqsList, ligands, tables, β)
        N_beads = length(tables.R) ÷ 2
        @test size(DMS) == (5, N_beads, 2, 2, 3)
        @test size(ΔF) == (2, 2, 3)
        net = buildNetwork(seqsList[3][2], tables)
        @test ΔF[:, 2, 3] ≈ computeBindingEnergies(net, ligands, tables, β)
        @test DMS[:, :, :, 2, 3] ≈ computeBindingDMS(net, ligands, tables, β)
    end

    @testset "stress ensembles take one population or a set of replicates" begin
        # the analysis notebooks broadcast over ensembles, which passes each
        # ensemble's replicate → sequence → beads list to these methods
        tables, ligands = testSetup()
        seqsList = [[randSeq(tables.R, tables.K) for _ in 1:2] for _ in 1:3]
        st = computeStressEnsemble(seqsList, tables)
        bs = computeBindingStressEnsemble(seqsList, ligands, tables)
        @test size(st) == size(bs) == (2, 3)
        @test st[:, 3] == computeStressEnsemble(seqsList[3], tables)
        @test bs[:, 3] == computeBindingStressEnsemble(seqsList[3], ligands, tables)
        @test computeMaxOverlapEnsemble(seqsList, ligands, tables)[:, 3] ≈
              computeMaxOverlapEnsemble(seqsList[3], ligands, tables)
        @test generalist([-0.1, 0.2]) == -0.2
    end
end
