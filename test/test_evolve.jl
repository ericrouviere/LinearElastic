# Tests for fitness assays, selection and population evolution.
#
#     julia --project=. test/test_evolve.jl
#
# Pins down the Assay scoring functions, the exponential fitness map, and that
# seeded ensembles are exactly reproducible and independent of worker scheduling.

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")); io=devnull)
using Test, Distributed, LinearElastic, Random, Statistics

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
        @test computeFitness(Binding1(), [1.0, 0.25, 0.5]) == 0.75
        @test computeFitness(Binding2(), [1.0, 0.25, 0.5]) == 0.5
    end

    @testset "folded fraction" begin
        @test foldedFraction(-10.0, 0.0, 1.0) ≈ 1 / (1 + exp(-10))
        @test foldedFraction(0.0, 0.0, 1.0) == 0.5
        @test foldedFraction(10.0, 0.0, 2.0) ≈ 1 / (1 + exp(20))
        @test foldedFraction(-120.0, Inf, 1.0) == 1 # no unfolding by default
        @test logFoldedFraction(1e4, 0.0, 1.0) ≈ -1e4 # no overflow far above the threshold
        @test logFoldedFraction(-1e4, 0.0, 1.0) == 0
        f, lnp = [0.1, -0.3, 0.05], [0.0, -2.0, -50.0]
        w = fits2weights(f, 4.0, lnp)
        @test maximum(w) == 1
        @test w[3] / w[1] ≈ exp(4.0 * (f[3] - f[1]) + lnp[3] - lnp[1])
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
        # without a pool, a temporary CachingPool of all workers is used
        parDefault = evolvePopEnsemble(ligands, tables, β, evo, seeds; k_init_vec, N_seqs2Keep=5)
        DMS_par, ΔF_par = computeBindingDMSEnsemble(pops, ligands, tables, β)
        rmprocs(workers())
        @test par == pops
        @test parDefault == pops
        DMS_ser, ΔF_ser = computeBindingDMSEnsemble(pops, ligands, tables, β)
        @test DMS_par ≈ DMS_ser && ΔF_par ≈ ΔF_ser
    end

    @testset "copy-on-write selection keeps every network consistent" begin
        tables, ligands = testSetup()
        evo = EvoParams(30, 200, 2e-2, 10, 200.0)
        nets = [buildRandNetwork(tables) for _ in 1:evo.P]
        F0_initial = [computeFreeEnergies(n, ligands, tables, β)[1] for n in nets]
        _, sampled, F0Bucket = evolvePop!(nets, ligands, tables, β, evo; N_seqs2Keep=10)
        @test length(IdSet{Network}(nets)) == evo.P # no two individuals share a network
        for n in [nets; sampled]
            fresh = buildNetwork(copy(n.seq), tables)
            @test n.G ≈ fresh.G rtol=1e-6
            @test n.u ≈ fresh.u rtol=1e-6
            @test n.lpdH ≈ fresh.lpdH
        end
        # F0Bucket[:, t] is the population at the start of generation t
        @test size(F0Bucket) == (evo.P, evo.N)
        @test F0Bucket[:, 1] == F0_initial
    end

    @testset "static selection" begin
        tables, ligands = testSetup()
        evo = EvoParams(40, 200, 1e-2, 1, 300.0)
        nets = [buildRandNetwork(tables) for _ in 1:evo.P]
        fitBucket, sampled, _ = evolvePopStatic!(nets, ligands, tables, β, evo, Binding1(); N_seqs2Keep=10)
        @test size(fitBucket) == (evo.P, 1, evo.N)
        @test length(sampled) == 10
        # selection on ligand 1 improves binding of ligand 1
        @test mean(fitBucket[:, 1, end]) > mean(fitBucket[:, 1, 1])
        seeds = UInt.([5, 6])
        pops = evolvePopStaticEnsemble(ligands, tables, β, evo, Binding1(), seeds; N_seqs2Keep=5)
        @test pops[2] == evolvePopStatic(ligands, tables, β, evo, Binding1(), seeds[2]; N_seqs2Keep=5)
    end

    @testset "unfolding threshold" begin
        tables, ligands = testSetup()
        Random.seed!(7)
        F0_rand = [computeFreeEnergies(buildRandNetwork(tables), ligands, tables, β)[1] for _ in 1:200]
        F_unfolded = sort(F0_rand)[20] # only about 10% of random sequences fold
        evo = EvoParams(40, 100, 1e-2, 1, 1.0) # weak selection on binding
        foldedAtEnd(F_unf) = begin
            Random.seed!(3)
            nets = [buildRandNetwork(tables) for _ in 1:evo.P]
            _, _, F0Bucket = evolvePopStatic!(nets, ligands, tables, β, evo, Binding1();
                                              sampleLastHalf=false, F_unfolded=F_unf)
            mean(foldedFraction.(F0Bucket[:, end], F_unfolded, β))
        end
        @test foldedAtEnd(F_unfolded) > 0.9 # selection drives the population to fold
        @test foldedAtEnd(Inf) < 0.5 # without the threshold it does not
        # the default F_unfolded = Inf is exactly the evolution without unfolding
        @test evolvePop(ligands, tables, β, EvoParams(10, 64, 1e-2, 4, 50.0), UInt(1)) ==
              evolvePop(ligands, tables, β, EvoParams(10, 64, 1e-2, 4, 50.0), UInt(1); F_unfolded=Inf)
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
