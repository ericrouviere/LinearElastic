# Tests for the network physics: the fast low-rank updates used for mutations
# and ligand perturbations, checked against building the quadratic energy from
# scratch.
#
#     julia --project=. test/test_physics.jl
#
# The incremental path (Woodbury update of G = pinv(H) and the matrix
# determinant lemma for log pdet(H)) is what makes evolution and deep
# mutational scans fast, so these tests pin it to the direct computation.

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")); io=devnull)
using Test, LinearElastic, LinearAlgebra, Random, Statistics

const β = 1.0

function testTables(; seed=1)
    Random.seed!(seed)
    return makeTables(7, 7, 0.15, 5, 1e-2, 1e1, 0.5, 0.1) # nonzero contact energies
end

function testLigands(tables)
    # a single bond between two adjacent beads near the middle of the lattice
    i = size(tables.A, 1) ÷ 2
    j = findfirst(!iszero, tables.A[:, i])
    p(d) = Perturbation([(i, j)], [1.0], [d])
    return Ligands([p(0.0), p(0.2), p(-0.2)])
end

"""
Free energy of `seq` perturbed by `per`, built from scratch: the stiffness and
mismatch of each perturbed bond are changed in the diagonal bond matrices
before forming H = T'ΓT, h = T'ΓM·1 and c = ½·1'ΓM²1 + Σ ε (contact energies).
"""
function directFreeEnergy(seq, per::Perturbation, tables)
    R, A, K, D = tables.R, tables.A, tables.K, tables.D
    T = makeBeads2BondsMatrix(A, R)
    bondIndex = Dict{Tuple{Int,Int},Int}()
    k, d = Float64[], Float64[]
    for j in 1:size(A, 1)-1, i in j+1:size(A, 1)
        if A[i, j] == 1
            push!(k, K[seq[i], seq[j]])
            push!(d, D[seq[i], seq[j]])
            bondIndex[(i, j)] = bondIndex[(j, i)] = length(k)
        end
    end
    for (b, Δk, Δd) in zip(per.bonds, per.Δk, per.Δd)
        l = bondIndex[b]
        k[l] += Δk
        d[l] += Δd
    end
    H = T' * Diagonal(k) * T
    h = T' * (k .* d)
    c = sum(k .* d .^ 2) / 2 + contactEnergy(seq, tables)
    G = pinv(H)
    vals = eigvals(Symmetric(H))
    lpdH = sum(log, vals[4:end]) # 3 zero modes in 2D
    N_dims = length(h)
    return -dot(h, G, h) / 2 + c - (N_dims * log(2pi / β) - lpdH) / (2β)
end

"Sum of the contact energies `E[s_i, s_j]` over the bonds of the network."
function contactEnergy(seq, tables)
    A, E = tables.A, tables.E
    return sum(E[seq[i], seq[j]] for j in 1:size(A, 1)-1 for i in j+1:size(A, 1) if A[i, j] == 1)
end

"Copy of `tables` with the contact-energy table replaced by `E`."
withContactEnergies(tables, E) = Tables(tables.R, tables.A, tables.K, tables.D, E)

@testset "LinearElastic physics" begin

    @testset "mutations: incremental update matches a fresh build" begin
        tables = testTables()
        net = buildRandNetwork(tables)
        Random.seed!(2)
        for _ in 1:30
            mutateBead!(net, tables)
        end
        G, lpdH, h, c = getNetParams(net.seq, tables)
        @test net.G ≈ G rtol=1e-8
        @test net.lpdH ≈ lpdH rtol=1e-10
        @test net.h ≈ h rtol=1e-10
        @test net.c ≈ c rtol=1e-10
    end

    @testset "updateNetParams! refreshes every cached parameter" begin
        tables = testTables()
        net = buildRandNetwork(tables)
        net.seq .= randSeq(tables.R, tables.K)
        updateNetParams!(net, tables)
        G, lpdH, h, c = getNetParams(net.seq, tables)
        @test net.G == G
        @test net.lpdH == lpdH
        @test net.h == h
        @test net.c == c
    end

    @testset "ligand perturbations: low-rank free energy matches direct" begin
        tables = testTables()
        ligands = testLigands(tables)
        net = buildRandNetwork(tables)
        for per in ligands.list
            @test computeFreeEnergy(net, per, tables, β) ≈ directFreeEnergy(net.seq, per, tables) rtol=1e-9
        end
        # a perturbation that both stiffens and stretches two bonds at once
        b1, b2 = ligands.list[1].bonds[1], (ligands.list[1].bonds[1][1], findlast(!iszero, tables.A[:, ligands.list[1].bonds[1][1]]))
        per = Perturbation([b1, b2], [0.5, 2.0], [0.1, -0.3])
        @test computeFreeEnergy(net, per, tables, β) ≈ directFreeEnergy(net.seq, per, tables) rtol=1e-9
    end

    @testset "deep mutational scan matches rebuilding each mutant" begin
        tables = testTables()
        ligands = testLigands(tables)
        net = buildRandNetwork(tables)
        dms = computeBindingDMS(net, ligands, tables, β)
        ΔF_wt = computeBindingEnergies(buildNetwork(net.seq, tables), ligands, tables, β)
        N_types = size(tables.K, 1)
        for j in eachindex(net.seq), t in 1:N_types
            if t == net.seq[j]
                @test all(iszero, dms[t, j, :])
                continue
            end
            seq = copy(net.seq)
            seq[j] = t
            ΔF = computeBindingEnergies(buildNetwork(seq, tables), ligands, tables, β)
            @test dms[t, j, :] ≈ ΔF .- ΔF_wt atol=1e-9
        end
    end

    @testset "contact energies: any symmetric table E" begin
        tables = testTables()
        ligands = testLigands(tables)
        q = size(tables.K, 1)
        E = [Float64(a * b) - 3.0 * (a == b) for a in 1:q, b in 1:q] # chosen, not random
        tablesE = withContactEnergies(tables, E)
        tables0 = withContactEnergies(tables, zeros(q, q))
        net = buildRandNetwork(tablesE)
        seq = copy(net.seq)

        # the energy constant c includes the contact energy of every bond
        _, _, c = seq2NetParams(seq, tablesE)
        _, _, c0 = seq2NetParams(seq, tables0)
        @test c ≈ c0 + contactEnergy(seq, tablesE)
        # so the free energy shifts by exactly the contact energy
        @test computeFreeEnergy(buildNetwork(seq, tablesE), β) ≈
              computeFreeEnergy(buildNetwork(seq, tables0), β) + contactEnergy(seq, tablesE)

        # mutations update the contact energy incrementally
        Random.seed!(5)
        for _ in 1:30
            mutateBead!(net, tablesE)
        end
        @test net.c ≈ getNetParams(net.seq, tablesE)[4] rtol=1e-10

        # binding does not change the sequence, so contact energies cancel in ΔF
        @test computeBindingEnergies(buildNetwork(seq, tablesE), ligands, tablesE, β) ≈
              computeBindingEnergies(buildNetwork(seq, tables0), ligands, tables0, β)
        @test computeBindingDMS(buildNetwork(seq, tablesE), ligands, tablesE, β) ≈
              computeBindingDMS(buildNetwork(seq, tables0), ligands, tables0, β)

        # a non-symmetric table would make the full build and the mutation update disagree
        Easym = copy(E); Easym[1, 2] += 1
        @test_throws ArgumentError buildNetwork(seq, withContactEnergies(tables, Easym))
        @test_throws ArgumentError buildNetwork(seq, withContactEnergies(tables, zeros(q, q + 1)))
    end

    @testset "computeBindingRandSeqs returns (ΔF, std, mean)" begin
        tables = testTables()
        ligands = testLigands(tables)
        Random.seed!(3)
        ΔF, s, m = computeBindingRandSeqs(200, ligands, tables, β)
        centroid = mean(ΔF, dims=1)
        @test s ≈ std(sqrt.(sum((ΔF .- centroid) .^ 2, dims=2)))
        @test m ≈ mean(centroid)
    end
end
