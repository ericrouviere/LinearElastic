# Test the network Hessian against the Python package Rigidpy
# (github.com/VardaHagh/Rigidpy), as in notebooks/unit_tests/unit_tests.ipynb.
#
#     julia --project=. test/test_rigidpy.jl
#
# Rigidpy is read from the git-ignored clone in Rigidpy/ at the repo root (it does
# not need to be pip-installed); the tests are skipped if the clone is missing.
# Python and numpy/scipy come from PyCall's Python, which PyPlot already requires.

using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(REPO_ROOT; io=devnull)
using Test, LinearElastic, LinearAlgebra, Random
import PyPlot
const PyCall = PyPlot.PyCall

const RIGIDPY_DIR = joinpath(REPO_ROOT, "Rigidpy")

"""
Import `rigidpy` from the clone at `RIGIDPY_DIR`, or return `nothing` if it is missing.
"""
function importRigidpy()
    isdir(joinpath(RIGIDPY_DIR, "rigidpy")) || return nothing
    sys = PyCall.pyimport("sys")
    RIGIDPY_DIR in sys."path" || pushfirst!(PyCall.PyVector(sys."path"), RIGIDPY_DIR)
    return PyCall.pyimport("rigidpy")
end

"""
Hessian of the network `seq` computed by Rigidpy: the bead coordinates, the
bonds of `tables.A` and each bond's stiffness `K[seq[i], seq[j]]`. Bonds are
listed in the same order as in `makeBeads2BondsMatrix`; Rigidpy is 0-indexed.
"""
function rigidpyHessian(rp, seq, tables)
    A, N = tables.A, size(tables.A, 1)
    coordinates = permutedims(reshape(tables.R, 2, :)) # (N, 2): one row per bead
    bonds = [(i, j) for j in 1:N-1 for i in j+1:N if A[i, j] == 1]
    bondMat = [b[m] - 1 for b in bonds, m in 1:2]
    k = [tables.K[seq[i], seq[j]] for (i, j) in bonds]
    framework = rp.framework(coordinates, bondMat; k)
    return framework.hessianMatrix(), framework.rigidityMatrix(), k
end

rp = importRigidpy()

@testset "Hessian vs Rigidpy" begin
    if rp === nothing
        @warn "Rigidpy clone not found at $RIGIDPY_DIR; skipping the Rigidpy comparison."
        @test_skip false
    else
        @testset "random stiffnesses, seed $seed" for seed in 1:3
            Random.seed!(seed)
            tables = makeTables(11, 11, 0.15, 5, 1e-2, 1e1, 0.0, 0.1)
            seq = randSeq(tables.R, tables.K)
            H, _ = seq2NetParams(seq, tables)
            H_rp, R_rp, k = rigidpyHessian(rp, seq, tables)
            @test size(H_rp) == size(H)
            @test H ≈ H_rp
            # the beads-to-bonds matrix is Rigidpy's rigidity matrix, up to the
            # sign convention of each bond vector, so H = T'ΓT = R'ΓR
            T = makeBeads2BondsMatrix(tables.A, tables.R)
            @test abs.(T) ≈ abs.(R_rp)
            @test H ≈ R_rp' * Diagonal(k) * R_rp
        end

        @testset "unit stiffness (the notebook's comparison)" begin
            Random.seed!(4)
            tables = makeTables(11, 11, 0.1, 5, 1.0, 1.0, 0.0, 0.0)
            seq = randSeq(tables.R, tables.K)
            H, _ = seq2NetParams(seq, tables)
            _, R_rp, _ = rigidpyHessian(rp, seq, tables)
            @test H ≈ R_rp' * R_rp
        end
    end
end
