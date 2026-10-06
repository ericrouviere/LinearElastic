
"""
Quantities that are fixed over evolution: bead positions `R = [x1, y1, x2, y2, ...]`,
adjacency matrix `A`, and the symmetric bead-type tables `K` (bond stiffness), `D`
(bond rest-length mismatch) and `E` (bond contact energy): a bond between beads of
types `a` and `b` has stiffness `K[a, b]`, mismatch `D[a, b]` and contact energy `E[a, b]`.
"""
struct Tables
    R::Vector{Float64}
    A::Matrix{Float64}
    K::Matrix{Float64}
    D::Matrix{Float64}
    E::Matrix{Float64}
end


"""
One evolving network: its sequence `seq` and the cached parameters of its quadratic
energy `U(r) = ½ r'Hr - h'r + c`, namely `G = pinv(H)`, `lpdH = log pdet(H)`, `h` and
`c`, together with the ground-state displacement `u = G h` and `hGh = h'G h`, which
make perturbed free energies cheap (`perturbedFreeEnergy`). The fields `δG`, `δlpdH`,
`δh` and `δc` are scratch space for `computePerturbedNetParams!`.
"""
mutable struct Network
    seq::Vector{Int}
    G::Matrix{Float64} # G = H^+
    δG::Matrix{Float64}
    lpdH::Float64 # log(pdet(H))
    δlpdH::Float64
    h::Vector{Float64}
    δh::Vector{Float64}
    c::Float64
    δc::Float64
    u::Vector{Float64} # ground-state displacement G h
    hGh::Float64 # h'G h
end

"""
Change in stiffness `Δk` and rest-length mismatch `Δd` of each bond in `bonds` (pairs
of bead indices). Models ligand binding, and is also used to apply mutations.
"""
struct Perturbation # 
    bonds::Vector{Tuple{Int64, Int64}}
    Δk::Vector{Float64}
    Δd::Vector{Float64}
end

"""
Ligand states: `list[1]` is the solvent (unbound reference) and `list[2:end]` are the
bound states.
"""
struct Ligands
    list::Vector{Perturbation}
end

"""
Parameters of the evolutionary process: population size `P`, number of generations
`N`, mutation rate per site per generation `μ`, environmental timescale `τ`
(generations between environment switches) and selection strength `α`.
"""
mutable struct EvoParams
    P::Int # population size
    N::Int # number of generations
    μ::Float64 # mutation rate per site per generation
    τ::Int # environmental timescale
    α::Float64 # selection strength
end
