
struct Tables
    # Stores variables that dont change over evolution.
    R::Vector{Float64}
    A::Matrix{Float64}
    K::Matrix{Float64}
    D::Matrix{Float64}
    E::Matrix{Float64}
end


mutable struct Network
    # stores variables that change over evolution
    seq::Vector{Int}
    G::Matrix{Float64} # G = H^+
    δG::Matrix{Float64}
    lpdH::Float64 # log(pdet(H))
    δlpdH::Float64
    h::Vector{Float64}
    δh::Vector{Float64}
    c::Float64
    δc::Float64
end

struct Perturbation # 
    bonds::Vector{Tuple{Int64, Int64}}
    Δk::Vector{Float64}
    Δd::Vector{Float64}
end

struct Ligands
    list::Vector{Perturbation}
end

mutable struct EvoParams
    # Parameters of the evolutionary process.
    P::Int # population size
    N::Int # number of generations
    μ::Float64 # mutation rate per site per generation
    τ::Int # environmental timescale
    α::Float64 # selection strength
end
