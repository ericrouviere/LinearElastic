

"""
    computeBindingDMS(net, ligands, tables, β)

Deep mutational scan of the binding energies of `net`. `DMS[t, j, k]` is the change
`ΔΔF` in the binding energy to ligand `k` when bead `j` is mutated to type `t` (zero
for the wild-type type). Returns an `N_types × N_beads × N_ligands` array; `net` is
not modified.
"""
function computeBindingDMS(net::Network, ligands::Ligands, tables::Tables, β::Number)
    net = buildNetwork(copy(net.seq), tables) # parameters from scratch
    ΔF_wt = computeBindingEnergies(net, ligands, tables, β)
    N_types = size(tables.K, 1)
    N_beads = length(net.seq)
    N_ligs = length(ligands.list)
    DMS = zeros(N_types, N_beads, N_ligs-1)
    F = zeros(N_ligs)
    for j in 1:N_beads, type_mut in 1:N_types
        type_mut == net.seq[j] && continue
        # each mutant in each ligand state is one low-rank perturbation of the wild type
        mbonds, mΔk, mΔd, _ = mutationPerturbation(net, tables, j, type_mut)
        for (s, per) in enumerate(ligands.list)
            bonds, Δk, Δd = mergePerturbations(mbonds, mΔk, mΔd, per.bonds, per.Δk, per.Δd)
            F[s] = perturbedFreeEnergy(net, bonds, Δk, Δd, tables, β)
        end
        for l in 1:N_ligs-1
            DMS[type_mut, j, l] = (F[l+1] - F[1]) - ΔF_wt[l]
        end
    end
    return DMS # an array of 3 dims
end

"""
Double mutational scan of the binding energies of `net`. Returns `S, D`: the
single-mutant `ΔΔF` (as `computeBindingDMS`) and `D[t1, j1, t2, j2, k]`, the `ΔΔF` to
ligand `k` of the double mutant (bead `j1` to type `t1`, bead `j2` to type `t2`),
symmetric in the two mutations.
"""
function computeBindingDCS(net::Network, ligands::Ligands, tables::Tables, β::Number)
    net = deepcopy(net)
    net_wt = deepcopy(net)
    updateNetParams!(net, tables)
    ΔF_wt = computeBindingEnergies(net, ligands, tables, β)

    N_types, N_beads, N_ligs = size(tables.K, 1), length(net.seq), length(ligands.list)
    S = zeros(N_types, N_beads, N_ligs-1)
    D = zeros(N_types, N_beads, N_types, N_beads, N_ligs-1)

    # Double mutants (symmetric assignment, includes singles)
    for j1 in 1:N_beads
        wt1 = net.seq[j1]
        for t1 in 1:N_types
            mutateBead!(net, tables, j1, t1)

            # Compute and store single mutant effect only if not WT
            if t1 != wt1
                ΔF1 = computeBindingEnergies(net, ligands, tables, β)
                ΔΔF1 = ΔF1 .- ΔF_wt
                S[t1, j1, :] = ΔΔF1
            end

            for j2 in j1+1:N_beads
                wt2 = net.seq[j2]
                for t2 in 1:N_types
                    mutateBead!(net, tables, j2, t2)
                    ΔF = computeBindingEnergies(net, ligands, tables, β)
                    ΔΔF = ΔF .- ΔF_wt
                    D[t1, j1, t2, j2, :] = ΔΔF
                    D[t2, j2, t1, j1, :] = ΔΔF
                    mutateBead!(net, tables, j2, wt2)
                end
            end
            networkCopy!(net, net_wt)
        end
    end

    return S, D
end


"""
    countMutSensitivity(dms, stepsVec, normalizedThresh; scale=0.0)

Mutational-sensitivity statistics of one deep mutational scan `dms`
(`N_types × N_beads × N_ligands`). A mutation is sensitive if the norm of its `ΔΔF`
over the ligands exceeds `normalizedThresh * scale`; `scale = 0` uses the mean `|ΔΔF|`
of mutations within one bond of the binding site. `stepsVec[j]` is the number of
bonds between bead `j` and the binding site.

Returns `(N_sensMuts, N_sensPos, evolvability, maxDist, scale)`: the number of
sensitive mutations, the number of positions with at least one, the root-mean-square
norm of `ΔΔF` over all single mutants (the evolvability), the largest distance of a
sensitive mutation from the binding site, and the scale used.
"""
function countMutSensitivity(dms::AbstractArray{T,3},
                             stepsVec::Vector,
                             normalizedThresh::Number; # sensitivity threshold relative to `scale`. eg 10%.
                             scale::Number=0.0 # use 0.0 for relative scale
                             ) where T <: Number
    N_types, N_beads, N_ligs = size(dms)

    if scale==0.0 # use active site mutations for scale
        scale = mean(abs, dms[:, 0 .<= stepsVec .<= 1,:])
    end

    thresh = normalizedThresh * scale
    N_sensMuts = N_sensPos = maxDist = 0
    acc = 0.0
    for j in axes(dms,2) # length of sequence
        pos = false
        dist = stepsVec[j]
        for i in axes(dms,1) # loop over number of bead types
            Δϕ = norm(dms[i,j,:])
            Δϕ2 = Δϕ^2
           # Δϕ2 = dms[i,j,1]^2 + dms[i,j,2]^2
            acc += Δϕ2
            if Δϕ > thresh
                N_sensMuts += 1
                pos = true
                dist > maxDist && (maxDist = dist)
            end
        end
        pos && (N_sensPos +=1) 
    end
    N_muts = (N_types-1)*N_beads
    evolvability = sqrt(acc/N_muts)
    return N_sensMuts, N_sensPos, evolvability, maxDist, scale
end


"""
Same as `countMutSensitivity(dms, stepsVec, normalizedThresh)`, with the distances
computed from the binding-site beads `siteBeads` and the adjacency matrix `adjMat`.
"""
function countMutSensitivity(dms::AbstractArray{T,3},
                             siteBeads::Vector,
                             adjMat::Matrix,
                             normalizedThresh::Number; # sensitivity threshold relative to `scale` eg 10%.
                             scale::Number=0.0 # std of ΔF cloud of random seqs
                            ) where T <: Number
    stepsVec = getStepsVector(adjMat, siteBeads)
    return countMutSensitivity(dms, stepsVec, normalizedThresh; scale)
end


"""
Matrix of graph distances (number of bonds) between every pair of beads of the network
with adjacency matrix `A`.
"""
function getStepsMat(A::Matrix{T}) where T<:Number
    full, i = false, 0
    B = similar(A)
    bucket = BitMatrix[]
    while !full
        B = A^i .!= 0
        push!(bucket, copy(B))
        full = sum( .!B ) == 0
        i += 1
    end
    C = reduce(+, bucket)
    C = maximum(C) .+ 1 .- C 
    C = C .- I(size(A,1))
    return C
end

"""
Graph distance from each bead to the nearest bead in `siteBeads` (the binding site).
"""
function getStepsVector(A::Matrix, siteBeads::AbstractVector)
    stepsMat = getStepsMat(A)
    X = stepsMat[:,siteBeads]
    return vec(minimum(X, dims=2))
end

"""
Beads of the binding site: the beads of the bonds perturbed by the first ligand state.
All states are assumed to act on the same bonds.
"""
function getSiteBeads(ligands::Ligands)
    return unique(vcat(collect.(ligands.list[1].bonds)...))
end


##################################################################


"""
Return networks with binding energies similar to `ΔF_target`.
"""
function getSeqWithTargetBinding(ΔF_target::Vector,
                                 ligands::Ligands,
                                 tables::Tables,
                                 β::Number,
                                 similarityThresh::Number) # binding energy scale

    assay = TargetBinding(ΔF_target...)
    N_iter = 5000
    ζ = 1e4
    ϕ_terminate = -similarityThresh
    net_similar = buildRandNetwork(tables)
    times, fitnesses, networks = evolveMC(net_similar,
                                          ligands,
                                          tables,
                                          assay,
                                          N_iter,
                                          ζ,
                                          β;
                                          Δt=0,
                                          ϕ_terminate);
    if fitnesses[end] < ϕ_terminate
        println("Did not find sequence with desired function")
        return false, networks[end]
    else 
        return true, networks[end]
    end
end


"""
Compare allostery and evolvablity scores with net with same function. 
"""
function compareSeqWithSameBinding(ΔF_target::Vector,
                                   dms::Array{T, 3},
                                   ligands::Ligands,
                                   tables::Tables,
                                   β::Number;
                                   similarityThresh::Number=1,
                                   actSiteBeads::Vector=[1],
                                   mutscale::Number=0.01,
                                   normalizedThresh::Number=1.0
                                  ) where T # binding energy scale

    # THIS IS REALLY GROSS CODE
    bool, net_ref = getSeqWithTargetBinding(ΔF_target, ligands, tables, β, similarityThresh)
    if bool==false
        return false, (0.0, 0.0, 0.0, 0.0)
    end

    dms_ref = computeBindingDMS(net_ref, ligands, tables, β)
    nm, np, e, md, _ = countMutSensitivity(dms, actSiteBeads, tables.A, normalizedThresh; scale=mutscale)
    nmr, npr, er, mdr, _ = countMutSensitivity(dms_ref, actSiteBeads, tables.A, normalizedThresh; scale=mutscale)
    answers = ( (nm-nmr)/nmr, (np-npr)/npr, (e-er)/er, (md-mdr)/mdr )
    return bool, answers
end

"""
Run `compareSeqWithSameBinding` for each pair of target binding energies and DMS, and
return the relative differences `(dnm, dnp, de, dmd)` of the comparisons that found a
reference sequence.
"""
function compareSeqWithSameBinding(ΔF_target_vec::Vector{Vector{T}},
                                   dms_vec::Vector{Array{T, 3}},
                                   ligands::Ligands,
                                   tables::Tables,
                                   β::Number;
                                   similarityThresh::Number=0.01,
                                   actSiteBeads::Vector=[1],
                                   mutscale::Number=0.01,
                                   normalizedThresh::Number=1.0
                                  ) where T # binding energy scale

    @assert length(ΔF_target_vec) == length(dms_vec)
    dnm = Float64[]
    dnp = Float64[]
    de = Float64[]
    dmd = Float64[]

    for i in eachindex(dms_vec)
        bool, answers = compareSeqWithSameBinding(ΔF_target_vec[i], dms_vec[i],
                          ligands, tables, β; similarityThresh, actSiteBeads,
                          mutscale, normalizedThresh)
        if bool
            push!(dnm, answers[1])
            push!(dnp, answers[2])
            push!(de,  answers[3])
            push!(dmd, answers[4])
        end
    end    
    return dnm, dnp, de, dmd
end


"""
Compute the norm of the vector h, which describes the 
internal stress in a network.
"""
function computeStress(seq::Vector, tables::Tables)
    H, h, c = seq2NetParams(seq, tables)
    return norm(h)
end


"""
Compute the dot product of h and b, the network stress vector
and the binding mode, respectively.
"""
function computeBindingStress(seq::Vector, ligands::Ligands, tables::Tables)

    function computeBindingMode(per::Perturbation, tables::Tables)
        R, bonds = tables.R, per.bonds
        N_m = length(bonds)
        T_m = Matrix{Float64}(undef, N_m, length(R)) ###################
        for i in 1:N_m
            bead_1, bead_2 = bonds[i]
            T_m[i,:] = computeBondElongationMode(R, bead_1,  bead_2);
        end   
        return vec(sum(T_m, dims=1))
    end
    H, h, c = seq2NetParams(seq, tables)
    b = computeBindingMode(ligands.list[1], tables)
    return abs(h'b)
end


"""
calculated the change in the ground state between the solvent
bound and ligand bound conditions.
"""
function computeConfChange(net::Network, ligands::Ligands, tables::Tables)
    N  = length(tables.R)
    N_ligs = length(ligands.list)
    GS = zeros(N, N_ligs)
    for i in eachindex(ligands.list)
        per = ligands.list[i]
        computePerturbedNetParams!(net, per, tables)
        GS[:,i] = (net.G + net.δG) * (net.h + net.δh)
    end
    return ΔR = GS[:,2:end] .- GS[:,1];
end

"""
Overlaps of the normalized ligand-induced conformational changes (`computeConfChange`)
with the eigenvectors of `G = pinv(H)`, the network's normal modes. One column per
ligand.
"""
function computeOverlaps(net::Network, ligands::Ligands, tables::Tables)
    ΔR = computeConfChange(net::Network, ligands::Ligands, tables::Tables)
    ΔR .= ΔR ./ reshape(norm.(eachcol(ΔR)), 1, :)
    eig = eigen(net.G)
    return eig.vectors'ΔR
end

"""
Same as `computeOverlaps(net, ligands, tables)`, for the network of sequence `seq`.
"""
function computeOverlaps(seq::Vector, ligands::Ligands, tables::Tables)
    return computeOverlaps(buildNetwork(seq, tables), ligands, tables)
end

"""
Largest absolute overlap of a ligand-induced conformational change with a normal mode
of the network of `seq`.
"""
computeMaxOverlap(seq::Vector, ligands::Ligands, tables::Tables) = maximum(abs, computeOverlaps(seq, ligands, tables))


"""
    generalist(ΔF)

Generalist score `min(-ΔF[1], -ΔF[2])` from the binding energies `ΔF` to the two
ligands (as stored in `binding_energies.jld2`). Applied per sequence with
`mapslices(generalist, ΔF, dims=1)`.
"""
generalist(ΔF::AbstractVector) = min(-ΔF[1], -ΔF[2])
