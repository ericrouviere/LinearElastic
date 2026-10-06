

"""
    computeFreeEnergy(G, lpdH, h, c, β)

Free energy of the quadratic energy landscape `U(r) = ½ r'Hr - h'r + c` at inverse
temperature `β`, given `G = pinv(H)` and `lpdH = log pdet(H)`:

    F = -½ h'Gh + c - (n log(2π/β) - lpdH) / 2β,

the ground-state energy plus the free energy of thermal fluctuations, `n = length(h)`.
"""
function computeFreeEnergy(G::Matrix, lpdH::Number, h::Vector, c::Number, β::Number; N_sdims::Integer=2)
    N_dims = length(h)
    GSE = - 0.5 * dot(h, G, h) + c # ground state energy
    FLE = -0.5*( N_dims*log(2*pi/β) - lpdH ) / β # Flucutation energy
    return GSE + FLE
end

"""
Free energy of the unperturbed network `net`.
"""
function computeFreeEnergy(net::Network, β::Number)
    return computeFreeEnergy(net.G, net.lpdH, net.h, net.c, β)
end

"""
    computeFreeEnergy(net, per, tables, β)

Free energy of `net` perturbed by `per` (e.g. a bound ligand), from a low-rank update
of the cached parameters (`computePerturbedNetParams!`). Overwrites the scratch fields
`net.δG` and `net.δh` with `G + δG` and `h + δh`.
"""
function computeFreeEnergy(net::Network, per::Perturbation, tables::Tables, β)
    computePerturbedNetParams!(net, per, tables)
    # here I am using the allocated memory of δG and δh to store G+δG and h+δh.
    net.δG .+= net.G
    net.δh .+= net.h
    return computeFreeEnergy(net.δG, net.lpdH+net.δlpdH, net.δh, net.c+net.δc, β)
end

"""
Free energies of `net` in every state of `ligands`, solvent first.
"""
function computeFreeEnergies(net::Network, ligands::Ligands, tables::Tables, β::Number)
    # THIS CAN BE OPTIMIZEDF FOR <10% IMPROVEMENT. 
    return computeFreeEnergy.([net], ligands.list, [tables], β) 
end


# -----------------------------------------------------------------------------

"""
Log pseudo-determinant of the square matrix `A`: the sum of the logs of its
eigenvalues, leaving out the `N_sdims(N_sdims + 1)/2` zero modes of rigid-body motion
(3 in 2D). Slow; warns if the smallest eigenvalue is not zero.
"""
function logpdet(A::AbstractMatrix{T}; N_sdims::Integer=2) where T
    @assert size(A,1) == size(A,2)
    N_zeromodes = Int(N_sdims* (N_sdims+1) / 2)
    vals = eigvals(A)
    if abs(vals[1]) > 1e-13
        println("Something Wrong with spectrum of Matrix")
    end
    
    lpd = 0.0
    @inbounds for i in N_zeromodes+1:length(vals)
        lpd += log(real(vals[i]))
    end
    return lpd
    #vals = abs.(real.(vals))
    #return sum(log.(vals[N_zeromodes+1:end]))
end


"""
Log pseudo-determinant from a precomputed SVD `J` of a symmetric positive
semidefinite matrix: the sum of the logs of the singular values, leaving out the
smallest `N_sdims(N_sdims + 1)/2` (the zero modes).
"""
function logpdet(J::SVD{T, T, Matrix{T}}; N_sdims::Integer=2) where T
    N_zeromodes = Int(N_sdims* (N_sdims+1) / 2)
    if abs(J.S[end]) > 1e-13
        println("Something Wrong with spectrum of Matrix")
    end
    
    lpd = 0.0
    @inbounds for i in 1:length(J.S)-N_zeromodes
        lpd += log(J.S[i])
    end
    return lpd
end

"""
Moore-Penrose pseudoinverse from a precomputed SVD `J`, adapted from
`LinearAlgebra.pinv`: singular values below `max(rtol * maximum(J.S), atol)` are
treated as zero.
"""
function pinv(J::SVD{T, T, Matrix{T}, Vector{T}};
              atol::Real = 0.0,
              rtol::Real=eps(real(float(oneunit(T))))*min(size(J.U,1), size(J.Vt,2))*iszero(atol)
             ) where T
    # code from source in LinearAlgrebra
    tol2        = max(rtol*maximum(J.S), atol)
    Stype       = eltype(J.S)
    Sinv        = copy(J.S)
    #Sinv        = fill!(similar(A, Stype, length(J.S)), 0)
    index       = J.S .> tol2
    Sinv[index] .= pinv.(view(J.S, index))
    return J.Vt' * (Diagonal(Sinv) * J.U')
end

# -----------------------------------------------------------------------------

"""
    makeBeads2BondsMatrix(A, r)

Beads-to-bonds matrix `T` of the network with adjacency `A` and positions `r`, so that
the Hessian is `H = T'ΓT`. Row `l` holds the unit vector of bond `l` (positive on its
first bead, negative on its second), so `T * u` is the elongation of every bond under
the displacement `u`. Bonds are ordered as `for j in 1:N-1, i in j+1:N`.
"""
function makeBeads2BondsMatrix(A::Matrix, r::Vector) 
    d = 2 # spatial dimensions
    N_bonds = Int(sum(A)/2)
    N_dims = length(r)
    N_beads = size(A,1)
    T = zeros(N_bonds, N_dims)
    l = 0
    for j in 1:N_beads-1, i in j+1:N_beads
        if A[i,j] == 1
            l += 1
            x1 = r[d*i-1]
            y1 = r[d*i] 
            x2 = r[d*j-1]
            y2 = r[d*j]
            Δr = sqrt( (x1-x2)^2 + (y1-y2)^2 )
            T[l, d*i-1] = (x1-x2)/Δr
            T[l, d*i] = (y1-y2)/Δr
            T[l, d*j-1] = (x2-x1)/Δr
            T[l, d*j] = (y2-y1)/Δr
        end
    end
    return T
end

"""
Unit vector `(nx, ny)` along the bond between `bead_1` and `bead_2`, pointing from
`bead_2` to `bead_1`.
"""
function computeBondVector(R::Vector, bead_1::Integer, bead_2::Integer; N_sdims=2)
    @assert N_sdims==2 # MAKE MORE GENERAL
    xi = R[ N_sdims*(bead_1-1) + 1]
    yi = R[ N_sdims*(bead_1-1) + 2]
    xj = R[ N_sdims*(bead_2-1) + 1]
    yj = R[ N_sdims*(bead_2-1) + 2]
    Δx = xi-xj
    Δy = yi-yj
    inv_dist = 1/sqrt(Δx^2 + Δy^2)
    # n = inv_dist * [Δx, Δy] # the the unit bond vector
    return Δx*inv_dist, Δy*inv_dist
end

"""
Displacement mode `b` (length `2N`) that stretches the bond between `bead_1` and
`bead_2`: `+n` on `bead_1` and `-n` on `bead_2`, with `n` from `computeBondVector`.
It is the bond's row of `makeBeads2BondsMatrix`.
"""
function computeBondElongationMode(R::Vector, bead_1::Integer, bead_2::Integer; N_sdims=2)
    @assert N_sdims==2 # MAKE MORE GENERAL
    nx, ny = computeBondVector(R, bead_1, bead_2; N_sdims)
    b = zeros(length(R))
    b[N_sdims*(bead_1-1)+1] = nx
    b[N_sdims*(bead_1-1)+2] = ny
    b[N_sdims*(bead_2-1)+1] = -nx
    b[N_sdims*(bead_2-1)+2] = -ny
    return b
end
 
# -----------------------------------------------------------------------------

"""
    seq2NetParams(seq, tables)

Parameters `H, h, c` of the quadratic energy `U(r) = ½ r'Hr - h'r + c` of sequence
`seq`, the sum over bonds of `½ k_l (e_l - d_l)²` with `e_l` the elongation of bond
`l`. Bond `l` between beads `i` and `j` has stiffness `k_l = K[seq[i], seq[j]]` and
mismatch `d_l = D[seq[i], seq[j]]`; with `Γ = diag(k)`, `M = diag(d)` and the
beads-to-bonds matrix `T`:

    H = T'ΓT,   h = T'ΓM·1,   c = ½·1'ΓM²1.
"""
function seq2NetParams(seq, tables::Tables)
    R, A, K, D, E = tables.R, tables.A, tables.K, tables.D, tables.E
    N_beads = size(A,1)
    N_bonds = Int(sum(A)/2)
    T = makeBeads2BondsMatrix(A, R)
    Γ = Diagonal(zeros(N_bonds)) # diagonal Stiffness matrix
    M = Diagonal(zeros(N_bonds))
   # C = zeros(N_bonds, N_bonds)
    l = 0
    @inbounds for j in 1:N_beads-1, i in j+1:N_beads
        if A[i,j] == 1
            l += 1
            Γ[l,l] = K[seq[i], seq[j]]
            M[l,l] = D[seq[i], seq[j]]
            #C[l,l] = E[seq[i], seq[j]]
        end
    end
    ones_vec = ones(N_bonds)
    H = T'Γ*T
    h = T'*(Γ*M)*ones_vec
    c = 0.5*ones_vec'*(Γ*M*M)*ones_vec
    return H, h, c
end


"""
Cached parameters `G = pinv(H)`, `lpdH = log pdet(H)`, `h` and `c` of sequence `seq`,
from a single SVD of `H`.
"""
function getNetParams(seq::Vector, tables::Tables)
    H, h, c = seq2NetParams(seq, tables)
    J = svd(Symmetric(H))
    return pinv(J), logpdet(J), h, c
end

"""
Recompute all cached parameters of `net` from scratch for its current sequence,
discarding the round-off accumulated by incremental mutations.
"""
function updateNetParams!(net::Network, tables::Tables)
    G, lpdH, h, c = getNetParams(net.seq, tables)
    net.G .= G
    net.lpdH = lpdH
    net.h .= h
    net.c = c
    return nothing
end

# -----------------------------------------------------------------------------


"""
    computePerturbedNetParams!(net, per, tables)

Change in the energy parameters of `net` caused by the perturbation `per`, stored in
`net.δG`, `net.δlpdH`, `net.δh` and `net.δc`; the cached parameters are unchanged.
With the unit vectors of the perturbed bonds as the rows of `T_m`, the Hessian changes
by `T_m' diag(Δk) T_m`, so by the Woodbury identity and the matrix determinant lemma

    δG = -G T_m' J⁻¹ T_m G,   δlpdH = log(det(diag(Δk)) det(J)),   J = diag(Δk)⁻¹ + T_m G T_m',

over the bonds with `Δk ≠ 0`. `h` and `c` change by the terms of the perturbed bonds.
This is a rank-`m` update for `m` perturbed bonds instead of a new SVD.
"""
function computePerturbedNetParams!(net::Network, per::Perturbation, tables::Tables)


    function updateCurvature!(net::Network, G::Matrix, T::Matrix, Δk::Vector) 
        Y = T*Symmetric(G) # compute and update change in G
        J = Y*T'
        @inbounds for i in eachindex(Δk)
            J[i,i] += 1/Δk[i]
        end
        net.δlpdH = log( prod(Δk) * det(J) ) # compute change in log determinate
        J .*= -1.0
        mul!( net.δG, Y', J\Y)
        return nothing
    end

    @assert length(per.bonds) == length(per.Δk) == length(per.Δd)
    # return the perturbed network parameters 
    seq, G, h, lpdH = net.seq, net.G, net.h, net.lpdH
    A, R, K, D, E = tables.A, tables.R, tables.K, tables.D, tables.E
    bonds, Δk, Δd = per.bonds, per.Δk, per.Δd
    N_m = length(Δk)

    # Build local spring matricies
    T_m = Matrix{Float64}(undef, N_m, length(R)) ###################
    k = Vector{Float64}(undef, N_m) 
    d = Vector{Float64}(undef, N_m) 
    v = Vector{Float64}(undef, N_m) 
    @inbounds for i in 1:N_m
        bead_1, bead_2 = bonds[i]
        beadType_1 = seq[bead_1]
        beadType_2 = seq[bead_2]
        a = A[bead_1, bead_2]
        T_m[i,:] = computeBondElongationMode(R, bead_1,  bead_2);
        k[i] =  a * K[beadType_1, beadType_2]
        d[i] = a * D[beadType_1, beadType_2]
    end
    
    # compute update to c
    net.δc = 0.0
    @inbounds for i in 1:N_m
        net.δc += k[i]*d[i]*Δd[i] + 0.5*k[i]*Δd[i]^2 + 0.5*Δk[i]*d[i]^2 + Δk[i]*d[i]*Δd[i] + 0.5*Δk[i]*Δd[i]^2
    end

    # compute update to h
    v .= Δk .* d .+ k .* Δd .+ Δk .* Δd
    mul!(net.δh, T_m', v)

    # update δG and δlpdh
    nzinds = findall(x-> x!=0, Δk) # find non zero stiffnesses
    N_nz = length(nzinds)
    if N_nz < N_m
        T_nz = T_m[nzinds, :]
        Δk_nz = Δk[nzinds]
        updateCurvature!(net, G, T_nz, Δk_nz)
    else
        updateCurvature!(net, G, T_m, Δk)
    end
    return nothing
end


"""
Pseudo-determinant of the diagonal matrix `A`: the product of its diagonal entries of
magnitude above `thresh`. Used by `computePerturbedNetParams_old`.
"""
function pdetdiagm(A::Matrix{T}, thresh=1e-12) where T<:Real
    pdet = one(eltype(A))
    for i in 1:size(A,1)
        if abs(A[i,i]) > thresh
            pdet *= A[i,i]
        end
    end
    return pdet
end

"""
Bead positions in the ground state without and with the perturbation `per`:
`R0 = R + Gh` and `R1 = R + (G + δG)(h + δh)`.
"""
function computePerturbedGroundState(net::Network, per::Perturbation, tables::Tables)
    computePerturbedNetParams!(net, per, tables)
    R0 = net.G * net.h + tables.R
    R1 = (net.G + net.δG) * (net.h + net.δh) + tables.R
    return R0, R1
end

# ------------------------------------------------------------------

"""
    computeBindingRandSeqs(N_seq, ligands, tables, β)

Binding energies of `N_seq` uniformly random sequences to the two ligands of
`ligands` (solvent plus two bound states). Returns `ΔF, s, m`: the `N_seq × 2`
binding energies, the spread of the cloud (standard deviation of the distances of
the points to its centroid), and the mean of the centroid's two coordinates.
"""
function computeBindingRandSeqs(N_seq::Int, ligands::Ligands, tables::Tables, β)
    function computeMeanStdOfPoints(points::Matrix{Float64})
        # Ensure the input is a 1000x2 matrix
        @assert size(points, 2) == 2 "Input must be a Nx2 matrix where each row is a point in 2D space."
        centroid = mean(points, dims=1)
        distances = sqrt.(sum((points .- centroid).^2, dims=2))
        return mean(centroid), std(distances)
    end

    @assert length(ligands.list) == 3
    ΔF = zeros(N_seq, length(ligands.list)-1)
    net = buildRandNetwork(tables)
    for i in 1:N_seq
        net.seq .= randSeq(tables.R, tables.K)
        updateNetParams!(net, tables)
        energies = computeFreeEnergies(net, ligands, tables, β)
        ΔF[i,1] = energies[2] - energies[1]
        ΔF[i,2] = energies[3] - energies[1]
    end
    m, s = computeMeanStdOfPoints(ΔF)
    return ΔF, s, m
end

"""
Binding free energies `ΔF_k = F_k - F_1` of `net` to each bound state `k ≥ 2` of
`ligands`, relative to the solvent state `ligands.list[1]`. Negative means binding.
"""
function computeBindingEnergies(net::Network, ligands::Ligands, tables::Tables, β::Number)
    F = computeFreeEnergy.([net], ligands.list, [tables], β) 
    return F[2:end] .- F[1]
end


############################################################
### Old Code ###############################################
############################################################


"""
Fill `H` with the Hessian of a network with positions `r` and bond stiffness matrix
`SM` (`SM[i, j]` is the stiffness of the bond between beads `i` and `j`, 0 if none),
built bead pair by bead pair. An independent check of `seq2NetParams`, used by the
unit-test notebook.
"""
function computeHessian!(H::AbstractMatrix{Float64},
                         r::AbstractVector{Float64},
                         SM::AbstractMatrix{Float64})
    N = size(SM,1) # Number of beads in network
    @inbounds for j in 1:N
        sxx = 0.0
        syy = 0.0 
        sxy = 0.0
        j2 = 2j
        xj = r[j2-1]
        yj = r[j2]
        @fastmath for i in 1:N
            i2 = 2i
            # fill H with zeros
            H[i2-1,j2-1] = 0.0
            H[i2,j2] = 0.0
            H[i2,j2-1] = 0.0
            H[i2-1,j2] = 0.0
            k = SM[i,j]
            if k > 0 && i != j
                dx = r[i2-1] - xj
                dy = r[i2] - yj
                c = -k / (dx^2 + dy^2)
                Exx = c * dx^2 
                Eyy = c * dy^2 
                Exy = c * dx*dy 
                sxx -= Exx
                syy -= Eyy
                sxy -= Exy
                H[i2-1,j2-1] = Exx
                H[i2,j2] = Eyy
                H[i2,j2-1] = Exy
                H[i2-1,j2] = Exy
            end
        end
        H[j2-1,j2-1] = sxx
        H[j2,j2] = syy
        H[j2,j2-1] = sxy
        H[j2-1,j2] = sxy
    end
    return nothing
end

"""
Allocating version of `computeHessian!`.
"""
function computeHessian(r::Vector, SM::Matrix)
    N = size(SM,1)
    H = Matrix{Float64}(undef, 2N, 2N)
    computeHessian!(H,r,SM)
    return H
end

"""
Slow reference implementation of `computePerturbedNetParams!` with explicit
pseudoinverses. Returns `δG, δlpdH, δh, δc` instead of storing them. Used by the
unit-test notebook.
"""
function computePerturbedNetParams_old(net::Network, per::Perturbation, tables::Tables)
    seq, G, h, lpdH = net.seq, net.G, net.h, net.lpdH
    A, R, K, D, E = tables.A, tables.R, tables.K, tables.D, tables.E
    nzinds = findall(x-> x!=0, per.Δk) # have to remove zero bonds  #################### 
    bonds, Δk, Δd = per.bonds[nzinds], per.Δk[nzinds], per.Δd[nzinds] #################

    N_m = length(bonds) # 
    #N_beads = length(seq)
    
    # Build local spring matricies
    T_m = zeros(N_m, length(R)) ###################
    Γ_m = zeros(N_m, N_m) #################
    δΓ_m = zeros(N_m, N_m) ################
    M_m = zeros(N_m, N_m) #################
    δM_m = zeros(N_m, N_m) ################
    ones_vec = ones(N_m) ################
    for i in 1:N_m
        bead_1, bead_2 = bonds[i]
        beadType_1 = seq[bead_1]
        beadType_2 = seq[bead_2]
        T_m[i,:] = computeBondElongationMode(R, bead_1,  bead_2);
        Γ_m[i,i] = A[bead_1, bead_2] * K[beadType_1, beadType_2]
        δΓ_m[i,i] = Δk[i]
        M_m[i,i] = A[bead_1, bead_2] * D[beadType_1, beadType_2]
        δM_m[i,i] = Δd[i]
    end
    # compute update to c
    δc = 0.5*ones_vec'*( 2Γ_m*M_m*δM_m + Γ_m*δM_m^2 + δΓ_m*M_m^2 + 2δΓ_m*M_m*δM_m + δΓ_m*δM_m^2 )*ones_vec  ##########
    # compute update to h
    δh = T_m'*( δΓ_m*M_m + Γ_m*δM_m + δΓ_m*δM_m)*ones_vec ##############
    # compute change in  G
    δG = -G*T_m'* pinv(pinv(δΓ_m) + T_m*G*T_m')*T_m*G ############3
    # compute and update the determinants
    δlpdH = log( pdetdiagm(δΓ_m) * det( pinv(δΓ_m) + T_m*G*T_m') )
    return δG, δlpdH, δh, δc
end

