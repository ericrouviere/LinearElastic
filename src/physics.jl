
function computeFreeEnergy(H::Matrix, h::Vector, c::Number, β::Number; N_sdims::Integer=2)
    # compute the free energy of energy landscape with the form
    # U(r) = 1/2 r'Hr + h'r + c
    N_dims = length(h)
    GSE = - 0.5 * dot(h, pinv(H), h) + c # ground state energy
    lpd = logpdet(H; N_sdims) # this should be log(det(Hessian))
    FLE = -0.5*( N_dims*log(2*pi/β) - lpd ) / β # Flucutation energy
    F = GSE + FLE
    return F
end

function computeFreeEnergy(G::Matrix, lpdH::Number, h::Vector, c::Number, β::Number; N_sdims::Integer=2)
    # compute the free energy of energy landscape with the form
    # U(r) = 1/2 r' G^+ r + h'r + c
    # G is the psuedo inverse of the hessian H
    # lpdH is the log(pdet(H)) 
    N_dims = length(h)
    GSE = - 0.5 * dot(h, G, h) + c # ground state energy
    FLE = -0.5*( N_dims*log(2*pi/β) - lpdH ) / β # Flucutation energy
    return GSE + FLE
end

function computeFreeEnergy(net::Network, β::Number)
    return computeFreeEnergy(net.G, net.lpdH, net.h, net.c, β)
end

function computeFreeEnergy(net::Network, per::Perturbation, tables::Tables, β)
    # return the free energy of the perturbed system.
    computePerturbedNetParams!(net, per, tables)
    # here I am using the allocated memory of δG and δh to store G+δG and h+δh.
    net.δG .+= net.G
    net.δh .+= net.h
    return computeFreeEnergy(net.δG, net.lpdH+net.δlpdH, net.δh, net.c+net.δc, β)
end

function computeFreeEnergies(net::Network, ligands::Ligands, tables::Tables, β::Number)
    # THIS CAN BE OPTIMIZEDF FOR <10% IMPROVEMENT. 
    # return the free energy of the system perturbed by each ligand perturbation
    return computeFreeEnergy.([net], ligands.list, [tables], β) 
end


# -----------------------------------------------------------------------------

function logpdet(A::AbstractMatrix{T}; N_sdims::Integer=2) where T
    # return the log of the pseudo determinant of A
    # This is a very slow function!!
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


function logpdet(J::SVD{T, T, Matrix{T}}; N_sdims::Integer=2) where T
    # return the log of the pseudo determinant of A
    # This is a very slow function!!
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

function pinv(J::SVD{T, T, Matrix{T}, Vector{T}};
              atol::Real = 0.0,
              rtol::Real=eps(real(float(oneunit(T))))*min(size(J.U,1), size(J.Vt,2))*iszero(atol)
             ) where T
    # compute the pinv using the precomputed SVD object J.
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

function makeBeads2BondsMatrix(A::Matrix, r::Vector) 
    # return the beads to bonds matrix T of the
    # Hessian decomposion H = T' * Γ * T.
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

function computeBondVector(R::Vector, bead_1::Integer, bead_2::Integer; N_sdims=2)
    # return the unit vector pointing from bead_1 to bead_2
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

function computeBondElongationMode(R::Vector, bead_1::Integer, bead_2::Integer; N_sdims=2)
    # return the mode of the network that corresponds to the elongation 
    # of bond specified between beads i and j.
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

function seq2NetParams(seq, tables::Tables)
    #  Convert a sequence to the pargameters of the quadratice energy
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

function getNetParams_old(seq::Vector, tables::Tables)
    # the old way to compute the network parameters
    H, h, c = seq2NetParams(seq, tables)
    return pinv(H), logpdet(H), h, c
end

function getNetParams(seq::Vector, tables::Tables)
    # compute the network parameters using 1 svd.
    H, h, c = seq2NetParams(seq, tables)
    J = svd(Symmetric(H))
    return pinv(J), logpdet(J), h, c
end

function updateNetParams!(net::Network, tables::Tables)
    # compute the values of the network parameters and update the struct
    G, lpdH, h, c = getNetParams(net.seq, tables)
    net.G .= G
    net.lpdH = lpdH
    net.h .= h
    net.c = c
    return nothing
end

# -----------------------------------------------------------------------------


#function computePerturbedNetParams!(net::Network, per::Perturbation, tables::Tables)
#    
#    @assert length(per.bonds) == length(per.Δk) == length(per.Δd)
#    # return the perturbed network parameters 
#    seq, G, h, lpdH = net.seq, net.G, net.h, net.lpdH
#    A, R, K, D, E = tables.A, tables.R, tables.K, tables.D, tables.E
#
#    
#
#    nzinds = findall(x-> x!=0, per.Δk) # have to remove zero bonds 
#    bonds, Δk, Δd = per.bonds[nzinds], per.Δk[nzinds], per.Δd[nzinds] 
#    N_m = length(bonds) 
#    # Build local spring matricies
#    T_m = Matrix{Float64}(undef, N_m, length(R)) ###################
#    Y = similar(T_m)
#    k = Vector{Float64}(undef, N_m) 
#    d = Vector{Float64}(undef, N_m) 
#    v = Vector{Float64}(undef, N_m) 
#    @show N_m
#    @inbounds for i in 1:N_m
#        bead_1, bead_2 = bonds[i]
#        beadType_1 = seq[bead_1]
#        beadType_2 = seq[bead_2]
#        a = A[bead_1, bead_2]
#        T_m[i,:] = computeBondElongationMode(R, bead_1,  bead_2);
#        k[i] =  a * K[beadType_1, beadType_2]
#        d[i] = a * D[beadType_1, beadType_2]
#    end
#    # compute update to c
#    net.δc = 0.0
#    @inbounds for i in 1:N_m
#        net.δc += k[i]*d[i]*Δd[i] + 0.5*k[i]*Δd[i]^2 + 0.5*Δk[i]*d[i]^2 + Δk[i]*d[i]*Δd[i] + 0.5*Δk[i]*Δd[i]^2
#    end
#
#    # compute update to h
#    v .= Δk .* d .+ k .* Δd .+ Δk .* Δd
#    @show v
#    mul!(net.δh, T_m', v)
#    @show net.δh
#    
#    # compute and update change in  G
#    mul!(Y,T_m, Symmetric(G))
#    J = Y*T_m'
#    @inbounds for i in 1:N_m
#        J[i,i] += 1/Δk[i]
#    end
#
#    # compute change in log determinate
#    net.δlpdH = log( prod(Δk) * det(J) )
#    J .*= -1.0
#    mul!( net.δG, Y', J\Y)
#    return nothing
#end


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


function pdetdiagm(A::Matrix{T}, thresh=1e-12) where T<:Real
    # return the pseudodeterminant of sqare diagonal matrix A.
    # this function used in slow version of computePerturbedNetParams
    pdet = one(eltype(A))
    for i in 1:size(A,1)
        if abs(A[i,i]) > thresh
            pdet *= A[i,i]
        end
    end
    return pdet
end

function computePerturbedGroundState(net::Network, per::Perturbation, tables::Tables)
    # return the unperturbed and perturbed ground states
    computePerturbedNetParams!(net, per, tables)
    R0 = net.G * net.h + tables.R
    R1 = (net.G + net.δG) * (net.h + net.δh) + tables.R
    return R0, R1
end

# ------------------------------------------------------------------

function computeBindingRandSeqs(N_seq::Int, ligands::Ligands, tables::Tables, β)
    # generate N_seq random sequences and compute their binding energies.
    # returns ΔF, s, m: the binding energies, the std of the distance to the
    # centroid of the ΔF cloud, and the mean of the centroid.
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

function computeBindingEnergies(net::Network, ligands::Ligands, tables::Tables, β::Number)
    # return the binding free energy of the system perturbed by each ligand perturbation
    # assuming first perturbation in ligands is the solvent (reference)
    F = computeFreeEnergy.([net], ligands.list, [tables], β) 
    return F[2:end] .- F[1]
end



############################################################
### Old Code ###############################################
############################################################


function computeHessian!(H::AbstractMatrix{Float64},
                         r::AbstractVector{Float64},
                         SM::AbstractMatrix{Float64})
    # compute the hessian matrix of the network
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
#computeHessian!(Q::Network) = computeHessian!(Q.H,Q.r,Q.SM)
#
function computeHessian(r::Vector, SM::Matrix)
    N = size(SM,1)
    H = Matrix{Float64}(undef, 2N, 2N)
    computeHessian!(H,r,SM)
    return H
end

function computePerturbedNetParams_old(net::Network, per::Perturbation, tables::Tables)
    # return the perturbed network parameters 
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

