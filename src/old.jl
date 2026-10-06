#function computeHessian(r::Vector{Float64},
#                        SM::Matrix{Float64})
#    # compute Hessian matrix of Hamiltonian
#    # returns hessian
#    N = size(SM,1)
#    H = zeros(2N, 2N)
#    for i in 1:N
#        sxx = 0.0
#        syy = 0.0 
#        sxy = 0.0
#        i2 = 2i
#        xi = r[i2-1]
#        yi = r[i2]
#        for j in 1:N
#            k = SM[i,j]
#            if k > 0 && i != j
#                j2 = 2j
#                dx = xi - r[j2-1]
#                dy = yi - r[j2]
#                c = -k / (dx^2 + dy^2)
#                Exx = c * dx^2 
#                Eyy = c * dy^2 
#                Exy = c * dx*dy 
#                sxx -= Exx
#                syy -= Eyy
#                sxy -= Exy
#                H[i2-1,j2-1] = Exx
#                H[i2,j2] = Eyy
#                H[i2,j2-1] = Exy
#                H[i2-1,j2] = Exy
#            end
#        end
#        H[i2-1,i2-1] = sxx
#        H[i2,i2] = syy
#        H[i2,i2-1] = sxy
#        H[i2-1,i2] = sxy
#    end
#    return H
#end


#function computeHessian!(H::Matrix,
#                         r::Vector,
#                         SM::Matrix)
#    # compute the hessian matrix of the network
#    N = size(SM,1)
#    for i in 1:N
#        sxx = 0.0
#        syy = 0.0 
#        sxy = 0.0
#        i2 = 2i
#        xi = r[i2-1]
#        yi = r[i2]
#        for j in 1:N
#            
#            j2 = 2j
#            # fill H with zeros
#            H[i2-1,j2-1] = 0.0
#            H[i2,j2-1] = 0.0
#            H[i2-1,j2] = 0.0
#            H[i2,j2] = 0.0
#
#            k = SM[i,j]
#            if k > 0 && i != j
#                dx = xi - r[j2-1]
#                dy = yi - r[j2]
#                c = -k / (dx^2 + dy^2)
#                Exx = c * dx^2 
#                Eyy = c * dy^2 
#                Exy = c * dx*dy 
#                sxx -= Exx
#                syy -= Eyy
#                sxy -= Exy
#                H[i2-1,j2-1] = Exx
#                H[i2,j2-1] = Exy
#                H[i2-1,j2] = Exy
#                H[i2,j2] = Eyy
#            end
#        end
#        H[i2-1,i2-1] = sxx
#        H[i2,i2-1] = sxy
#        H[i2-1,i2] = sxy
#        H[i2,i2] = syy
#    end
#    return nothing
#end

#function sampleSequences!(seqBucket, nets, sampleLastHalf, t, gen2StartSampling, numSeqs2KeepPerGen)
#    # sample sequence from the last half of evolution.
#    if sampleLastHalf && t >= gen2StartSampling
#        inds = randperm(P)[1:numSeqs2KeepPerGen]
#        # [copy(nets[i].seq) for i in inds]
#        selectedSeqs = [copy(nets[i].seq) for i in inds]
#        append!(seqBucket, selectedSeqs)
#    end
#end


#function samplingDetails(N, P, numSeqs2Keep)
#    # some book keeping for sampling sequence from the population
#    # in evolvePop.
#    numGen2Sample = Int(floor(N / 2))
#    maxSeqsSampleSize = P * numGen2Sample
#    if numSeqs2Keep > maxSeqsSampleSize
#        println("You are trying to keep too many sequences\n"*
#                "Setting numSeqs2Keep = maxSeqsSampleSize")
#        numSeqs2Keep = maxSeqsSampleSize
#    end
#    numSeqs2KeepPerGen = Int(ceil(numSeqs2Keep / numGen2Sample))
#    gen2StartSampling = Int(ceil(N/2))
#    return numGen2Sample,maxSeqsSampleSize,numSeqs2KeepPerGen, gen2StartSampling,numSeqs2Keep
#end


#function computeResponse_LM( H::FloatMatrix, site::IntVector, strain::FloatVector)
#    # compute the response via Lagrange Multiplier.
#    # dr0 is the imposed displacment.
#    # H is the hessian matrix.
#
#    numDims = size(H,1)
#    @assert length(site) == length(strain)
#    n = length(site)
#
#    # build constraints vector
#    b = zeros(numDims+n) # 
#    b[numDims+1:end] = strain
#
#    B = zeros(numDims+n,numDims+n)
#    B[1:numDims,1:numDims] .= H
#    B[site+1:end, site] = I(n)
#    B[site, site+1:end] = I(n);
#
#    # solve Linear Response problem
#    y = B \ b
#
#    # construct full diplacement
#    dr = y[1:numDims]
#
#    # compute force
#    f = H * dr
#
#    # calculate energy
#    E = f ⋅ dr
#
#    return E, dr, f
#end
#
#
##function computePerturbedNetParams!(net::Network, per::Perturbation, tables::Tables)
#    
#    @assert length(per.bonds) == length(per.Δk) == length(per.Δd)
#    # return the perturbed network parameters 
#    seq, G, h, lpdH = net.seq, net.G, net.h, net.lpdH
#    A, R, K, D, E = tables.A, tables.R, tables.K, tables.D, tables.E
#    nzinds = findall(x-> x!=0, per.Δk) # have to remove zero bonds 
#    bonds, Δk, Δd = per.bonds[nzinds], per.Δk[nzinds], per.Δd[nzinds] 
#    N_m = length(bonds) 
#    # Build local spring matricies
#    T_m = Matrix{Float64}(undef, N_m, length(R)) ###################
#    k = Vector{Float64}(undef, N_m) 
#    d = Vector{Float64}(undef, N_m) 
#    v = Vector{Float64}(undef, N_m) 
#    J = zeros(N_m, N_m) 
#    for i in 1:N_m
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
#    for i in 1:N_m
#        net.δc += k[i]*d[i]*Δd[i] + 0.5*k[i]*Δd[i]^2 + 0.5*Δk[i]*d[i]^2 + Δk[i]*d[i]*Δd[i] + 0.5*Δk[i]*Δd[i]^2
#    end
#
#    # compute update to h
#    v .= Δk .* d .+ k .* Δd .+ Δk .* Δd
#    mul!(net.δh, T_m', v)
#
#    # compute and update change in  G
#    for i in 1:N_m
#        J[i,i]= 1/Δk[i]
#    end 
#    
#    Y = T_m*Symmetric(G) ############
#    #Y = T_m*G ############
#
#    #J0 =  lu(J + Y*T_m') 
#    J0 = J + Y*T_m'
#    mul!( net.δG, -Y', J0\Y)
#
#    # compute change in log determinate
#    det_Δk = prod(Δk)
#    net.δlpdH = log( det_Δk * det(J0) )
#    return nothing
#end

