# Commented-out code moved verbatim from src/physics.jl.
# Not included in the module; kept for reference.

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#function computeFreeEnergy(H::Matrix, h::Vector, c::Number, β::Number; N_sdims::Integer=2)
#    # compute the free energy of energy landscape with the form
#    # U(r) = 1/2 r'Hr + h'r + c
#    N_dims = length(h)
#    GSE = - 0.5 * dot(h, pinv(H), h) + c # ground state energy
#    lpd = logpdet(H; N_sdims) # this should be log(det(Hessian))
#    FLE = -0.5*( N_dims*log(2*pi/β) - lpd ) / β # Flucutation energy
#    F = GSE + FLE
#    return F
#end

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#function getNetParams_old(seq::Vector, tables::Tables)
#    # the old way to compute the network parameters
#    H, h, c = seq2NetParams(seq, tables)
#    return pinv(H), logpdet(H), h, c
#end

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

#computeHessian!(Q::Network) = computeHessian!(Q.H,Q.r,Q.SM)
#
