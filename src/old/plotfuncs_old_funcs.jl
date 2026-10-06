# Commented-out code moved verbatim from src/plotfuncs.jl.
# Not included in the module; kept for reference.

#plotNetwork!(ax, Q::Network) = plotNetwork!(ax, Q.r, Q.SM)

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#function plotNetwork(r::Vector, A::Matrix)
#    fig, ax = subplots()
#    plotNetwork!(ax, r, A)
#    return fig
#end

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#plotNetwork(tables::Tables) = plotNetwork(tables.R, tables.A)

#function plotNetworkColoredNodes!(ax, fig, r::Vector,
#                                  A::NumMatrix,
#                                  nodeColorVec::NumVector;
#                                  bondColor::Vector=[0.7,0.7,0.7],
#                                  bondWidth=3,
#                                  nodeSize=30,
#                                  cmap="Purples",
#                                  vmin=minimum(nodeColorVec),
#                                  vmax=maximum(nodeColorVec))
#
#    # plot a network, with the nodes painted according to nodeColorVec.
#    xy = r2xy(r)
#    x = xy[:,1]
#    y = xy[:,2]
#    N = size(A,1)
#    for i in 2:N, j in 1:i-1
#        if A[i,j] != 0
#            ax.plot([x[i],x[j]], [y[i],y[j]], lw=bondWidth, c=bondColor, zorder=0)
#        end
#    end
#    ax.scatter(x, y, nodeSize; c=nodeColorVec, zorder=2, cmap, vmin, vmax)
#    ax.axis("equal")
#    clb = addCustomColorBar!(ax, fig; cmap)
#    return clb
#end

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#    plotBondStiffness!(ax, fig, Q.r, Q.SM, Q.SM; bondWidth, nodeSizem, cmap)
#    return nothing

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#"""
#Plot population average negative binding energies over time.
#ΔF1 and ΔF2 should be binding energies, ie negative for tight binding
#"""
#function plotBindingEvoTraces(ΔF1::Vector,ΔF2::Vector; 
#                              τ=0, μ=0, P=0, selParams=(0,0))
#    fig, ax = subplots( figsize=(10,2))
#    ax.plot(-ΔF1, lw=1, "r", label="Lig 1")
#    ax.plot(-ΔF2, lw=1, "b", label="Lig 2")
#    ax.set_xlabel("generations", fontsize=14)
#    ax.set_ylabel(L"\langle -\Delta F \rangle_{\rm pop}", fontsize=14)
#    ax.set_xlim([0, length(ΔF1)])
#    ax.legend()
#    ax.axhline(y=0, color="grey", linestyle="--", linewidth=0.8, zorder=0)
#    ax.set_title(L"\tau=" * @sprintf("%.3g", τ) * ", " *
#             L"\mu=" * @sprintf("%.3g", μ) * ", " *
#             L"P=" * @sprintf("%.3g", P) * ", " *
#             L"a=" * @sprintf("%.3g", selParams[1]) * ", " *
#             L"x_0=" * @sprintf("%.3g", selParams[2]))
#    return fig
#end

# UNUSED: not called from exps/*/scripts or notebooks (commented out, not deleted)
#"""
#Plot the number of dead, specialist and generalist sequences over time for
#an evolutionary trajectory.

#fits[:,:,i] should show the population negative binding energies at generation i.
#"""
#function plotGeneralistSpecialistEvoTraces(fits::AbstractArray{T,3}) where T


#    function binByQuadrants(fits::Array{T,3}) where T
#        # Classify each sequence as dead specialist, or generalist
#        #counts[1,1] is the number of "dead" sequences
#        #counts[2,1] is the number of ligand 1 specialist sequences
#        #counts[1,2] is the number of ligand 1 specialist sequences
#        #counts[2,2] is the number of generalist sequences

#        function kernel(fits::Matrix)
#            counts = zeros(2,2)
#            for i in axes(fits, 1)
#                x, y = fits[i,:] .> 0
#                counts[x+1, y+1] += 1
#            end
#            return counts
#        end

#        @assert size(fits, 2) == 2
#        fullCounts = zeros(2,2,size(fits,3))
#        for i in 1:size(fits,3)
#           fullCounts[:,:,i] = kernel(fits[:,:,i])
#        end
#        return fullCounts
#    end
#    counts = binByQuadrants(fits);
#    fig, ax = subplots( figsize=(10,2))
#    ax.plot(counts[1,1,:], lw=0.5, color="grey")
#    ax.plot(counts[2,1,:], lw=0.5, color="red")
#    ax.plot(counts[1,2,:], lw=0.5, color="blue")
#    ax.plot(counts[2,2,:], lw=1, color="green")
#    #ax.set_xlabel("generations", fontsize=14)
#    #ax.set_ylabel(L"\langle -\Delta F \rangle_{\rm pop}", fontsize=14)
#    ax.set_xlim([0, size(fits, 3)])
#    #ax.legend()
#    ax.axhline(y=0, color="grey", linestyle="--", linewidth=0.8, zorder=0)
#    return fig, counts
#end
