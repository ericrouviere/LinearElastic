function plotNetwork!(ax, r::Vector, adjacencyMatrix::Matrix;
                      bondColor::Vector=[0.7,0.7,0.7],
                      bondWidth=3,
                      nodeSize=30)
    # Plots network, modifies figure axis object in place
    # uses PyPlot, LinearAlgebra
    #xy = r2xy(r)
    xy = reshape(r, 2, :)'
    x = xy[:,1]
    y = xy[:,2]
    N = size(adjacencyMatrix,1)
    A = LowerTriangular(adjacencyMatrix)
    for i =1:N
        whichj = findall(A[i,:].>0)
        for j in whichj
            if A[i,j] != 0
                ax.plot([x[i],x[j]], [y[i],y[j]], lw=bondWidth, c=bondColor, zorder=0)
            end
        end
    end
    ax.scatter(x,y, nodeSize, c="black",zorder=2)
    ax.axis("equal")
    ax.axis("off")
    return nothing
end

#plotNetwork!(ax, Q::Network) = plotNetwork!(ax, Q.r, Q.SM)

function plotNetwork(r::Vector, A::Matrix)
    fig, ax = subplots()
    plotNetwork!(ax, r, A)
    return fig
end

plotNetwork(tables::Tables) = plotNetwork(tables.R, tables.A)

function plotDisplacement!(ax, r, dr; zorder=1,
                         width=0.004, color="r", angles="xy", scale_units="xy", scale=0.5)
    xy = reshape(r, 2, :)'
    x,y = xy[:,1], xy[:,2]
    uv = reshape(dr, 2, :)'
   # uv = r2xy(dr)
    u,v = uv[:,1], uv[:,2]
    ax.quiver(x,y,u,v, zorder=zorder, width = width, color=color,
              scale=scale, angles=angles, scale_units=scale_units)
    return nothing
end

function plotNetworkColoredBonds!(ax, fig, R::Vector,
                                  A::Matrix,
                                  bondColorMat::Matrix;
                                  bondWidth=3,
                                  nodeSize=30,
                                  cmap="Purples")

    # plot a network, with the bond painted according to bondColorMat.
    # the floats in bondColorMat should be between 0 and 1.
    cm = get_cmap(cmap)
    #xy = r2xy(r)
    xy = reshape(R, 2, :)'
    x = xy[:,1]
    y = xy[:,2]
    N = size(A,1)
    for i in 2:N, j in 1:i-1
        if A[i,j] != 0
            bondColor = cm(bondColorMat[i,j])
            ax.plot([x[i],x[j]], [y[i],y[j]], lw=bondWidth, c=bondColor, zorder=0)
        end
    end
    ax.scatter(x,y, nodeSize, c="black",zorder=2)
    ax.axis("equal")
    clb = addCustomColorBar!(ax, fig; cmap)
    return clb
end

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


function plotBondStiffness!(ax, fig, net::Network, tables::Tables, SM::Matrix;
        bondWidth=3, nodeSize=30, cmap="Purples",
        k_max = log(10, maximum(tables.K)),
        k_min = log(10, minimum(tables.K)))
    
    A, K, R = tables.A, tables.K, tables.R
    seq = net.seq
    N = size(A, 1)

    bondColorMat = copy(SM)
   # k_min, k_max = log.(10, extrema(filter(x-> x!=0, SM)))
    # Normalize the bondColorMat to 0,1.
    for i in 2:N, j in 1:i-1
        if A[i,j] != 0
            bondColorMat[i,j] = bondColorMat[j,i] =  (log(10,SM[i,j]) - k_min) / (k_max - k_min)
        end
    end
    clb = plotNetworkColoredBonds!(ax, fig, R, A, bondColorMat; bondWidth, nodeSize, cmap)
    clb.set_label("Stiffness,  " *L"\log_{10} k", labelpad=-10, fontsize=16)
    clb.set_ticks([0,1])
    clb.set_ticklabels([round(k_min), round(k_max)])
    ax.axis("equal")
    ax.axis("off")
    return nothing
    plotBondStiffness!(ax, fig, Q.r, Q.SM, Q.SM; bondWidth, nodeSizem, cmap)
    return nothing
end


function plotBondStiffness!(ax, fig,
                            net::Network,
                            tables::Tables,
                            bondWidth=3,
                            nodeSize=30,
                            cmap="Purples",
                            k_max = log(10, maximum(tables.K)),
                            k_min = log(10, minimum(tables.K)))

    A, K, R = tables.A, tables.K, tables.R
    seq = net.seq
    N = size(A, 1)
    SM = zeros(N,N)
    for i in 1:N, j in 1:N
        SM[i,j] = tables.A[i,j] * tables.K[seq[i], seq[j]]
    end
    plotBondStiffness!(ax, fig, net, tables,  SM; bondWidth, nodeSize, cmap, k_max, k_min)
    return nothing
end




function addCustomColorBar!(ax::PyObject,
                            fig::Figure;
                            cmap="inferno",
                            colorRange=(0,1),
                            width=1,
                            height=0.8)

    # add a colobar for a figure without changing the size of the axes
    # and without needing a mappable object. Use this function when
    # you know the smallest and largest values that will be colored
    # and 

    norm = matplotlib.colors.Normalize(vmin=colorRange[1], vmax=colorRange[2]) 
    sm = matplotlib.cm.ScalarMappable(cmap=get_cmap(cmap), norm=norm)
    bb = ax.get_position()
    x0, y0, h, w = bb.x0, bb.y0, bb.height, bb.width
    
    final_x0 = x0+w + w*0.01
    final_y0 =  y0 + bb.height * (1-height) / 2
    final_w = w*0.05 * width
    final_h = bb.height * height

    cax = fig.add_axes([final_x0, final_y0, final_w, final_h])
    clb = fig.colorbar(sm, cax=cax) 
    return clb
end

# ------------------------------------------------------------------------


"""
Plot the binding free energies points in negative binding space for many networks.

ΔF stores binding energies to ligand 1 in col 1 and ligand 2 in col 2.
"""
function plotDensityBindingSpaceScatter(ΔF::Matrix; buffer=0.04) 
    fig, ax = subplots()
    mm,MM = extrema( [-ΔF[:]..., 0.0] )
    ax.scatter(-ΔF[:,1], -ΔF[:,2], 1, "k")
    ax.set_xlim([mm-buffer, MM+buffer])
    ax.set_ylim([mm-buffer, MM+buffer])
    ax.set_xlabel(L"-\Delta F_1", fontsize=16)
    ax.set_ylabel(L"-\Delta F_2", fontsize=16)
    ax.axvline(x=0, color="k", linestyle="--")
    ax.axhline(y=0, color="k", linestyle="--")
    ax.set_aspect("equal")
    return fig
end



"""
Plot population average negative binding energies over time.
ΔF1 and ΔF2 should be binding energies, ie negative for tight binding
"""
function plotBindingEvoTraces(ΔF1::Vector,ΔF2::Vector; 
                              τ=0, μ=0, P=0, selParams=(0,0))
    fig, ax = subplots( figsize=(10,2))
    ax.plot(-ΔF1, lw=1, "r", label="Lig 1")
    ax.plot(-ΔF2, lw=1, "b", label="Lig 2")
    ax.set_xlabel("generations", fontsize=14)
    ax.set_ylabel(L"\langle -\Delta F \rangle_{\rm pop}", fontsize=14)
    ax.set_xlim([0, length(ΔF1)])
    ax.legend()
    ax.axhline(y=0, color="grey", linestyle="--", linewidth=0.8, zorder=0)
    ax.set_title(L"\tau=" * @sprintf("%.3g", τ) * ", " *
             L"\mu=" * @sprintf("%.3g", μ) * ", " *
             L"P=" * @sprintf("%.3g", P) * ", " *
             L"a=" * @sprintf("%.3g", selParams[1]) * ", " *
             L"x_0=" * @sprintf("%.3g", selParams[2]))
    return fig
end

"""
Plot trajectory of population mean negative binding energy for an evolutionary trajectory. 
ΔF1 and ΔF2 should be binding energies, ie negative for tight binding
"""
function plotEvoPathBindingSpace(ΔF1::Vector,ΔF2::Vector; lw=0.5)
    t = 1:length(ΔF1)  # Time or index
    norm = (t .- minimum(t)) ./ (maximum(t) - minimum(t))
    cmap = PyPlot.get_cmap("viridis")
    fig, ax = subplots()
    for i in 1:length(t)-1
        color = cmap(norm[i])  # Get the RGBA color for the current segment
        ax.plot(-ΔF1[i:i+1], -ΔF2[i:i+1], lw=lw, color=color)
    end
    mm, MM = extrema([-ΔF1[:]..., -ΔF2[:]..., 0.0])
    ax.set_xlim([mm - 0.01, MM+0.01])
    ax.set_ylim([mm - 0.01, MM+0.01])
    ax.set_xlabel(L"-\Delta F_1", fontsize=16)
    ax.set_ylabel(L"-\Delta F_2", fontsize=16)
    ax.axvline(x=0, lw=0.5, color="k", linestyle="--", zorder=0)
    ax.axhline(y=0, lw=0.5, color="k", linestyle="--", zorder=0)
    ax.set_aspect("equal")
    return fig
end


"""
Make a png showing population's distribution in negative binding space.
save to file.
"""
function makePopulationFrame(fits::Matrix, i, k, dir_path; scale = 0.05)
    @assert dir_path[end] == '/'
    # for i in 1:N_frames
    fig, ax = subplots()
    #ax = axs[i,j]
    ax.scatter(fits[:,1], fits[:,2], 5, color="b", alpha=0.1)
    ax.set_aspect("equal")
    ax.set_xlim([-scale, scale])
    ax.set_ylim([-scale, scale])
    ax.axvline(x=0, lw=0.5, color="k", linestyle="--", zorder=0)
    ax.axhline(y=0, lw=0.5, color="k", linestyle="--", zorder=0)
    ax.axis("off")
    ax.text(-0.9scale, 0.9scale, "t="*"$i")
    ax.text(0.8scale, -0.1scale, L"-\Delta F_1")
    ax.text(-0.1scale, 0.8scale, L"-\Delta F_2", rotation=90)

    if k == 1
        ax.annotate("", xy=(0.7 * scale, -0.9 * scale), xytext=(0.5 * scale, -0.9 * scale),
            arrowprops=Dict(
                "facecolor" => "red",
                "edgecolor" => "red",
                "width" => 0.4,
                "headwidth" => 6))
    elseif k == 2
        ax.annotate("", xy=(-0.9 * scale, 0.5 * scale), xytext=(-0.9 * scale, 0.3 * scale),
            arrowprops=Dict(
                "facecolor" => "red",
                "edgecolor" => "red",
                "width" => 0.4,
                "headwidth" => 6))
    end
    fig.savefig(dir_path*@sprintf("%04d.png", i), bbox_inches="tight", dpi=500)
    close(fig)
    PyPlot.matplotlib.pyplot.close("all")
    return nothing
end


"""
make movie of the population dynamics in negative binding space.
"""
function makePopulationDynamicsMovieFrames(fitBucket::Array{T, 3},
                                           dir_path::String,
                                           τ::Number;
                                           N_frames = size(fitBucket, 3),
                                           numEnv=2,
                                           k_init=2) where T
    k = k_init
    mx = maximum(abs, fitBucket)
    fits = [fitBucket[:,:,i] for i in 1:N_frames]
    envs = zeros(Int, N_frames)
    for i in 1:N_frames
        envs[i] = k
        i % τ == 0 && (k = k % numEnv + 1)
    end
    
    foreach( (x,y,z) -> makePopulationFrame(x,y,z, dir_path; scale=mx), fits, 1:length(fits), envs) 


    # ffmpeg -framerate 5 -i %04d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p populations.mp4
    run(`ffmpeg -framerate 5 -i $(dir_path)%04d.png $(dir_path)populations.mp4 -y`)
    run(`open $(dir_path)populations.mp4`)
    return nothing
end


"""
Plot the number of dead, specialist and generalist sequences over time for
an evolutionary trajectory.

fits[:,:,i] should show the population negative binding energies at generation i.
"""
function plotGeneralistSpecialistEvoTraces(fits::AbstractArray{T,3}) where T


    function binByQuadrants(fits::Array{T,3}) where T
        # Classify each sequence as dead specialist, or generalist
        #counts[1,1] is the number of "dead" sequences
        #counts[2,1] is the number of ligand 1 specialist sequences
        #counts[1,2] is the number of ligand 1 specialist sequences
        #counts[2,2] is the number of generalist sequences

        function kernel(fits::Matrix)
            counts = zeros(2,2)
            for i in axes(fits, 1)
                x, y = fits[i,:] .> 0
                counts[x+1, y+1] += 1
            end
            return counts
        end

        @assert size(fits, 2) == 2
        fullCounts = zeros(2,2,size(fits,3))
        for i in 1:size(fits,3)
           fullCounts[:,:,i] = kernel(fits[:,:,i])
        end
        return fullCounts
    end
    counts = binByQuadrants(fits);
    fig, ax = subplots( figsize=(10,2))
    ax.plot(counts[1,1,:], lw=0.5, color="grey")
    ax.plot(counts[2,1,:], lw=0.5, color="red")
    ax.plot(counts[1,2,:], lw=0.5, color="blue")
    ax.plot(counts[2,2,:], lw=1, color="green")
    #ax.set_xlabel("generations", fontsize=14)
    #ax.set_ylabel(L"\langle -\Delta F \rangle_{\rm pop}", fontsize=14)
    ax.set_xlim([0, size(fits, 3)])
    #ax.legend()
    ax.axhline(y=0, color="grey", linestyle="--", linewidth=0.8, zorder=0)
    return fig, counts
end

"""
plot all single mutants in binding space.
"""
function plotDMS(dms::Array{T,3}) where T <: Number
    xy = -reshape(dms, :, 2)
    fig, ax = subplots()
    ax.axis("square")  # Equal aspect ratio
    ax.set_xlabel(L"\Delta (-\Delta F_1)")
    ax.set_ylabel(L"\Delta (-\Delta F_2)")
    lim = 1.1 * maximum(abs, xy)  # Symmetric limits
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.scatter(xy[:, 1], xy[:, 2], s=2, alpha=0.3, color="k")
    ax.scatter(0, 0, s=20, color="red", zorder=3)
    return fig  # Return figure and axes for further customization
end

function plotMutSensVsDist!(ax, dms::Array{T,3}, siteBeads::Vector, tables::Tables; color="k") where T <: Number
    dmss = abs.(vcat([dms[:, :, i] for i in 1:size(dms, 3)]...))
    stepsVec = getStepsVector(tables.A, siteBeads);
    x = repeat(stepsVec, size(dmss, 1))
    y = log.(10, dmss'[:] .+ eps())
    inds = y .> -14
    y = y[inds]
    x = x[inds]
    ax.scatter(x, y, 15, color=color, alpha=0.1)
    ax.set_ylabel("log10 abs mut sensitivity", fontsize=14)
    ax.set_xlabel("distance from active site",fontsize=14)
    ax.set_ylim([minimum(y), maximum(y) + 0.5])
    return nothing
end

function plotMutSensVsDist(dms::Array{T,3}, siteBeads::Vector, tables::Tables; color="k") where T <: Number
    fig, ax = subplots(figsize=(4,4))
    plotMutSensVsDist!(ax, dms, siteBeads, tables; color)
    return fig
end



















