"""
    convertEvoParams(S::Dict)

Build an `EvoParams` from a settings dict with keys "P", "N", "μ", "τ" and "α".
"""
function convertEvoParams(S::Dict)
    return EvoParams(S["P"], S["N"], Float64(S["μ"]), S["τ"], Float64(S["α"]))
end

# convert between postions as a vector r = [x1,y1,x2,y2,...]
# to position as a matrix xy = [x1 y1; x1 y1;..]

function networkCopy!(net_des::Network, net_src::Network)
    # copy one network into another without allocating memory.  
    copyto!(net_des.seq, net_src.seq)
    copyto!(net_des.G, net_src.G)
    copyto!(net_des.h, net_src.h)
    net_des.lpdH = net_src.lpdH
    net_des.c = net_src.c
    return nothing
end

