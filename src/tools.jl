"""
    convertEvoParams(S::Dict)

Build an `EvoParams` from a settings dict with keys "P", "N", "μ", "τ", "α" and,
optionally, "F_unfolded" (default `Inf`, as in settings written before it existed).
"""
function convertEvoParams(S::Dict)
    return EvoParams(S["P"], S["N"], Float64(S["μ"]), S["τ"], Float64(S["α"]),
                     Float64(get(S, "F_unfolded", Inf)))
end

# convert between postions as a vector r = [x1,y1,x2,y2,...]
# to position as a matrix xy = [x1 y1; x1 y1;..]

"""
Copy the sequence and cached parameters of `net_src` into `net_des` without
allocating (the scratch fields are not copied).
"""
function networkCopy!(net_des::Network, net_src::Network)
    copyto!(net_des.seq, net_src.seq)
    copyto!(net_des.G, net_src.G)
    copyto!(net_des.h, net_src.h)
    copyto!(net_des.u, net_src.u)
    net_des.lpdH = net_src.lpdH
    net_des.c = net_src.c
    net_des.hGh = net_src.hGh
    return nothing
end

