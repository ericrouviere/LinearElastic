
"""
    buildStructure(W, L, offset)

Bead positions and adjacency matrix of a triangular lattice `W` beads wide and `L`
layers long (`L` must be odd). Beads less than 1.1 lattice spacings apart are bonded;
then every coordinate is displaced by Gaussian noise of standard deviation `offset`.
Returns `R = [x1, y1, x2, y2, ...]` and the adjacency matrix `A`.
"""
function buildStructure(W::Int, L::Int, offset::Number)
    @assert isodd(L)
    numLayers = Int((L-1)/2)
    stride = 2W-1
    N_beads = stride * numLayers + W - 1
    v1 = [cos(pi/6), sin(pi/6)] # basis vector 1 for triangular matrix
    v2 = [0,1] # basis vector 2 for triangular matrix
    V = [v1 v2] # change of basis matrix. From lattice basis to standard basis
    layer = (V * [zeros(W-1) 0:W-2 ; ones(W) -1:W-2]')'
    layer_x = layer[:,1]
    layer_y = layer[:,2]

    px = zeros(stride*numLayers + W-1 )
    py = zeros(stride*numLayers + W -1)
    for i in 1:numLayers    
        px[(i-1)*stride+1:(i)*stride] .= layer_x .+ (i-1)*2*sin(pi/3)
        py[(i-1)*stride+1:(i)*stride] .= layer_y
    end
    # add last layer
    px[end-W+2:end] .= layer_x[1:W-1] .+ (numLayers)*2*sin(pi/3)
    py[end-W+2:end] .= layer_y[1:W-1]

    # build Adjacency Matrix
    A = zeros(N_beads,N_beads);
    for i in 2:N_beads
        for j in 1:i-1
            dx = px[i] - px[j]
            dy = py[i] - py[j]
            dist = sqrt(dx^2 + dy^2)
            if dist < 1.1
                A[i,j] = A[j,i] = 1.
            end
        end
    end
    x = copy(px) .+ offset .* randn(N_beads)
    y = copy(py) .+ offset .* randn(N_beads)
    R = reshape([x y]', :)
    return R, A
end


"""
    buildSpringTable(N_types, k_min, k_max; binary=false, N_soft=0)

Symmetric `N_types × N_types` table of bond stiffnesses: `K[a, b]` is the spring
constant of a bond between beads of types `a` and `b`. Entries are log-uniform and
rescaled to span exactly `[k_min, k_max]`. With `binary=true`, all entries are `k_max`
except `N_soft` random entries set to `k_min`.
"""
function buildSpringTable(N_types::Int, k_min::Number, k_max::Number; binary=false, N_soft=0)

    if binary # make all but N_soft entries stiff
        K = k_max*ones(N_types,N_types)
        inds = findall(LowerTriangular(K) .> 0)
        rp = randperm(length(inds))
        softInds = inds[rp[1:N_soft]]
        K[softInds] .= k_min
    else # pull K from log uniform distribution
        log_diff = log(10, k_max) - log(10, k_min)
        randos = rand(N_types, N_types)
        randos =  Float64.(Symmetric(randos, :L))
        m, M = extrema(randos)
        randos .-= m
        randos ./= (M-m)
        K = randos
        K = 10 .^ (randos .* log_diff .+ log(10,k_min))
        # The Old Stuff ---------------------------------
        # log_diff = log(10, k_max) - log(10, k_min)
        # randos = rand(N_types, N_types)
        # K = 10 .^ (randos .* log_diff .+ log(10,k_min))
        # -----------------------------------------------
    end
    K = Float64.(Symmetric(K, :L))
    return K
end

"""
Symmetric table `E[a, b]` of the contact energy of a bond between beads of types `a`
and `b`, Gaussian with standard deviation `energyVar`. Each bond adds its contact
energy to the network's energy, independent of the displacement.
"""
function buildContactEnergyTable(N_types::Int, energyVar::Number)
    return Float64.(Symmetric( energyVar*randn(N_types, N_types) ))
end

"""
Symmetric table `D[a, b]` of the rest-length mismatch of a bond between beads of
types `a` and `b`, Gaussian with standard deviation `mismatchVar`.
"""
function buildMismatchTable(N_types::Int, mismatchVar::Number)
    return Float64.(Symmetric( mismatchVar*randn(N_types, N_types) ))
end


"""
    makeTables(W, L, offset, N_types, k_min, k_max, energyVar, mismatchVar; binary=false, N_soft=0)

Build the `Tables` of an experiment: a perturbed `W × L` triangular lattice
(`buildStructure`) and random bead-type tables of stiffness (`buildSpringTable`),
mismatch (`buildMismatchTable`) and contact energy (`buildContactEnergyTable`).
Draws from the global RNG; seed it first for a reproducible network.

To use chosen values for a table instead, build the struct directly, e.g.
`Tables(t.R, t.A, t.K, t.D, E)` for a symmetric `N_types × N_types` contact-energy
matrix `E`.
"""
function makeTables(W::Int,
                   L::Int,
                   offset::Number,
                   N_types::Int,
                   k_min::Number,
                   k_max::Number,
                   energyVar::Number,
                   mismatchVar::Number;
                   binary=false,
                   N_soft=0)
    R, A = buildStructure(W, L, offset)
    K = buildSpringTable(N_types, k_min, k_max; binary, N_soft)
    D = buildMismatchTable(N_types, mismatchVar)
    E = buildContactEnergyTable(N_types, energyVar)
    return Tables(R, A, K, D, E)
end

"""
Uniformly random sequence of bead types: one type in `1:size(K, 1)` for each bead
of the positions `R`.
"""
function randSeq(R, K)
    N_beads = Int(length(R)/2)
    N_types = size(K, 1)
    seq = rand(1:N_types, N_beads)
    return seq
end

"""
Build the `Network` of sequence `seq`, computing its energy parameters from scratch.
"""
function buildNetwork(seq::Vector, tables::Tables)
    G, lpdH, h, c = getNetParams(seq, tables)
    return Network(seq, G, similar(G), lpdH, 0.0, h, similar(h), c, 0.0)
end

"""
Build the `Network` of a uniformly random sequence.
"""
function buildRandNetwork(tables::Tables)
    seq = randSeq(tables.R, tables.K)
    return buildNetwork(seq, tables)
end
