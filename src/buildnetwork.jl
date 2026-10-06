
function buildStructure(W::Int, L::Int, offset::Number)
    # return a structure `r` and Adacency matrix `A` 
    # offset controlles how much points are perturbed of tri lattice.    
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


function buildSpringTable(N_types::Int, k_min::Number, k_max::Number; binary=false, N_soft=0)
    # generate a Matrix K_ij where i and j are nodes types and
    # and K_ij is the spring constant that connect nodes of type i
    # and type j.

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

function buildContactEnergyTable(N_types::Int, energyVar::Number)
    # Generate a Matrix E_ij where i and j are bead types and
    # and E_ij is the contact energy that connect nodes of type i and j
    return Float64.(Symmetric( energyVar*randn(N_types, N_types) ))
end

function buildMismatchTable(N_types::Int, mismatchVar::Number)
    # Generate a Matrix D_ij where i and j are bead types and
    # and D_ij is the contact energy that connect nodes of type i and j
    return Float64.(Symmetric( mismatchVar*randn(N_types, N_types) ))
end


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
    # make a table object.
    R, A = buildStructure(W, L, offset)
    K = buildSpringTable(N_types, k_min, k_max; binary, N_soft)
    D = buildMismatchTable(N_types, mismatchVar)
    E = buildContactEnergyTable(N_types, energyVar)
    return Tables(R, A, K, D, E)
end

function randSeq(R, K)
    # generate a random sequence of node types.
    N_beads = Int(length(R)/2)
    N_types = size(K, 1)
    seq = rand(1:N_types, N_beads)
    return seq
end

function buildNetwork(seq::Vector, tables::Tables)
    # return an updated network struct
    G, lpdH, h, c = getNetParams(seq, tables)
    return Network(seq, G, similar(G), lpdH, 0.0, h, similar(h), c, 0.0)
end

function buildRandNetwork(tables::Tables)
    seq = randSeq(tables.R, tables.K)
    return buildNetwork(seq, tables)
end
