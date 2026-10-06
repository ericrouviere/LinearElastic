

function computeFitness(seq, landscape)
    J, h = landscape.J, landscape.h
    N = length(seq)
    energy = 0.0

    # Field contribution
    for i in 1:N
        energy -= h[seq[i], i]
    end

    # Pairwise coupling contribution
    for i in 1:N-1
        for j in i+1:N
            energy -= J[seq[i], seq[j], i, j]
        end
    end
    return energy
end

