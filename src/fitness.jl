
# Selective-pressure scoring functions on raw ligand free energies.
# F0 is the solvent (unbound) free energy; the rest are bound states.
"""
Binding score `F0 - F1` from the solvent and bound free energies (positive for a binder).
"""
binding(F0, F1) = F0 - F1 # strong binders > 0

"""
Generalist score `min(F0 - F1, F0 - F2)`: binding to both ligands.
"""
doubleBinding(F0, F1, F2) = min(binding(F0, F1), binding(F0, F2))

"""
Allosteric cooperativity `E10 + E01 - E11 - E00` from the free energies of the
unbound, singly bound and doubly bound states.
"""
allostery(E00, E10, E01, E11) = E10 + E01 - E11 - E00

"""
Fitness `-max(|ΔF1 - ΔF1_target|, |ΔF2 - ΔF2_target|)`, which rises to zero as both
binding energies approach their targets.
"""
function targetBinding(ΔF1::Number, ΔF2::Number, ΔF1_target::Number, ΔF2_target::Number)
    err1 = abs(ΔF1 - ΔF1_target)
    err2 = abs(ΔF2 - ΔF2_target)
    return -max(err1, err2)
end


# One immutable type per selective pressure. Selection dispatches on the assay
# instance instead of branching on a string, so adding a pressure is just a new
# struct plus a `computeFitness` method (no central if/elseif to edit).
"""
Selective pressure scored by `computeFitness`. Subtypes are passed as values, e.g.
`DoubleBinding()`.
"""
abstract type Assay end

"Binding to one ligand: `F0 - F1`. Needs 2 ligand states (solvent, ligand)."
struct Binding <: Assay end

"Binding to ligand 1, `F0 - F1`, from three ligand states (solvent, ligand 1, ligand 2)."
struct Binding1 <: Assay end

"Binding to ligand 2, `F0 - F2`, from three ligand states (solvent, ligand 1, ligand 2)."
struct Binding2 <: Assay end

"""
Simultaneous binding to two ligands, the generalist score `min(F0 - F1, F0 - F2)`.
Needs 3 ligand states (solvent, ligand 1, ligand 2).
"""
struct DoubleBinding <: Assay end

"""
Allosteric cooperativity `E10 + E01 - E11 - E00` between two sites.
Needs 4 ligand states (unbound, site 1, site 2, both).
"""
struct Allostery <: Assay end

"""
    TargetBinding(ΔF1_target, ΔF2_target)

Fitness rises to zero as the two binding energies approach the targets.
Needs 3 ligand states (solvent, ligand 1, ligand 2).
"""
struct TargetBinding <: Assay
    ΔF1_target::Float64
    ΔF2_target::Float64
end

# Scalar fitness from precomputed ligand free energies, dispatched on the assay.
"""
Fitness under `Binding` from the free energies `[F0, F1]`.
"""
computeFitness(::Binding, energies) = binding(energies...)

"""
Fitness under `Binding1` from the free energies `[F0, F1, F2]`.
"""
computeFitness(::Binding1, energies) = binding(energies[1], energies[2])

"""
Fitness under `Binding2` from the free energies `[F0, F1, F2]`.
"""
computeFitness(::Binding2, energies) = binding(energies[1], energies[3])

"""
Fitness under `DoubleBinding` from the free energies `[F0, F1, F2]`.
"""
computeFitness(::DoubleBinding, energies) = doubleBinding(energies...)

"""
Fitness under `Allostery` from the free energies `[E00, E10, E01, E11]`.
"""
computeFitness(::Allostery, energies) = allostery(energies...)

"""
Fitness under `TargetBinding` from the free energies `[F0, F1, F2]`.
"""
computeFitness(assay::TargetBinding, energies) =
    targetBinding(energies[2] - energies[1], energies[3] - energies[1],
                  assay.ΔF1_target, assay.ΔF2_target)

"""
    computeFitness(net, ligands, tables, assay::Assay, β)

Scalar fitness of `net` under selective pressure `assay` (e.g. `DoubleBinding()`).
Computes the ligand free energies and dispatches on `assay` to the matching
scoring function.
"""
function computeFitness(net::Network, ligands::Ligands, tables::Tables,
                        assay::Assay, β::Number)
    return computeFitness(assay, computeFreeEnergies(net, ligands, tables, β))
end


# Folding. A sequence folds if the free energy F of its network in the solvent state
# (λ0) is at most the threshold F_unfolded, and is unfolded otherwise.

"""
    isFolded(F, F_unfolded)

Whether a protein whose network has free energy `F` in the solvent state is folded:
`F <= F_unfolded`, a hard stability threshold. Always true for `F_unfolded = Inf`.
In population evolution an unfolded individual has selection weight 0: it cannot
reproduce, whatever its binding energy.
"""
isFolded(F::Real, F_unfolded::Real) = F <= F_unfolded
