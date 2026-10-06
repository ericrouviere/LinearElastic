
# Selective-pressure scoring functions on raw ligand free energies.
# F0 is the solvent (unbound) free energy; the rest are bound states.
binding(F0, F1) = F0 - F1 # strong binders > 0
doubleBinding(F0, F1, F2) = min(binding(F0, F1), binding(F0, F2))
allostery(E00, E10, E01, E11) = E10 + E01 - E11 - E00

function targetBinding(ΔF1::Number, ΔF2::Number, ΔF1_target::Number, ΔF2_target::Number)
    # return a fitness function that increased (to zero) as
    # the binding energies approach the target binding energies.
    err1 = abs(ΔF1 - ΔF1_target)
    err2 = abs(ΔF2 - ΔF2_target)
    return -max(err1, err2)
end


# One immutable type per selective pressure. Selection dispatches on the assay
# instance instead of branching on a string, so adding a pressure is just a new
# struct plus a `computeFitness` method (no central if/elseif to edit).
abstract type Assay end

"Binding to one ligand: `F0 - F1`. Needs 2 ligand states (solvent, ligand)."
struct Binding <: Assay end

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
computeFitness(::Binding, energies) = binding(energies...)
computeFitness(::DoubleBinding, energies) = doubleBinding(energies...)
computeFitness(::Allostery, energies) = allostery(energies...)
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
