#pragma once

#include <core/array_types.h>
#include <core/math_utils.h>
#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/molecule_dynamic.h>

struct HydrogenBond {
    AtomIdx acc_idx = 0;
    AtomIdx don_idx = 0;
    AtomIdx hyd_idx = 0;
    float strength  = 0;
};

namespace hydrogen_bond {
DynamicArray<HydrogenBondAcceptor> compute_acceptors(const Element in_elements[], i64 count);
DynamicArray<HydrogenBondDonor>    compute_donors(const MoleculeStructure& mol);
DynamicArray<HydrogenBond>         compute_bonds(Array<const HydrogenBondDonor> in_donors, Array<const HydrogenBondAcceptor> in_acceptors, const soa_vec3 in_pos,
                                         float dist_cutoff = 3.f, float angle_cutoff = math::deg_to_rad(20.f));

/*
void compute_bonds_trajectory(HydrogenBondTrajectory* hbt, const MoleculeDynamic& dyn, float dist_cutoff, float angle_cutoff);
HydrogenBondTrajectory compute_bonds_trajectory(const MoleculeDynamic& dyn, float dist_cutoff, float angle_cutoff);
*/

}  // namespace hydrogen_bond
