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
};

struct HydrogenBondTrajectory {
    DynamicArray<HydrogenBond> bond_data{};
    DynamicArray<Array<HydrogenBond>> frame_bonds{};
};

namespace hydrogen_bond {
i32 compute_acceptors(DynamicArray<HydrogenBondAcceptor>* acceptors, Array<const Element> elements);
DynamicArray<HydrogenBondAcceptor> compute_acceptors(Array<const Element> elements);

i32 compute_donors(DynamicArray<HydrogenBondAcceptor>* donors, Array<const Element> elements, Array<const ResIdx> residue_indices, Array<const Residue> residues, Array<const Bond> covalent_bonds);
DynamicArray<HydrogenBondDonor> compute_donors(Array<const Element> elements, Array<const ResIdx> residue_indices, Array<const Residue> residues, Array<const Bond> covalent_bonds);

i32 compute_bonds(DynamicArray<HydrogenBond>* bonds, Array<const HydrogenBondDonor> donors, Array<const HydrogenBondAcceptor> acceptors, const float* atom_pos_x, const float* atom_pos_y,
                    const float* atom_pos_z, float dist_cutoff = 3.f, float angle_cutoff = math::deg_to_rad(20.f));
DynamicArray<HydrogenBond> compute_bonds(Array<const HydrogenBondDonor> donors, Array<const HydrogenBondAcceptor> acceptors, const float* atom_pos_x, const float* atom_pos_y, const float* atom_pos_z,
                                         float dist_cutoff = 3.f, float angle_cutoff = math::deg_to_rad(20.f));

/*
void compute_bonds_trajectory(HydrogenBondTrajectory* hbt, const MoleculeDynamic& dyn, float dist_cutoff, float angle_cutoff);
HydrogenBondTrajectory compute_bonds_trajectory(const MoleculeDynamic& dyn, float dist_cutoff, float angle_cutoff);
*/

}  // namespace hydrogen_bond
