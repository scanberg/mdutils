#include "mol/hydrogen_bond.h"
#include <core/common.h>
#include <core/log.h>
#include <core/spatial_hash.h>
#include <mol/trajectory_utils.h>

namespace hydrogen_bond {

// Computes the potential donors given a set of atom labels.
// OH and NH atoms are assumed to be donors if the concecutive atom is marked with 'H' for Hydrogen.
i32 compute_donors(DynamicArray<HydrogenBondDonor>* donors, Array<const Element> elements, Array<const ResIdx> residue_indices, Array<const Residue> residues, Array<const Bond> covalent_bonds) {
    ASSERT(donors);
    i32 pre_count = (i32)donors->size();
    for (i32 i = 0; i < (i32)elements.size(); i++) {
        if (elements[i] == Element::H) {
            // get all bonds with atom i
            const auto& res = residues[residue_indices[i]];
            for (const auto& bond : covalent_bonds.subarray(res.bond_idx.beg, res.bond_idx.end - res.bond_idx.beg)) {
                if (i == bond.idx[0] || i == bond.idx[1]) {
                    const i32 j = bond.idx[0] != i ? bond.idx[0] : bond.idx[1];
                    const auto elem = elements[j];
                    if (elem == Element::O || elem == Element::N || elem == Element::F) {
                        donors->push_back({j, i});
                        break;
                    }
                }
            }
        }
    }

    /*
ASSERT(donors);
int32 pre_count = (int32)donors->count;
const int32 num_labels = (int32)labels.count;
for (int32 i = 0; i < num_labels; i++) {
    if (compare_n(labels[i], "OH", 2) || compare_n(labels[i], "NH", 2)) {
        if (i + 1 < num_labels && compare_n(labels[i + 1], "H", 1)) {
            donors->push_back({i, i + 1});
        }
        if (i + 2 < num_labels && compare_n(labels[i + 2], "H", 1)) {
            donors->push_back({i, i + 2});
        }
    }
}
    */

    return (i32)donors->size() - pre_count;
}

DynamicArray<HydrogenBondDonor> compute_donors(Array<const Element> elements, Array<const ResIdx> residue_indices, Array<const Residue> residues, Array<const Bond> covalent_bonds) {
    DynamicArray<HydrogenBondDonor> donors;
    compute_donors(&donors, elements, residue_indices, residues, covalent_bonds);
    return donors;
}

// Computes the potential acceptors given a set of atom elements.
// This essentially just a filter on atom element which extracts Oxygen and Nitrogen
i32 compute_acceptors(DynamicArray<HydrogenBondAcceptor>* acceptors, Array<const Element> elements) {
    ASSERT(acceptors);
    const i32 pre_count = (i32)acceptors->size();
    for (i32 i = 0; i < (i32)elements.count; i++) {
        if (elements[i] == Element::O || elements[i] == Element::N || elements[i] == Element::F) {
            acceptors->push_back(i);
        }
    }
    return (i32)acceptors->size() - pre_count;
}

DynamicArray<HydrogenBondAcceptor> compute_acceptors(Array<const Element> elements) {
    DynamicArray<HydrogenBondAcceptor> acceptors;
    compute_acceptors(&acceptors, elements);
    return acceptors;
}

// Computes hydrogen bonds given a certain set of potential donors, acceptors and atomic positions from a frame.
// The distance cutoff sets the distance from bonds to potential acceptors.
//

i32 compute_bonds(DynamicArray<HydrogenBond>* bonds, Array<const HydrogenBondDonor> donors, Array<const HydrogenBondAcceptor> acceptors, const float* atom_pos_x, const float* atom_pos_y,
                    const float* atom_pos_z, float dist_cutoff, float angle_cutoff) {
    const i32 num_acceptors = (i32)acceptors.count;
    if (!num_acceptors) return 0;

    DynamicArray<float> acceptor_pos_x(num_acceptors);
    DynamicArray<float> acceptor_pos_y(num_acceptors);
    DynamicArray<float> acceptor_pos_z(num_acceptors);
    DynamicArray<AtomIdx> acceptor_idx(num_acceptors);

    for (i32 i = 0; i < num_acceptors; i++) {
        acceptor_pos_x[i] = atom_pos_x[acceptors[i]];
        acceptor_pos_y[i] = atom_pos_y[acceptors[i]];
        acceptor_pos_z[i] = atom_pos_z[acceptors[i]];
        acceptor_idx[i] = acceptors[i];
    }

    i32 pre_count = (i32)bonds->size();
    spatialhash::Frame frame = spatialhash::compute_frame(acceptor_pos_x.data(), acceptor_pos_y.data(), acceptor_pos_z.data(), num_acceptors, vec3(dist_cutoff));
    for (const auto& don : donors) {
        vec3 donor_pos_xyz = {atom_pos_x[don.donor_idx], atom_pos_y[don.donor_idx], atom_pos_z[don.donor_idx]};
        vec3 hydro_pos_xyz = {atom_pos_x[don.hydro_idx], atom_pos_y[don.hydro_idx], atom_pos_z[don.hydro_idx]};
        spatialhash::for_each_within(frame, hydro_pos_xyz, dist_cutoff, [bonds, &donor_pos_xyz, &hydro_pos_xyz, &acceptor_idx, &don, angle_cutoff](i32 idx, const vec3& pos) {
            AtomIdx g_idx = acceptor_idx[idx];
            if (g_idx == don.donor_idx) return;
            const vec3 a = hydro_pos_xyz - donor_pos_xyz;
            const vec3 b = pos - hydro_pos_xyz;
            if (math::angle(a, b) < angle_cutoff) {
                bonds->push_back({g_idx, don.donor_idx, don.hydro_idx});
            }
        });
    }
    return (i32)bonds->size() - pre_count;
}

DynamicArray<HydrogenBond> compute_bonds(Array<const HydrogenBondDonor> donors, Array<const HydrogenBondAcceptor> acceptors, const float* atom_pos_x, const float* atom_pos_y, const float* atom_pos_z,
                                         float dist_cutoff, float angle_cutoff) {
    DynamicArray<HydrogenBond> bonds;
    compute_bonds(&bonds, donors, acceptors, atom_pos_x, atom_pos_y, atom_pos_z, dist_cutoff, angle_cutoff);
    return bonds;
}

/*
void compute_bonds_trajectory(HydrogenBondTrajectory* hbt, const MoleculeDynamic& dyn, float max_dist, float max_angle) {
    ASSERT(hbt);
    for (int32 i = 0; i < dyn.trajectory.num_frames; i++) {
        Array<HydrogenBond> frame_bonds{(HydrogenBond*)(hbt->bond_data.end()), int64(0)};
        Array<const vec3> atom_positions = get_trajectory_positions(dyn.trajectory, i);
        frame_bonds.count = compute_bonds(&hbt->bond_data, dyn.molecule.hydrogen_bond.donors, dyn.molecule.hydrogen_bond.acceptors, atom_positions, max_dist, max_angle);
    }
}

HydrogenBondTrajectory compute_bonds_trajectory(const MoleculeDynamic& dyn, float max_dist, float max_angle) {
    HydrogenBondTrajectory hbt;
    compute_bonds_trajectory(&hbt, dyn, max_dist, max_angle);
    return hbt;
}
*/

}  // namespace hydrogen_bond
