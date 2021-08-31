#include "mol/hydrogen_bond.h"
#include <core/common.h>
#include <core/log.h>
#include <core/spatial_hash.h>
#include <mol/trajectory_utils.h>

namespace hydrogen_bond {

// Computes the potential donors given a set of atom labels.
// OH and NH atoms are assumed to be donors if the concecutive atom is marked with 'H' for Hydrogen.
DynamicArray<HydrogenBondDonor> compute_donors(const MoleculeStructure& mol) {
    DynamicArray<HydrogenBondDonor> donors;
    for (i32 i = 0; i < (i32)mol.atom.count; i++) {
        if (mol.atom.element[i] == Element::H) {
            // get all bonds with atom i
            const auto& res_bond_range = mol.residue.bond.complete[mol.atom.res_idx[i]];
            const Array<Bond> res_bonds = {mol.covalent_bond.bond + res_bond_range.beg, res_bond_range.ext()};
                for (const auto& bond : res_bonds){
                if (i == bond.idx[0] || i == bond.idx[1]) {
                    const i32 j = bond.idx[0] != i ? bond.idx[0] : bond.idx[1];
                    const auto elem = mol.atom.element[j];
                    if (elem == Element::O || elem == Element::N || elem == Element::F) {
                        donors.push_back({j, i});
                        break;
                    }
                }
            }
        }
    }
    return donors;
}

// Computes the potential acceptors given a set of atom elements.
// This essentially just a filter on atom element which extracts Oxygen and Nitrogen
DynamicArray<HydrogenBondAcceptor> compute_acceptors(const Element in_element[], i64 count) {
    DynamicArray<HydrogenBondAcceptor> acceptors;
    for (i64 i = 0; i < count; i++) {
        if (in_element[i] == Element::O || in_element[i] == Element::N || in_element[i] == Element::F) {
            acceptors.push_back(i);
        }
    }
    return acceptors;
}

// Computes hydrogen bonds given a certain set of potential donors, acceptors and atomic positions from a frame.
// The distance cutoff sets the distance from bonds to potential acceptors.
//

DynamicArray<HydrogenBond> compute_bonds(Array<const HydrogenBondDonor> donors, Array<const HydrogenBondAcceptor> acceptors, soa_vec3 in_position, float dist_cutoff, float angle_cutoff) {
    DynamicArray<HydrogenBond> bonds;

    const i32 num_acceptors = (i32)acceptors.count;
    if (!num_acceptors) return bonds;

    DynamicArray<float> acceptor_pos_x(num_acceptors);
    DynamicArray<float> acceptor_pos_y(num_acceptors);
    DynamicArray<float> acceptor_pos_z(num_acceptors);
    DynamicArray<AtomIdx> acceptor_idx(num_acceptors);

    for (i32 i = 0; i < num_acceptors; i++) {
        acceptor_pos_x[i] = in_position.x[acceptors[i]];
        acceptor_pos_y[i] = in_position.y[acceptors[i]];
        acceptor_pos_z[i] = in_position.z[acceptors[i]];
        acceptor_idx[i] = acceptors[i];
    }

    spatialhash::Frame frame = spatialhash::compute_frame(acceptor_pos_x.data(), acceptor_pos_y.data(), acceptor_pos_z.data(), num_acceptors, vec3(dist_cutoff));
    for (const auto& don : donors) {
        const vec3 donor_pos_xyz = {in_position.x[don.donor_idx], in_position.y[don.donor_idx], in_position.z[don.donor_idx]};
        const vec3 hydro_pos_xyz = {in_position.x[don.hydro_idx], in_position.y[don.hydro_idx], in_position.z[don.hydro_idx]};
        spatialhash::for_each_within(frame, hydro_pos_xyz, dist_cutoff,
                                     [&bonds, &donor_pos_xyz, &hydro_pos_xyz, &acceptor_idx, &don, angle_cutoff](i32 idx, const vec3& pos) {
                                         AtomIdx g_idx = acceptor_idx[idx];
                                         if (g_idx == don.donor_idx) return;
                                         const vec3 a = hydro_pos_xyz - donor_pos_xyz;
                                         const vec3 b = pos - hydro_pos_xyz;
                                         if (math::angle(a, b) < angle_cutoff) {
                                             bonds.push_back({g_idx, don.donor_idx, don.hydro_idx, 0.0f});
                                         }
                                     });
    }
    return bonds;
}

}  // namespace hydrogen_bond
