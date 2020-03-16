#include "molecule_structure.h"
#include <core/spatial_hash.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>
#include <mol/element_utils.h>

// 64-byte alignment for 512-bit vectorization (AVX512 and beyond!)
#define ALIGNMENT 64

// Computes covalent bonds between a set of atoms with given positions and elements.
// The approach is inspired by the technique used in NGL (https://github.com/arose/ngl)
inline bool covelent_bond_heuristic(const vec3& p0, Element e0, const vec3& p1, Element e1) {
    const float d = element::covalent_radius(e0) + element::covalent_radius(e1);
    const float d_max = d + 0.3f;
    const float d_min = d - 0.5f;
    const float d2 = math::distance2(p0, p1);
    return (d_min * d_min) < d2 && d2 < (d_max * d_max);
}

static DynamicArray<Bond> compute_covalent_bonds(MoleculeStructure& mol) {
    constexpr float max_covelent_bond_length = 4.0f;
    DynamicArray<Bond> bonds;
    spatialhash::Frame frame;

    for (i64 ri = 0; ri < mol.residue.count - 1; ri++) {
        const AtomRange ri_range = mol.residue.atom_range[ri];

        if (ri > 0) {
            // Include potential shared bonds from previous residue
            mol.residue.bond.complete[ri].beg = mol.residue.bond.intra[ri - 1].end;
            mol.residue.bond.intra[ri].beg = mol.residue.bond.intra[ri].end = mol.residue.bond.complete[ri].end =
                mol.residue.bond.complete[ri - 1].end;
        } else {
            mol.residue.bond.complete[ri] = {(BondIdx)bonds.size(), (BondIdx)bonds.size()};
            mol.residue.bond.intra[ri] = {(BondIdx)bonds.size(), (BondIdx)bonds.size()};
        }

        // Compute internal bonds
        {
            spatialhash::compute_frame(&frame, mol.atom.position.x + mol.residue.atom_range[ri].beg,
                                       mol.atom.position.y + mol.residue.atom_range[ri].beg, mol.atom.position.z + mol.residue.atom_range[ri].beg,
                                       ri_range.ext(), vec3(max_covelent_bond_length));
            for (AtomIdx i = ri_range.beg; i < ri_range.end; i++) {
                const vec3& atom_i_pos = get_atom_position(mol, i);
                spatialhash::for_each_within(
                    frame, atom_i_pos, max_covelent_bond_length, [&bonds, &mol, atom_i_pos, i, offset = ri_range.beg](int j, const vec3& atom_j_pos) {
                        (void)atom_j_pos;
                        j += offset;
                        if (i < j && covelent_bond_heuristic(atom_i_pos, mol.atom.element[i], atom_j_pos, mol.atom.element[j])) {
                            bonds.push_back({{i, j}});
                        }
                    });
            }
            mol.residue.bond.intra[ri].end = (BondIdx)bonds.size();
        }

        // Compute bonds to next
        {
            const i64 rj = ri + 1;
            const AtomRange rj_range = mol.residue.atom_range[rj];

            spatialhash::compute_frame(&frame, mol.atom.position.x + rj_range.beg, mol.atom.position.y + rj_range.beg,
                                       mol.atom.position.z + rj_range.beg, rj_range.ext(), vec3(max_covelent_bond_length));
            for (AtomIdx i = ri_range.beg; i < ri_range.end; i++) {
                const vec3& atom_i_pos = get_atom_position(mol, i);
                spatialhash::for_each_within(
                    frame, atom_i_pos, max_covelent_bond_length, [&bonds, &mol, atom_i_pos, i, offset = rj_range.beg](int j, const vec3& atom_j_pos) {
                        (void)atom_j_pos;
                        j += offset;
                        if (i < j && covelent_bond_heuristic(atom_i_pos, mol.atom.element[i], atom_j_pos, mol.atom.element[j])) {
                            bonds.push_back({{i, j}});
                        }
                    });
            }
            mol.residue.bond.complete[ri].end = (BondIdx)bonds.size();
        }
    }

    return bonds;
}

template <i64 N>
static bool match(const Label& lbl, const char (&cstr)[N]) {
    for (i64 i = 0; i < N; i++) {
        if (to_lower(lbl[i]) != to_lower(cstr[i])) return false;
    }
    return true;
}

struct BackboneData {
    DynamicArray<BackboneSegment> segments;
    DynamicArray<BackboneSequence> sequences;
};

static BackboneData compute_backbone(const MoleculeStructure& mol) {
    constexpr i32 min_atom_count = 4;  // Must contain at least 4 atoms to be considered as an amino acid.

    BackboneData data{};
    BackboneSequence seq{};

    for (i64 ci = 0; ci < mol.chain.count; ci++) {
        for (i64 ri = mol.chain.residue_range[ci].beg; ri < mol.chain.residue_range[ci].end; ri++) {
            BackboneSegment seg = {-1, -1, -1, -1, ri};
            if (mol.residue.atom_range[ri].ext() >= min_atom_count) {
                // find atoms
                for (i32 i = mol.residue.atom_range[ri].beg; i < mol.residue.atom_range[ri].end; i++) {
                    const auto& lbl = mol.atom.label[i];
                    if (seg.ca_idx == -1 && match(lbl, "CA")) seg.ca_idx = i;
                    if (seg.n_idx  == -1 && match(lbl, "N"))  seg.n_idx  = i;
                    if (seg.c_idx  == -1 && match(lbl, "C"))  seg.c_idx  = i;
                    if (seg.o_idx  == -1 && match(lbl, "O"))  seg.o_idx  = i;
                }

                // Could not match "O"
                if (seg.o_idx == -1) {
                    // Pick first atom containing O after C atom
                    for (i32 i = seg.c_idx; i < mol.residue.atom_range[ri].end; i++) {
                        const auto& lbl = mol.atom.label[i];
                        if (lbl[0] == 'o' || lbl[0] == 'O') {
                            seg.o_idx = i;
                            break;
                        }
                    }
                }
            }

            if ((seg.ca_idx != -1) && (seg.n_idx != -1) && (seg.c_idx != -1) && (seg.o_idx != -1) && valid_segment(seg)) {
                data.segments.push_back(seg);
                seq.end++;
            } else {
                if (seq.ext() > 1) {
                    data.sequences.push_back(seq);
                }
                seq = {(SegIdx)data.segments.size(), (SegIdx)data.segments.size()};
            }
        }
        if (seq.ext() > 1) {
            data.sequences.push_back(seq);
            seq = {(SegIdx)data.segments.size(), (SegIdx)data.segments.size()};
        }
    }

    return data;
}

bool init_molecule_structure(MoleculeStructure* mol, const MoleculeStructureDescriptor& desc) {
    ASSERT(mol);
    free_molecule_structure(mol);

    if (desc.num_atoms > 0) {
        // Allocate Aligned data (@NOTE: Is perhaps not necessary as trajectory data is not aligned anyways...)
        const i64 aligned_size = (desc.num_atoms * sizeof(float) + ALIGNMENT) * 5;
        void* aligned_mem = ALIGNED_MALLOC(aligned_size, ALIGNMENT);
        if (!aligned_mem) return false;
        memset(aligned_mem, 0, aligned_size);

        mol->atom.count = desc.num_atoms;
        mol->atom.position.x = (float*)aligned_mem;
        mol->atom.position.y = (float*)get_next_aligned_adress(mol->atom.position.x + desc.num_atoms, ALIGNMENT);
        mol->atom.position.z = (float*)get_next_aligned_adress(mol->atom.position.y + desc.num_atoms, ALIGNMENT);
        mol->atom.radius = (float*)get_next_aligned_adress(mol->atom.position.z + desc.num_atoms, ALIGNMENT);
        mol->atom.mass = (float*)get_next_aligned_adress(mol->atom.radius + desc.num_atoms, ALIGNMENT);

        const i64 other_size = desc.num_atoms * (sizeof(Element) + sizeof(Label) + sizeof(ResIdx) + sizeof(ChainIdx));
        void* other_mem = MALLOC(other_size);
        if (!other_mem) return false;
        memset(other_mem, 0, other_size);

        mol->atom.element = (Element*)other_mem;
        mol->atom.label = (Label*)(mol->atom.element + mol->atom.count);
        mol->atom.res_idx = (ResIdx*)(mol->atom.label + mol->atom.count);
        mol->atom.chain_idx = (ChainIdx*)(mol->atom.res_idx + mol->atom.count);

        if (desc.atoms) {
            for (i32 i = 0; i < desc.num_atoms; i++) {
                mol->atom.position.x[i] = desc.atoms[i].x;
                mol->atom.position.y[i] = desc.atoms[i].y;
                mol->atom.position.z[i] = desc.atoms[i].z;

                mol->atom.element[i] = desc.atoms[i].element;
                if (mol->atom.element[i] == Element::Unknown) {
                    mol->atom.element[i] = get_element_from_string(desc.atoms[i].label);
                }

                mol->atom.radius[i] = element::vdw_radius(mol->atom.element[i]);
                mol->atom.mass[i] = element::atomic_mass(mol->atom.element[i]);
                mol->atom.label[i] = desc.atoms[i].label;
                mol->atom.res_idx[i] = desc.atoms[i].residue_index;
                mol->atom.chain_idx[i] = INVALID_CHAIN_IDX;
            }
        }
    }

    if (desc.num_residues > 0) {
        const i64 mem_size = desc.num_residues * (sizeof(ResIdx) + sizeof(Label) + sizeof(AtomRange) + 2 * sizeof(BondRange));
        void* mem = MALLOC(mem_size);
        if (!mem) return false;
        memset(mem, 0, mem_size);

        mol->residue.count = desc.num_residues;
        mol->residue.id = (ResIdx*)mem;
        mol->residue.name = (Label*)(mol->residue.id + desc.num_residues);
        mol->residue.atom_range = (AtomRange*)(mol->residue.name + desc.num_residues);
        mol->residue.bond.complete = (BondRange*)(mol->residue.atom_range + desc.num_residues);
        mol->residue.bond.intra = (BondRange*)(mol->residue.bond.complete + desc.num_residues);
        if (desc.residues) {
            for (i32 i = 0; i < desc.num_residues; i++) {
                mol->residue.id[i] = desc.residues[i].id;
                mol->residue.name[i] = desc.residues[i].name;
                mol->residue.atom_range[i] = desc.residues[i].atom_range;
                // bond
            }
        }
    }

    // We ignore bonds for now and compute our own
    {
        DynamicArray<Bond> bonds = compute_covalent_bonds(*mol);
        mol->covalent_bond.bond = (Bond*)MALLOC(bonds.size_in_bytes());
        memcpy(mol->covalent_bond.bond, bonds.data(), bonds.size_in_bytes());
        mol->covalent_bond.count = bonds.size();
    }

    if (desc.num_chains > 0) {
        const i64 mem_size = desc.num_chains * (sizeof(Label) + sizeof(AtomRange) + sizeof(ResRange));
        void* mem = MALLOC(mem_size);
        memset(mem, 0, mem_size);

        mol->chain.count = desc.num_chains;
        mol->chain.id = (Label*)mem;
        mol->chain.atom_range = (AtomRange*)(mol->chain.id + mol->chain.count);
        mol->chain.residue_range = (ResRange*)(mol->chain.atom_range + mol->chain.count);
        if (desc.chains) {
            for (i32 i = 0; i < desc.num_chains; i++) {
                mol->chain.id[i] = desc.chains[i].id;
                mol->chain.residue_range[i] = desc.chains[i].residue_range;
            }
        }
    } else if (desc.num_residues > 0 && desc.residues) {
        // Generate artificial chains for every connected sequence of residues, ignore single residues e.g ext() == 1
        DynamicArray<ResRange> seq;
        ResRange range{0, 1};
        for (i64 i = 0; i < mol->residue.count - 1; i++) {
            if (ranges_overlap(mol->residue.bond.complete[i], mol->residue.bond.complete[i + 1])) {
                range.end++;
            } else {
                if (range.ext() > 1) {
                    range.end++;
                    seq.push_back(range);
                }
                range = {(ResIdx)i + 1, (ResIdx)i + 2};
            }
        }

        const i64 mem_size = seq.size() * (sizeof(Label) + sizeof(AtomRange) + sizeof(ResRange));
        void* mem = MALLOC(mem_size);

        mol->chain.count = seq.size();
        mol->chain.id = (Label*)mem;
        mol->chain.atom_range = (AtomRange*)(mol->chain.id + mol->chain.count);
        mol->chain.residue_range = (ResRange*)(mol->chain.atom_range + mol->chain.count);

        memset(mol->chain.id, 0, seq.size() * sizeof(Label));
        for (i64 i = 0; i < seq.size(); i++) {
            if (seq.size() < ('Z' - 'A')) {
                mol->chain.id[i] = 'A' + i;
            } else {
                snprintf(mol->chain.id[i].buffer, Label::MaxSize, "C%i", (int)i);
            }
            mol->chain.residue_range[i] = seq[i];
        }
    }

    // Fix missing info for chains
    for (i64 i = 0; i < mol->chain.count; i++) {
        // Compute atom ranges for chains
        auto rr = mol->chain.residue_range[i];
        if (rr.ext() > 0)
            mol->chain.atom_range[i] = {mol->residue.atom_range[rr.beg].beg, mol->residue.atom_range[rr.end - 1].end};

        // Set references atoms to chains
        for (i64 j = mol->chain.atom_range[i].beg; j < mol->chain.atom_range[i].end; j++) {
            mol->atom.chain_idx[j] = i;
        }
    }

    {
        BackboneData bb_data = compute_backbone(*mol);

        const i64 mem_size = bb_data.segments.size() * (sizeof(BackboneSegment) + sizeof(BackboneAngle) + sizeof(SecondaryStructure));
        void* mem = MALLOC(mem_size);
        mol->backbone.segment.count = bb_data.segments.size();
        mol->backbone.segment.segment = (BackboneSegment*)mem;
        mol->backbone.segment.angle = (BackboneAngle*)(mol->backbone.segment.segment + mol->backbone.segment.count);
        mol->backbone.segment.secondary_structure = (SecondaryStructure*)(mol->backbone.segment.angle + mol->backbone.segment.count);

        for (i64 i = 0; i < bb_data.segments.size(); i++) {
            mol->backbone.segment.segment[i] = bb_data.segments[i];
            mol->backbone.segment.angle[i] = {};
            mol->backbone.segment.secondary_structure[i] = {};
        }

        compute_backbone_angles(mol->backbone.segment.angle, mol->atom.position, mol->backbone.segment.segment, mol->backbone.segment.count);

        mol->backbone.sequence.segment_range = (SegRange*)MALLOC(bb_data.sequences.size_in_bytes());
        mol->backbone.sequence.count = bb_data.sequences.size();
        memcpy(mol->backbone.sequence.segment_range, bb_data.sequences.data(), bb_data.sequences.size_in_bytes());
    }

    if (desc.num_secondary_structures > 0) {
        // TODO: Implement configurations of secondary structure onto backbone
    }

    {
        auto acc = hydrogen_bond::compute_acceptors(mol->atom.element, mol->atom.count);
        auto don = hydrogen_bond::compute_donors(*mol);
        const i64 mem_size = acc.size_in_bytes() + don.size_in_bytes();
        void* mem = MALLOC(mem_size);
        mol->hydrogen_bond.donor.count = don.size();
        mol->hydrogen_bond.donor.data = (HydrogenBondDonor*)mem;
        memcpy(mol->hydrogen_bond.donor.data, don.data(), don.size_in_bytes());
        mol->hydrogen_bond.acceptor.count = acc.size();
        mol->hydrogen_bond.acceptor.data = (HydrogenBondAcceptor*)(mol->hydrogen_bond.donor.data + don.size());
        memcpy(mol->hydrogen_bond.acceptor.data, acc.data(), acc.size_in_bytes());
    }

    return true;
}

void free_molecule_structure(MoleculeStructure* mol) {
    ASSERT(mol);
    if (mol->atom.position.x) ALIGNED_FREE(mol->atom.position.x);
    if (mol->atom.element) FREE(mol->atom.element);
    if (mol->covalent_bond.bond) FREE(mol->covalent_bond.bond);
    if (mol->residue.id) FREE(mol->residue.id);
    if (mol->backbone.segment.segment) FREE(mol->backbone.segment.segment);
    if (mol->backbone.sequence.segment_range) FREE(mol->backbone.sequence.segment_range);

    *mol = {};
}
