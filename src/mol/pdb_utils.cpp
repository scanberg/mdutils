#include "pdb_utils.h"
#include <mol/element.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>
#include <mol/trajectory_utils.h>
#include <core/string_utils.h>
#include <core/log.h>

bool allocate_and_load_pdb_from_file(MoleculeDynamic* md, CString filename) {
    String txt = allocate_and_read_textfile(filename);
    defer { FREE(txt); };
    if (!txt) {
        LOG_ERROR("Could not read file: '%s'.", filename);
        return false;
    }

    free_molecule_structure(&md->molecule);
    free_trajectory(&md->trajectory);
    auto res = allocate_and_parse_pdb_from_string(md, txt);
    return res;
}

inline int char_to_digit(uint8 c) { return c - '0'; }

#include <ctype.h>
inline float fast_and_unsafe_str_to_float(CString str) {
    const uint8* c = str.beg();
    while (c != str.end() && *c == ' ') ++c;
    const uint8* end = c;
    while (end != str.end() && *end != ' ') ++end;
    if (c == end) return 0;

    float val = 0;
    float base = 1;
    float sign = 1;
    if (*c == '-') {
        sign = -1;
        c++;
    }
    while (c != end && *c != '.') {
        val *= 10;
        val += char_to_digit(*c);
        c++;
    }
    if (c == end) return sign * val;
    c++;
    while (c != end) {
        base *= 0.1f;
        val += char_to_digit(*c) * base;
        c++;
    }

    return sign * val;
}

bool allocate_and_parse_pdb_from_string(MoleculeDynamic* md, CString pdb_string) {
    free_molecule_structure(&md->molecule);
    free_trajectory(&md->trajectory);

    DynamicArray<vec3> positions;
    DynamicArray<Label> labels;
    DynamicArray<Element> elements;
    DynamicArray<ResIdx> residue_indices;
    DynamicArray<float> occupancies;
    DynamicArray<float> temp_factors;
    DynamicArray<Residue> residues;
    DynamicArray<Chain> chains;

    positions.reserve(1024);
    labels.reserve(1024);
    elements.reserve(1024);
    residue_indices.reserve(1024);
    occupancies.reserve(1024);
    temp_factors.reserve(1024);
    residues.reserve(128);
    chains.reserve(32);

    int current_res_id = -1;
    char current_chain_id = -1;
    int num_atoms = 0;
    int num_frames = 0;
    mat3 box(0);
    CString line;
    while ((line = extract_line(pdb_string))) {
        if (compare_n(line, "ATOM", 4) || compare_n(line, "HETATM", 6)) {
            vec3 pos;

            // SLOW AS 💩
            // sscanf(line.substr(30).data, "%8f%8f%8f", &pos.x, &pos.y, &pos.z);

            // FASTER 🚴
            // pos.x = to_float(line.substr(30, 8));
            // pos.y = to_float(line.substr(38, 8));
            // pos.z = to_float(line.substr(46, 8));

            // FASTEST? 🏎️💨
            pos.x = fast_and_unsafe_str_to_float(line.substr(30, 8));
            pos.y = fast_and_unsafe_str_to_float(line.substr(38, 8));
            pos.z = fast_and_unsafe_str_to_float(line.substr(46, 8));

            positions.push_back(pos);
            if (num_frames > 0) continue;

            labels.push_back(trim(line.substr(12, 4)));
            if (line.count > 60) {
                occupancies.push_back(to_float(line.substr(54, 6)));
            }
            if (line.count > 66) {
                temp_factors.push_back(to_float(line.substr(60, 6)));
            }

            // Try to determine element from optional element column first, then from label
            Element elem = Element::Unknown;
            if (line.count >= 78) {
                elem = element::get_from_string(line.substr(76, 2));
            }
            if (elem == Element::Unknown) {
                elem = element::get_from_string(labels.back());
            }
            elements.push_back(elem);

            auto res_id = to_int(line.substr(22, 4));
            char chain_id = line[21];

            // New Chain
            if (current_chain_id != chain_id && chain_id != ' ') {
                current_chain_id = chain_id;
                Chain chain;
                chain.res_idx.beg = (ResIdx)residues.size();
                chain.res_idx.end = (ResIdx)residues.size();
                chain.id = chain_id;
                chains.push_back(chain);
            }

            // New Residue
            if (res_id != current_res_id) {
                current_res_id = res_id;
                Residue res{};
                res.name = trim(line.substr(17, 3));
                res.id = res_id;
                res.chain_idx = (ChainIdx)(chains.size() - 1);
                res.atom_idx = {num_atoms, num_atoms};
                residues.push_back(res);
                if (chains.size() > 0) {
                    chains.back().res_idx.end++;
                }
            }
            residues.back().atom_idx.end++;

            residue_indices.push_back((ResIdx)(residues.size() - 1));

            // Add Atom
            num_atoms++;
        } else if (compare_n(line, "CRYST1", 6)) {
            vec3 dim(to_float(line.substr(6, 9)), to_float(line.substr(15, 9)), to_float(line.substr(24, 9)));
            vec3 angles(to_float(line.substr(33, 7)), to_float(line.substr(40, 7)), to_float(line.substr(47, 7)));
            // @NOTE: If we are given a zero dim, just use unit length
            if (dim == vec3(0)) dim = vec3(1);
            box[0].x = dim.x;
            box[1].y = dim.y;
            box[2].z = dim.z;
        } else if (compare_n(line, "ENDMDL", 6)) {
            num_frames++;
            // @TODO: Handle the case where the models are different and not consecutive frames of an animation.
        }
    }

    if (!md->molecule) {
        auto mol_pos = positions.sub_array(0, num_atoms);
        auto covalent_bonds = compute_covalent_bonds(residues, residue_indices, mol_pos, elements);
        auto backbone_segments = compute_backbone_segments(residues, labels);
        auto backbone_sequences = compute_backbone_sequences(backbone_segments, residues);
        auto backbone_angles = compute_backbone_angles(mol_pos, backbone_segments, backbone_sequences);
        auto donors = hydrogen_bond::compute_donors(elements, residue_indices, residues, covalent_bonds);
        auto acceptors = hydrogen_bond::compute_acceptors(elements);
        if (chains.size() == 0) {
            chains = compute_chains(residues);
        }

        init_molecule_structure(&md->molecule, num_atoms, (int32)covalent_bonds.count, (int32)residues.count, (int32)chains.count, (int32)backbone_segments.count, (int32)backbone_sequences.count,
                                (int32)donors.count, (int32)acceptors.count);

        // Copy data into molecule
        memcpy(md->molecule.atom.positions, mol_pos.data, mol_pos.size_in_bytes());
        memcpy(md->molecule.atom.elements, elements.data, elements.size_in_bytes());
        memcpy(md->molecule.atom.labels, labels.data, labels.size_in_bytes());
        memcpy(md->molecule.atom.residue_indices, residue_indices.data, residue_indices.size_in_bytes());

        memcpy(md->molecule.residues.data, residues.data, residues.size_in_bytes());
        memcpy(md->molecule.chains.data, chains.data, chains.size_in_bytes());
        memcpy(md->molecule.covalent_bonds.data, covalent_bonds.data, covalent_bonds.size_in_bytes());
        memcpy(md->molecule.backbone.segments.data, backbone_segments.data, backbone_segments.size_in_bytes());
        memcpy(md->molecule.backbone.angles.data, backbone_angles.data, backbone_angles.size_in_bytes());
        memcpy(md->molecule.backbone.sequences.data, backbone_sequences.data, backbone_sequences.size_in_bytes());
        memcpy(md->molecule.hydrogen_bond.donors.data, donors.data, donors.size_in_bytes());
        memcpy(md->molecule.hydrogen_bond.acceptors.data, acceptors.data, acceptors.size_in_bytes());
    }

    if (num_frames > 0) {
        init_trajectory(&md->trajectory, num_atoms, num_frames);
        // COPY POSITION DATA

        ASSERT(positions.count > 0);
        memcpy(md->trajectory.position_data.data, positions.data, sizeof(vec3) * positions.count);

        for (int i = 0; i < md->trajectory.num_frames; i++) {
            int index = i;
            float time = 0;
            Array<vec3> atom_positions{md->trajectory.position_data.beg() + i * num_atoms, num_atoms};
            md->trajectory.frame_buffer[i] = {index, time, box, atom_positions};
        }
    }

    return true;
}
