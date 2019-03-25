#include "pdb_utils.h"
#include <mol/element.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>
#include <mol/trajectory_utils.h>
#include <core/string_utils.h>
#include <core/log.h>

bool allocate_and_load_pdb_from_file(MoleculeDynamic* md, CString filename) {
    String txt = allocate_and_read_textfile(filename);
    defer { free_string(&txt); };
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

    DynamicArray<float> pos_x;
    DynamicArray<float> pos_y;
    DynamicArray<float> pos_z;
    DynamicArray<Label> labels;
    DynamicArray<Element> elements;
    DynamicArray<ResIdx> residue_indices;
    DynamicArray<float> occupancies;
    DynamicArray<float> temp_factors;
    DynamicArray<Residue> residues;
    DynamicArray<Chain> chains;

    pos_x.reserve(1024);
    pos_y.reserve(1024);
    pos_z.reserve(1024);
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
    while (pdb_string && (line = extract_line(pdb_string))) {
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

            pos_x.push_back(pos.x);
            pos_y.push_back(pos.y);
            pos_z.push_back(pos.z);

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
                chain.res_range = {(ResIdx)residues.size(), (ResIdx)residues.size()};
                chain.atom_range = {num_atoms, num_atoms};
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
                res.atom_range = {num_atoms, num_atoms};
                residues.push_back(res);
                if (chains.size() > 0) {
                    chains.back().res_range.end++;
                }
            }
            if (residues.size() > 0) residues.back().atom_range.end++;
            if (chains.size() > 0) chains.back().atom_range.end++;

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
        auto mol_pos_x = pos_x.subarray(0, num_atoms);
        auto mol_pos_y = pos_y.subarray(0, num_atoms);
        auto mol_pos_z = pos_z.subarray(0, num_atoms);

        auto masses = compute_atom_masses(elements);
        auto radii = compute_atom_radii(elements);
        auto covalent_bonds = compute_covalent_bonds(residues, mol_pos_x.data(), mol_pos_y.data(), mol_pos_z.data(), residue_indices.data(), elements.data(), num_atoms);
        auto backbone_segments = compute_backbone_segments(residues, labels);
        auto backbone_sequences = compute_backbone_sequences(backbone_segments, residues);
        auto backbone_angles = compute_backbone_angles(backbone_segments, backbone_sequences, mol_pos_x.data(), mol_pos_y.data(), mol_pos_z.data());
        auto donors = hydrogen_bond::compute_donors(elements, residue_indices, residues, covalent_bonds);
        auto acceptors = hydrogen_bond::compute_acceptors(elements);

        if (chains.size() == 0) {
            chains = compute_chains(residues);
        }

        init_molecule_structure(&md->molecule, num_atoms, (int32)covalent_bonds.count, (int32)residues.count, (int32)chains.count, (int32)backbone_segments.count, (int32)backbone_sequences.count,
                                (int32)donors.count, (int32)acceptors.count);

        // Copy data into molecule
        memcpy(md->molecule.atom.position.x, mol_pos_x.data(), mol_pos_x.size_in_bytes());
        memcpy(md->molecule.atom.position.y, mol_pos_y.data(), mol_pos_y.size_in_bytes());
        memcpy(md->molecule.atom.position.z, mol_pos_z.data(), mol_pos_z.size_in_bytes());
        memset(md->molecule.atom.velocity.x, 0, num_atoms * sizeof(float));
        memset(md->molecule.atom.velocity.y, 0, num_atoms * sizeof(float));
        memset(md->molecule.atom.velocity.z, 0, num_atoms * sizeof(float));
        memcpy(md->molecule.atom.radius, radii.data(), num_atoms * sizeof(float));
        memcpy(md->molecule.atom.mass, masses.data(), num_atoms * sizeof(float));
        memcpy(md->molecule.atom.element, elements.ptr, elements.size_in_bytes());
        memcpy(md->molecule.atom.label, labels.ptr, labels.size_in_bytes());
        memcpy(md->molecule.atom.res_idx, residue_indices.ptr, residue_indices.size_in_bytes());

        memcpy(md->molecule.residues.ptr, residues.ptr, residues.size_in_bytes());
        memcpy(md->molecule.chains.ptr, chains.ptr, chains.size_in_bytes());
        memcpy(md->molecule.covalent_bonds.ptr, covalent_bonds.ptr, covalent_bonds.size_in_bytes());
        memcpy(md->molecule.backbone.segments.ptr, backbone_segments.ptr, backbone_segments.size_in_bytes());
        memcpy(md->molecule.backbone.angles.ptr, backbone_angles.ptr, backbone_angles.size_in_bytes());
        memcpy(md->molecule.backbone.sequences.ptr, backbone_sequences.ptr, backbone_sequences.size_in_bytes());
        memcpy(md->molecule.hydrogen_bond.donors.ptr, donors.ptr, donors.size_in_bytes());
        memcpy(md->molecule.hydrogen_bond.acceptors.ptr, acceptors.ptr, acceptors.size_in_bytes());
    }

    if (num_frames > 0) {
        const float time_between_frames = 1.0f;
        init_trajectory(&md->trajectory, num_atoms, num_frames, time_between_frames, box);

        memcpy(md->trajectory.position_data.x, pos_x.data(), pos_x.size_in_bytes());
        memcpy(md->trajectory.position_data.y, pos_y.data(), pos_y.size_in_bytes());
        memcpy(md->trajectory.position_data.z, pos_z.data(), pos_z.size_in_bytes());
    }

    return true;
}

inline CString extract_next_model(CString& pdb_string) {
    CString beg_mdl = find_string(pdb_string, "\nMODEL ");
    if (beg_mdl) {
        pdb_string.count = pdb_string.end() - beg_mdl.end();
        pdb_string.ptr = beg_mdl.end();

        CString end_mdl = find_string(pdb_string, "\nENDMDL");
        if (end_mdl) {
            pdb_string.count = pdb_string.end() - end_mdl.end();
            pdb_string.ptr = end_mdl.end();
            return {beg_mdl.beg(), end_mdl.end()};
        }
    }

    return {};
}

bool extract_pdb_info(PdbInfo* info, CString pdb_string) {
	ASSERT(info);
	
	int32 num_atoms = 0;
    int32 num_residues = 0;
    int32 num_chains = 0;
    int32 num_frames = 0;

    uint32 curr_res_pattern = 0;
    uint8 curr_chain_pattern = 0;

    CString mdl_block = extract_next_model(pdb_string);
    if (mdl_block) {
        num_frames++;
    } else {
        mdl_block = pdb_string;
    }

    CString line;
    while (mdl_block && (line = extract_line(mdl_block))) {
        if (compare_n(line, "ATOM", 4) || compare_n(line, "HETATM", 6)) {
            const uint32 res_pattern = *(uint32*)(&line[22]);
            const uint8 chain_pattern = line[21];

            num_atoms++;
            if (res_pattern != curr_res_pattern) {
                num_residues++;
                curr_res_pattern = res_pattern;
            }
            if (chain_pattern != curr_chain_pattern) {
                num_chains++;
                curr_chain_pattern = chain_pattern;
            }
        }
    }

	info->frame_byte_ranges.push_back({ (int64)mdl_block.beg(), (int64)mdl_block.end() }); // Add frame byte range of model block
    while ((mdl_block = extract_next_model(pdb_string))) {
		info->frame_byte_ranges.push_back({ (int64)mdl_block.beg(), (int64)mdl_block.end() });
        num_frames++;
    }

    info->num_atoms = num_atoms;
    info->num_residues = num_residues;
    info->num_chains = num_chains;
    info->num_frames = num_frames;

    return true;
}
