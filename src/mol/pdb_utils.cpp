#include "pdb_utils.h"
#include <mol/element.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>
#include <mol/trajectory_utils.h>
#include <core/string_utils.h>
#include <core/log.h>

#include <ctype.h>
#include <stdio.h>

namespace pdb {

inline int char_to_digit(uint8 c) { return c - '0'; }

inline float fast_and_unsafe_str_to_float(CString str) {
    const uint8* c = str.beg();

    float val = 0;
    float base = 1;
    float sign = 1;

    while (c != str.end()) {
		if ('0' <= *c && *c <= '9') {
			val *= 10;
			val += char_to_digit(*c);
		}
		else if (*c == '-') {
			sign = -1;
		}
		else if (*c == '.') {
			c++;
			break;
		}
        c++;
    }
    while (c != str.end() && *c != ' ') {
		if ('0' <= *c && *c <= '9') {
			base *= 0.1f;
			val += char_to_digit(*c) * base;
		}
        c++;
    }
    return sign * val;
}

inline CString extract_next_model(CString& pdb_string) {
    CString beg_mdl = find_string(pdb_string, "MODEL ");
    if (beg_mdl) {
		CString tmp_mdl = { beg_mdl.end(),  pdb_string.end() - beg_mdl.end() };
        CString end_mdl = find_string(tmp_mdl, "ENDMDL"); // @NOTE: The more characters as a search pattern, the merrier
        if (end_mdl) {
			// @NOTE: Only modify pdb_string if we found a complete model block.
			pdb_string = { end_mdl.end(), pdb_string.end() - end_mdl.end() };
            return {beg_mdl.beg(), end_mdl.end()};
        }
    }

    return {};
}

inline void extract_position(float* x, float* y, float* z, CString line) {
	// SLOW 
	// sscanf(line.substr(30).ptr, "%8f%8f%8f", &pos.x, &pos.y, &pos.z);

	// FASTER 🚴
	// pos.x = to_float(line.substr(30, 8));
	// pos.y = to_float(line.substr(38, 8));
	// pos.z = to_float(line.substr(46, 8));

	// FASTEST? 🏎️💨
	*x = fast_and_unsafe_str_to_float(line.substr(30, 8));
	*y = fast_and_unsafe_str_to_float(line.substr(38, 8));
	*z = fast_and_unsafe_str_to_float(line.substr(46, 8));
}

inline void extract_simulation_box(mat3* box, CString line) {
	vec3 dim(to_float(line.substr(6, 9)), to_float(line.substr(15, 9)), to_float(line.substr(24, 9)));
	vec3 angles(to_float(line.substr(33, 7)), to_float(line.substr(40, 7)), to_float(line.substr(47, 7)));
	// @NOTE: If we are given a zero dim, just use unit length
	if (dim == vec3(0)) dim = vec3(1);
	(*box)[0].x = dim.x;
	(*box)[1].y = dim.y;
	(*box)[2].z = dim.z;
}

bool load_dynamic_from_file(MoleculeDynamic* md, CString filename) {
    String txt = allocate_and_read_textfile(filename);
    defer { free_string(&txt); };
    if (!txt) {
        LOG_ERROR("Could not read file: '%s'.", filename);
        return false;
    }

    free_molecule_structure(&md->molecule);
    free_trajectory(&md->trajectory);
    auto res = load_dynamic_from_string(md, txt);
    return res;
}

bool load_dynamic_from_string(MoleculeDynamic* md, CString pdb_string) {
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
			extract_position(&pos.x, &pos.y, &pos.z, line);

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

bool load_molecule_from_file(MoleculeStructure* mol, CString filename) {
    StringBuffer<256> zstr = filename; // Zero terminated
	FILE* file = fopen(zstr.cstr(), "rb");
    if (!file) {
        LOG_ERROR("Could not open file: %s", zstr.cstr());
        return false;
    }

    // @NOTE: We pray to the gods above and hope that one single frame will fit in this memory
    constexpr auto mem_size = MEGABYTES(16);
    void* mem = TMP_MALLOC(mem_size);
    defer{ TMP_FREE(mem); };

	auto bytes_read = fread(mem, 1, mem_size, file);
	CString pdb_str = { (uint8*)mem, (int64)bytes_read };
	CString mdl_str = extract_next_model(pdb_str);

    return load_molecule_from_string(mol, mdl_str);
}

bool load_molecule_from_string(MoleculeStructure* mol, CString pdb_string) {
    ASSERT(mol);
    free_molecule_structure(mol);

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
    mat3 box(1);
    CString line;
    while (pdb_string && (line = extract_line(pdb_string))) {
        if (compare_n(line, "ATOM", 4) || compare_n(line, "HETATM", 6)) {
            vec3 pos;
			extract_position(&pos.z, &pos.y, &pos.z, line);
            pos_x.push_back(pos.x);
            pos_y.push_back(pos.y);
            pos_z.push_back(pos.z);

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
            num_atoms++;
        } else if (compare_n(line, "ENDMDL", 6)) {
            break;
        }
    }

    auto masses = compute_atom_masses(elements);
    auto radii = compute_atom_radii(elements);
    auto covalent_bonds = compute_covalent_bonds(residues, pos_x.data(), pos_y.data(), pos_z.data(), residue_indices.data(), elements.data(), num_atoms);
    auto backbone_segments = compute_backbone_segments(residues, labels);
    auto backbone_sequences = compute_backbone_sequences(backbone_segments, residues);
    auto backbone_angles = compute_backbone_angles(backbone_segments, backbone_sequences, pos_x.data(), pos_y.data(), pos_z.data());
    auto donors = hydrogen_bond::compute_donors(elements, residue_indices, residues, covalent_bonds);
    auto acceptors = hydrogen_bond::compute_acceptors(elements);

    if (chains.size() == 0) {
        chains = compute_chains(residues);
    }

    init_molecule_structure(mol, num_atoms, (int32)covalent_bonds.count, (int32)residues.count, (int32)chains.count, (int32)backbone_segments.count, (int32)backbone_sequences.count,
                            (int32)donors.count, (int32)acceptors.count);

    // Copy data into molecule
    memcpy(mol->atom.position.x, pos_x.data(), pos_x.size_in_bytes());
    memcpy(mol->atom.position.y, pos_y.data(), pos_y.size_in_bytes());
    memcpy(mol->atom.position.z, pos_z.data(), pos_z.size_in_bytes());
    memset(mol->atom.velocity.x, 0, num_atoms * sizeof(float));
    memset(mol->atom.velocity.y, 0, num_atoms * sizeof(float));
    memset(mol->atom.velocity.z, 0, num_atoms * sizeof(float));
    memcpy(mol->atom.radius, radii.data(), num_atoms * sizeof(float));
    memcpy(mol->atom.mass, masses.data(), num_atoms * sizeof(float));
    memcpy(mol->atom.element, elements.ptr, elements.size_in_bytes());
    memcpy(mol->atom.label, labels.ptr, labels.size_in_bytes());
    memcpy(mol->atom.res_idx, residue_indices.ptr, residue_indices.size_in_bytes());

    memcpy(mol->residues.ptr, residues.ptr, residues.size_in_bytes());
    memcpy(mol->chains.ptr, chains.ptr, chains.size_in_bytes());
    memcpy(mol->covalent_bonds.ptr, covalent_bonds.ptr, covalent_bonds.size_in_bytes());
    memcpy(mol->backbone.segments.ptr, backbone_segments.ptr, backbone_segments.size_in_bytes());
    memcpy(mol->backbone.angles.ptr, backbone_angles.ptr, backbone_angles.size_in_bytes());
    memcpy(mol->backbone.sequences.ptr, backbone_sequences.ptr, backbone_sequences.size_in_bytes());
    memcpy(mol->hydrogen_bond.donors.ptr, donors.ptr, donors.size_in_bytes());
    memcpy(mol->hydrogen_bond.acceptors.ptr, acceptors.ptr, acceptors.size_in_bytes());

    return true;
}



bool extract_info(Info* info, CString pdb_string) {
    ASSERT(info);

    int32 num_atoms = 0;
    int32 num_residues = 0;
    int32 num_chains = 0;
    int32 num_frames = 0;

    CString mdl_block = extract_next_model(pdb_string);
    if (mdl_block) {
        num_frames++;
    } else {
        mdl_block = pdb_string;
    }

	uint32 curr_res_pattern = 0;
    uint8 curr_chain_pattern = 0;
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

    info->frame_byte_ranges.push_back({(int64)mdl_block.beg(), (int64)mdl_block.end()});  // Add frame byte range of model block
    while ((mdl_block = extract_next_model(pdb_string))) {
        info->frame_byte_ranges.push_back({(int64)mdl_block.beg(), (int64)mdl_block.end()});
        num_frames++;
    }

    info->num_atoms = num_atoms;
    info->num_residues = num_residues;
    info->num_chains = num_chains;
    info->num_frames = num_frames;

    return true;
}

bool init_dynamic_from_file(MoleculeDynamic* md, CString filename) {
	ASSERT(md);
	free_molecule_structure(&md->molecule);
	free_trajectory(&md->trajectory);
    LOG_NOTE("Loading pdb dynamic from file: %.*s", (int32)filename.size_in_bytes(), filename.cstr());

    StringBuffer<256> zstr = filename;
    FILE* file = fopen(zstr.cstr(), "rb");
    if (!file) {
        LOG_ERROR("Could not open file: %s", zstr.cstr());
        return false;
    }

	constexpr auto page_size = MEGABYTES(32);
    void* mem = TMP_MALLOC(2 * page_size);
    defer { TMP_FREE(mem); };
	uint8* page[2] = {(uint8*)mem, (uint8*)mem + page_size};

	auto bytes_read = fread(page[0], 1, 2 * page_size, file);
    int64 global_offset = 0;
	CString pdb_str = {page[0], (int64)bytes_read};
    CString mdl_str = extract_next_model(pdb_str);

	if (!mdl_str) {
        LOG_ERROR("Could not locate MODEL entry in Pdb file!");
        return false;
	}
    
    // @NOTE: Search space for CRYST1 containing global simulation box parameters
    mat3 sim_box(0);
    CString box_str = {page[0], mdl_str.beg() - page[0]};
    CString line;
    while ((line = extract_line(box_str))) {
        if (compare_n(line, "CRYST1", 6)) {
            extract_simulation_box(&sim_box, line);
            break;
        }
    }
    
    load_molecule_from_string(&md->molecule, mdl_str);

	DynamicArray<int64> offsets;
    do {
        offsets.push_back(global_offset + (mdl_str.ptr - page[0]));

		// @NOTE: Have we crossed the boundry to the second page
		if (mdl_str.ptr > page[1]) {
			// Copy contents of second page to first page and read in a new page...
			memcpy(page[0], page[1], page_size);
			bytes_read = fread(page[1], 1, page_size, file);
			
			// Modify pointers accordingly
			mdl_str.ptr -= page_size;
			pdb_str.ptr -= page_size;
			pdb_str.count += bytes_read;
			global_offset += page_size;
		}
	} while ((mdl_str = extract_next_model(pdb_str)));
    
    // Time between frames
    const float dt = 1.0f;
	init_trajectory(&md->trajectory, (int32)md->molecule.atom.count, (int32)offsets.size(), dt, sim_box);

	md->trajectory.file.handle = file;
	md->trajectory.file.path = allocate_string(filename);
	md->trajectory.file.tag = PDB_FILE_TAG;
	
    md->trajectory.num_frames = 0;
	md->trajectory.frame_offsets = allocate_array<int64>(offsets.size());
	memcpy(md->trajectory.frame_offsets.data(), offsets.data(), offsets.size_in_bytes());

    rewind(file);
	pdb::read_next_trajectory_frame(&md->trajectory);

	return true;
}

bool read_next_trajectory_frame(MoleculeTrajectory* traj) {
	ASSERT(traj);
	if (traj->file.handle == 0) return false;
	if (traj->file.tag != PDB_FILE_TAG) {
		LOG_ERROR("Wrong file tag for reading trajectory frame... Expected PDB_FILE_TAG");
		return false;
	}
	auto num_frames = traj->frame_offsets.count;
	if (traj->num_frames == num_frames) return false;
    
    const int i = traj->num_frames;
    const bool last_frame = (i == num_frames - 1);

    const auto num_bytes = last_frame ?
        (traj->frame_offsets[i] - traj->frame_offsets[i-1]) :
        (traj->frame_offsets[i+1] - traj->frame_offsets[i]);
	void* mem = TMP_MALLOC(num_bytes);
	defer{ TMP_FREE(mem); };

	FSEEK((FILE*)traj->file.handle, traj->frame_offsets[i], SEEK_SET);
    const auto bytes_read = fread(mem, 1, num_bytes, (FILE*)traj->file.handle);
    
	CString mdl_str = { (uint8*)mem, (int64)bytes_read };
	TrajectoryFrame* frame = traj->frame_buffer.ptr + i;

	// Read positions
	float* x = frame->atom_position.x;
	float* y = frame->atom_position.y;
	float* z = frame->atom_position.z;
	int32 atom_idx = 0;
	CString line;
	while (mdl_str && (line = extract_line(mdl_str))) {
		if (compare_n(line, "ATOM", 4) || compare_n(line, "HETATM", 6)) {
			extract_position(x + atom_idx, y + atom_idx, z + atom_idx, line);
			atom_idx++;
		}
		else if (compare_n(line, "CRYST1", 6)) {
			extract_simulation_box(&frame->box, line);
		}
	}

	traj->num_frames++;
	return true;
}

bool close_file_handle(MoleculeTrajectory* traj) {
	ASSERT(traj);
	if (traj->file.tag != PDB_FILE_TAG) {
		LOG_ERROR("Wrong file tag for closing file handle... Expected PDB_FILE_TAG");
		return false;
	}

	if (traj->file.handle) {
		fclose((FILE*)traj->file.handle);
		traj->file.handle = nullptr;
		return true;
	}
	return false;
}

}
