#include "pdb_utils.h"
#include <mol/element.h>
#include <mol/element_utils.h>
#include <mol/aminoacid.h>
#include <mol/aminoacid_utils.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>
#include <mol/trajectory_utils.h>
#include <core/string_utils.h>
#include <core/log.h>
#include <core/file.h>

#include <ctype.h>
#include <stdio.h>

namespace pdb {

inline CStringView extract_next_model(CStringView& pdb_string) {
    CStringView beg_mdl = find_string(pdb_string, "\nMODEL ");
    if (beg_mdl) {
        CStringView tmp_mdl = {beg_mdl.end(), pdb_string.end() - beg_mdl.end()};
        CStringView end_mdl = find_string(tmp_mdl, "\nENDMDL");  // @NOTE: The more characters as a search pattern, the merrier
        if (end_mdl) {
            // @NOTE: Only modify pdb_string if we found a complete model block.
            pdb_string = {end_mdl.end(), pdb_string.end() - end_mdl.end()};
            return {beg_mdl.beg(), end_mdl.end()};
        }
    }

    return {};
}

inline void extract_position(float* x, float* y, float* z, CStringView line) {
    // SLOW 🐢
    // sscanf(line.substr(30).ptr, "%8f%8f%8f", &pos.x, &pos.y, &pos.z);

    // FASTER 🚴
    // pos.x = to_float(line.substr(30, 8));
    // pos.y = to_float(line.substr(38, 8));
    // pos.z = to_float(line.substr(46, 8));

    // FASTEST 🏎️💨
    *x = str_to_float(line.substr(30, 8));
    *y = str_to_float(line.substr(38, 8));
    *z = str_to_float(line.substr(46, 8));
}

inline void extract_simulation_box(mat3* box, CStringView line) {
    vec3 dim(to_float(line.substr(6, 9)), to_float(line.substr(15, 9)), to_float(line.substr(24, 9)));
    vec3 angles(to_float(line.substr(33, 7)), to_float(line.substr(40, 7)), to_float(line.substr(47, 7)));
    // @NOTE: If we are given a zero dim, just use unit length
    if (dim == vec3(0)) dim = vec3(1);
    (*box)[0].x = dim.x;
    (*box)[1].y = dim.y;
    (*box)[2].z = dim.z;
}

inline void extract_element(Element* element, CStringView line) {
    Element elem = Element::Unknown;
    if (line.size() >= 78) {
        // @NOTE: Try optional atom element field, ignore case
        elem = get_element_from_string(line.substr(76, 2), true);
    }

    if (elem == Element::Unknown) {
        // @NOTE: Try to deduce from atom id
        const CStringView atom_id = line.substr(12, 4);
        const CStringView res_name = line.substr(17, 3);
        if (compare_n(atom_id, "CA", 2) && get_amino_acid_from_string(res_name) == AminoAcid::Unknown) {
            // @NOTE: Ambigous case where CA is probably calcium if not part of an amino acid
            elem = Element::Ca;
        } else {
            elem = get_element_from_string(atom_id);
        }
    }
    *element = elem;
}

inline void extract_helix_residue_indices(int* init_idx, int* term_idx, CStringView line) {
    ASSERT(init_idx);
    ASSERT(term_idx);
    ASSERT(line.length() >= 37);
    *init_idx = to_int32(line.substr(21, 4));
    *term_idx = to_int32(line.substr(33, 4));
}

inline void extract_sheet_residue_indices(int* init_idx, int* term_idx, CStringView line) {
    ASSERT(init_idx);
    ASSERT(term_idx);
    ASSERT(line.length() >= 37);
    *init_idx = to_int32(line.substr(22, 4));
    *term_idx = to_int32(line.substr(33, 4));
}

inline bool extract_trajectory_frame_data(TrajectoryFrame* frame, i32 num_atoms, CStringView mdl_str) {
    ASSERT(frame);
    float* x = frame->atom_position.x;
    float* y = frame->atom_position.y;
    float* z = frame->atom_position.z;
    i32 atom_idx = 0;
    CStringView line;
    while (mdl_str && (line = extract_line(mdl_str)) && atom_idx < num_atoms) {
        if (compare_n(line, "ATOM", 4) || compare_n(line, "HETATM", 6)) {
            extract_position(x + atom_idx, y + atom_idx, z + atom_idx, line);
            atom_idx++;
        } else if (compare_n(line, "CRYST1", 6)) {
            extract_simulation_box(&frame->box, line);
        }
    }
    return atom_idx == num_atoms;
}

bool load_molecule_from_file(MoleculeStructure* mol, CStringView filename) {
    FILE* file = fopen(filename, "rb");
    defer { fclose(file); };
    if (!file) {
        LOG_ERROR("Could not open file: %.*s", filename.length(), filename.cstr());
        return false;
    }

    // @NOTE: We pray to the gods above and hope that one single frame will fit in this memory
    constexpr auto mem_size = MEGABYTES(32);
    void* mem = TMP_MALLOC(mem_size);
    defer { TMP_FREE(mem); };

    auto bytes_read = fread(mem, 1, mem_size, file);
    CStringView pdb_str = {(const char*)mem, (i64)bytes_read};

    return load_molecule_from_string(mol, pdb_str);
}

bool load_molecule_from_string(MoleculeStructure* mol, CStringView pdb_string) {
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
    DynamicArray<Range<ResIdx>> helices;
    DynamicArray<Range<ResIdx>> sheets;

    constexpr auto atom_reserve_size = 4096;
    constexpr auto residue_reserve_size = 128;
    constexpr auto chain_reserve_size = 16;

    pos_x.reserve(atom_reserve_size);
    pos_y.reserve(atom_reserve_size);
    pos_z.reserve(atom_reserve_size);
    labels.reserve(atom_reserve_size);
    elements.reserve(atom_reserve_size);
    residue_indices.reserve(atom_reserve_size);
    occupancies.reserve(atom_reserve_size);
    temp_factors.reserve(atom_reserve_size);
    residues.reserve(residue_reserve_size);
    chains.reserve(chain_reserve_size);

    int current_res_id = -1;
    char current_chain_id = -1;
    int num_atoms = 0;
    mat3 box(1);
    CStringView line;
    while (pdb_string && (line = extract_line(pdb_string))) {
        if (compare_n(line, "ATOM", 4) || compare_n(line, "HETATM", 6)) {
            vec3 pos;
            extract_position(&pos.x, &pos.y, &pos.z, line);
            pos_x.push_back(pos.x);
            pos_y.push_back(pos.y);
            pos_z.push_back(pos.z);

            labels.push_back(trim(line.substr(12, 4)));

            if (line.size() > 60) {
                const auto [occupancy, success] = to_float(line.substr(54, 6));
                occupancies.push_back(success ? occupancy : 0.0f);
            }
            if (line.size() > 66) {
                const auto [temp, success] = to_float(line.substr(60, 6));
                temp_factors.push_back(success ? temp : 0.0f);
            }

            Element elem;
            extract_element(&elem, line);
            elements.push_back(elem);

            int res_id = to_int(line.substr(22, 4));
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
                // res.chain_idx = (ChainIdx)(chains.size() - 1);
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
        } else if (compare_n(line, "HELIX", 5)) {
            Range<ResIdx> helix;
            extract_helix_residue_indices(&helix.beg, &helix.end, line);
            helices.push_back(helix);
        } else if (compare_n(line, "SHEET", 5)) { 
            Range<ResIdx> sheet;
            extract_sheet_residue_indices(&sheet.beg, &sheet.end, line);
            sheets.push_back(sheet);
        } else if (compare_n(line, "ENDMDL", 6) || compare_n(line, "END", 3)) {
            break;
        }
    }

    auto masses = compute_atom_masses(elements);
    auto radii = compute_atom_radii(elements);
    auto covalent_bonds = compute_covalent_bonds(residues, pos_x.data(), pos_y.data(), pos_z.data(), elements.data(), num_atoms);
    auto sequences = compute_sequences(residues);
    auto backbone_segments = compute_backbone_segments(residues, labels);
    auto backbone_sequences = compute_backbone_sequences(backbone_segments, residues);
    auto donors = hydrogen_bond::compute_donors(elements, residue_indices, residues, covalent_bonds);
    auto acceptors = hydrogen_bond::compute_acceptors(elements);

    init_molecule_structure(mol, num_atoms, (i32)covalent_bonds.size(), (i32)residues.size(), (i32)chains.size(), (i32)sequences.size(),
                            (i32)backbone_segments.size(), (i32)backbone_sequences.size(), (i32)donors.size(), (i32)acceptors.size());

    for (i32 i = 0; i < num_atoms; i++) {
        mol->atom.res_idx[i] = -1;
        mol->atom.chain_idx[i] = -1;
        mol->atom.seq_idx[i] = -1;
    }

    for (SeqIdx seq_idx = 0; seq_idx < (SeqIdx)sequences.size(); seq_idx++) {
        for (AtomIdx i = sequences[seq_idx].atom_range.beg; i != sequences[seq_idx].atom_range.end; i++) {
            mol->atom.seq_idx[i] = seq_idx;
        }
    }

    for (ChainIdx chain_idx = 0; chain_idx < (ChainIdx)chains.size(); chain_idx++) {
        for (AtomIdx i = chains[chain_idx].atom_range.beg; i != chains[chain_idx].atom_range.end; i++) {
            mol->atom.chain_idx[i] = chain_idx;
        }
    }

    // Copy data into molecule
    memcpy(mol->atom.position.x, pos_x.data(), pos_x.size_in_bytes());
    memcpy(mol->atom.position.y, pos_y.data(), pos_y.size_in_bytes());
    memcpy(mol->atom.position.z, pos_z.data(), pos_z.size_in_bytes());
    memset(mol->atom.velocity.x, 0, num_atoms * sizeof(float));
    memset(mol->atom.velocity.y, 0, num_atoms * sizeof(float));
    memset(mol->atom.velocity.z, 0, num_atoms * sizeof(float));
    memcpy(mol->atom.radius, radii.data(), num_atoms * sizeof(float));
    memcpy(mol->atom.mass, masses.data(), num_atoms * sizeof(float));
    memcpy(mol->atom.element, elements.data(), elements.size_in_bytes());
    memcpy(mol->atom.label, labels.data(), labels.size_in_bytes());
    memcpy(mol->atom.res_idx, residue_indices.data(), residue_indices.size_in_bytes());

    memcpy(mol->residues.data(), residues.data(), residues.size_in_bytes());
    memcpy(mol->chains.data(), chains.data(), chains.size_in_bytes());
    memcpy(mol->sequences.data(), sequences.data(), sequences.size_in_bytes());
    memcpy(mol->covalent_bonds.data(), covalent_bonds.data(), covalent_bonds.size_in_bytes());
    memcpy(mol->backbone.segments.data(), backbone_segments.data(), backbone_segments.size_in_bytes());
    memcpy(mol->backbone.sequences.data(), backbone_sequences.data(), backbone_sequences.size_in_bytes());
    memcpy(mol->hydrogen_bond.donors.data(), donors.data(), donors.size_in_bytes());
    memcpy(mol->hydrogen_bond.acceptors.data(), acceptors.data(), acceptors.size_in_bytes());

    for (const auto& bb_seq : mol->backbone.sequences) {
        auto segments = get_backbone(*mol, bb_seq);
        compute_backbone_angles(mol->backbone.angles.data(), segments.data(), mol->atom.position.x, mol->atom.position.y, mol->atom.position.z,
                                segments.size());
    }

    return true;
}

bool load_trajectory_from_file(MoleculeTrajectory* traj, CStringView filename) {
    StringView pdb_string = allocate_and_read_textfile(filename);
    defer { free_string(&pdb_string); };
    if (!pdb_string) {
        LOG_ERROR("Could not load pdb file");
        return false;
    }
    return load_trajectory_from_string(traj, pdb_string);
}

bool load_trajectory_from_string(MoleculeTrajectory* traj, CStringView pdb_string) {
    ASSERT(traj);
    free_trajectory(traj);

    CStringView pdb_str = pdb_string;
    CStringView mdl_str = extract_next_model(pdb_str);
    if (!mdl_str) {
        LOG_NOTE("Supplied string does not contain MODEL entry and is therefore not a trajectory");
        return false;
    }

    MoleculeInfo info;
    extract_molecule_info(&info, mdl_str);

    if (info.num_atoms == 0) {
        LOG_ERROR("Could not determine number of atoms in trajectory");
        return false;
    }

    // @NOTE: Search space for CRYST1 containing global simulation box parameters
    mat3 sim_box(0);
    CStringView box_str = {pdb_string.beg(), mdl_str.beg() - pdb_string.beg()};
    CStringView line;
    while ((line = extract_line(box_str))) {
        if (compare_n(line, "CRYST1", 6)) {
            extract_simulation_box(&sim_box, line);
            break;
        }
    }

    DynamicArray<CStringView> model_entries;
    model_entries.reserve(1024);

    do {
        model_entries.push_back(mdl_str);
    } while (pdb_str && (mdl_str = extract_next_model(pdb_str)));

    // Time between frames
    const float dt = 1.0f;
    init_trajectory(traj, info.num_atoms, (i32)model_entries.size(), dt, sim_box);
    traj->num_frames = (i32)model_entries.size();

    for (i64 i = 0; i < model_entries.size(); i++) {
        TrajectoryFrame* frame = traj->frame_buffer.data() + i;
        extract_trajectory_frame_data(frame, traj->num_atoms, model_entries[i]);
    }

    return true;
}

bool extract_molecule_info(MoleculeInfo* info, CStringView pdb_string) {
    ASSERT(info);

    i32 num_atoms = 0;
    i32 num_residues = 0;
    i32 num_chains = 0;

    u32 curr_res_pattern = 0;
    char curr_chain_pattern = 0;

    CStringView line;
    while (pdb_string && (line = extract_line(pdb_string))) {
        if (compare_n(line, "ATOM", 4) || compare_n(line, "HETATM", 6)) {
            const u32 res_pattern = *(u32*)(&line[22]);
            const char chain_pattern = line[21];

            num_atoms++;
            if (res_pattern != curr_res_pattern) {
                num_residues++;
                curr_res_pattern = res_pattern;
            }
            if (chain_pattern != curr_chain_pattern) {
                num_chains++;
                curr_chain_pattern = chain_pattern;
            }
        } else if (compare_n(line, "ENDMDL", 6) || compare_n(line, "END", 3)) {
            break;
        }
    }

    info->num_atoms = num_atoms;
    info->num_residues = num_residues;
    info->num_chains = num_chains;

    return true;
}

static bool generate_cache(CStringView filename) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        LOG_ERROR("Could not open file '.*s", filename.length(), filename.beg());
        return false;
    }

    StringBuffer<512> cache_file = get_directory(filename);
    cache_file += "/";
    cache_file += get_file_without_extension(filename);
    cache_file += ".cache";

    constexpr CStringView pattern = "\nMODEL ";
    DynamicArray<FrameBytes> frame_bytes;

    for_each_pattern_found_in_file(filename, pattern, [&frame_bytes](i64 offset) {
        frame_bytes.push_back({(u64)offset, 0});
    });

    fseeki64(file, 0, SEEK_END);
    const u64 file_size = (u64)ftelli64(file);
    fclose(file);

    for (i64 i = 0; i < frame_bytes.size() - 1; i++) {
        frame_bytes[i].extent = frame_bytes[i + 1].offset - frame_bytes[i].offset;
    }
    frame_bytes.back().extent = file_size - frame_bytes.back().offset;

    u64 UID = generate_UID(filename);
    return write_trajectory_cache(UID, frame_bytes.data(), frame_bytes.size(), cache_file);
}

bool read_trajectory_num_frames(i32* num_frames, CStringView filename) {
    u64 UID;
    if ((UID = generate_UID(filename)) != INVALID_UID) {
        StringBuffer<512> cache_file = get_directory(filename);
        cache_file += "/";
        cache_file += get_file_without_extension(filename);
        cache_file += ".cache";

        u64 c_UID;
        i64 c_num_frames;
        if (read_trajectory_cache_header(&c_UID, &c_num_frames, cache_file) && (UID == c_UID)) {
            // Cache is valid
            *num_frames = (i32)c_num_frames;
            return true;
        }
    }

    i32 nframes = 0;
    bool res = for_each_pattern_found_in_file(filename, "\nMODEL ", [&nframes](i64 offset) {
        (void)offset;
        nframes++;
    });

    if (res) {
        *num_frames = nframes;
        return true;
    }
    
    return false;
}

bool read_trajectory_frame_bytes(FrameBytes* frame_bytes, CStringView filename) {
    u64 UID;
    if ((UID = generate_UID(filename)) != INVALID_UID) {
        StringBuffer<512> cache_file = get_directory(filename);
        cache_file += "/";
        cache_file += get_file_without_extension(filename);
        cache_file += ".cache";

        u64 c_UID;
        i64 c_num_frames;
        if (read_trajectory_cache_header(&c_UID, &c_num_frames, cache_file) && (UID == c_UID)) {
            // Cache is valid
            return read_trajectory_cache(frame_bytes, cache_file);
        } else {
            // Cache is invalid
            // Regenerate data
            if (generate_cache(filename)) {
                if (read_trajectory_cache_header(&c_UID, &c_num_frames, cache_file) && (UID == c_UID)) {
                    return read_trajectory_cache(frame_bytes, cache_file);
                }
            }
        }
    }
    return false;
}

bool read_trajectory_simulation_box(mat3* sim_box, CStringView filename) {
    FILE* file = fopen(filename, "rb");
    if (file) {
        constexpr i64 buf_size = MEGABYTES(1);
        char* buf = (char*)TMP_MALLOC(buf_size);
        ASSERT(buf);
        const i64 bytes_read = fread(buf, 1, buf_size, file);
        CStringView line = find_string({buf, bytes_read}, "CRYST1");
        if (line) {
            line.count = 54;
            extract_simulation_box(sim_box, line);
            return true;
        }
    }
    LOG_ERROR("Could not read file '.*s'", filename.length(), filename.beg());
    return false;
}

bool extract_trajectory_frame(TrajectoryFrame* frame, i32 num_atoms, Array<u8> data) {
    return extract_trajectory_frame_data(frame, num_atoms, {(char*)data.beg(), (char*)data.end()});
}

}  // namespace pdb
