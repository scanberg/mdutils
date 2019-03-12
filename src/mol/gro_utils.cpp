#include "gro_utils.h"
#include <core/string_utils.h>
#include <core/log.h>
#include <mol/element.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>

bool allocate_and_load_gro_from_file(MoleculeStructure* mol, CString filename) {
    String txt = allocate_and_read_textfile(filename);
    defer { FREE(txt); };
    if (!txt) {
        LOG_ERROR("Could not read file: '%.*s'.", filename.length(), filename);
        return false;
    }
    auto res = allocate_and_parse_gro_from_string(mol, txt);
    return res;
}

const char* get_format(CString line) {
    const char* format_n = "%5d%5c%5c%5d%8f%8f%8f%8f%8f%8f";        // narrow 8 characters per float
    const char* format_w = "%5d%5c%5c%5d%11f%11f%11f%11f%11f%11f";  // wide 11 characters per float, perhaps meant to be 10 with a space character.. who knows!

    // First float starts at offset 20, count
    if (line.length() < 20) {
        return nullptr;
    }

    const uint8* c = &line[20];
    while (c != line.end() && *c == ' ') c++;
    while (c != line.end() && *c != ' ') c++;

    auto len = c - (&line[20]);
    if (len < 10) {
        return format_n;
    } else {
        return format_w;
    }
}

bool allocate_and_parse_gro_from_string(MoleculeStructure* mol, CString gro_string) {
    CString header = extract_line(gro_string);
    CString length = extract_line(gro_string);
    (void)header;

    int num_atoms = to_int(length);
    if (num_atoms == 0) {
        return false;
    }

    int64 mem_size = (sizeof(float) * (3 + 3 + 1) + sizeof(Label) + sizeof(Element) + sizeof(ResIdx)) * num_atoms;
    void* mem = TMP_MALLOC(mem_size);
    defer { TMP_FREE(mem); };
    memset(mem, 0, mem_size);

    float* atom_pos_x = (float*)mem;
    float* atom_pos_y = (float*)(atom_pos_x + num_atoms);
    float* atom_pos_z = (float*)(atom_pos_y + num_atoms);
    float* atom_vel_x = (float*)(atom_pos_z + num_atoms);
    float* atom_vel_y = (float*)(atom_vel_x + num_atoms);
    float* atom_vel_z = (float*)(atom_vel_y + num_atoms);
    float* atom_radius = (float*)(atom_vel_z + num_atoms);
    Label* atom_label = (Label*)(atom_radius + num_atoms);
    Element* atom_element = (Element*)(atom_label + num_atoms);
    ResIdx* atom_res_idx = (ResIdx*)(atom_element + num_atoms);

    DynamicArray<Residue> residues;

    const char* format = get_format(peek_line(gro_string));
    if (!format) {
        LOG_ERROR("Could not identify internal format of gro file!");
        return false;
    }
    int res_count = 0;
    int cur_res = -1;
    StringBuffer<256> line;

    for (int i = 0; i < num_atoms; ++i) {
        float pos[3] = {0, 0, 0};
        float vel[3] = {0, 0, 0};
        int atom_idx, res_id;
        char atom_name[8] = {};
        char res_name[8] = {};

        line = extract_line(gro_string);
        int result = sscanf(line, format, &res_id, res_name, atom_name, &atom_idx, &pos[0], &pos[1], &pos[2], &vel[0], &vel[1], &vel[2]);
        if (result > 0) {
            if (cur_res != res_id) {
                cur_res = res_id;
                res_count = (int)residues.count;
                CString res_name_trim = trim(CString(res_name));
                Residue res{};
                res.name = res_name_trim;
                res.id = res_id;
                res.chain_idx = 0;
                res.atom_idx = {i, i};
                residues.push_back(res);
            }
            residues.back().atom_idx.end++;

            CString atom_name_trim = trim(CString(atom_name));
            CString element_str = atom_name_trim;

            if (is_amino_acid(residues.back())) {
                // If we have an amino acid, we can assume its an organic element with just one letter. C/N/H/O?
                element_str = element_str.substr(0, 1);
            }
            Element elem = element::get_from_string(element_str);

            // Convert from nm to ångström
            atom_pos_x[i] = pos[0] * 10.f;
            atom_pos_y[i] = pos[1] * 10.f;
            atom_pos_z[i] = pos[2] * 10.f;
            atom_vel_x[i] = vel[0] * 10.f;
            atom_vel_y[i] = vel[1] * 10.f;
            atom_vel_z[i] = vel[2] * 10.f;

            atom_label[i] = atom_name_trim;
            atom_element[i] = elem;
            atom_res_idx[i] = res_count;
        }
    }

    vec3 box{};
    line = extract_line(gro_string);
    sscanf(line, "%8f %8f %8f", &box.x, &box.y, &box.z);

    // Convert from nm to ångström
    box *= 10.f;

    compute_atom_radii(atom_radius, atom_element, num_atoms);
    auto covalent_bonds = compute_covalent_bonds(residues, atom_pos_x, atom_pos_y, atom_pos_z, atom_res_idx, atom_element, num_atoms);
    auto backbone_segments = compute_backbone_segments(residues, {atom_label, num_atoms});
    auto backbone_sequences = compute_backbone_sequences(backbone_segments, residues);
    auto backbone_angles = compute_backbone_angles(backbone_segments, backbone_sequences, atom_pos_x, atom_pos_y, atom_pos_z);
    auto chains = compute_chains(residues);
    auto donors = hydrogen_bond::compute_donors({atom_element, num_atoms}, {atom_res_idx, num_atoms}, residues, covalent_bonds);
    auto acceptors = hydrogen_bond::compute_acceptors({atom_element, num_atoms});

    for (ChainIdx c = 0; c < chains.count; c++) {
        for (auto i = chains[c].res_idx.beg; i < chains[c].res_idx.end; i++) {
            residues[i].chain_idx = c;
        }
    }

    init_molecule_structure(mol, num_atoms, (int32)covalent_bonds.size(), (int32)residues.size(), (int32)chains.size(), (int32)backbone_segments.size(), (int32)backbone_sequences.size(),
                            (int32)donors.size(), (int32)acceptors.size());

    // Copy data into molecule
    memcpy(mol->atom.position.x, atom_pos_x, sizeof(float) * num_atoms);
    memcpy(mol->atom.position.y, atom_pos_y, sizeof(float) * num_atoms);
    memcpy(mol->atom.position.z, atom_pos_z, sizeof(float) * num_atoms);

    memcpy(mol->atom.velocity.x, atom_vel_x, sizeof(float) * num_atoms);
    memcpy(mol->atom.velocity.y, atom_vel_y, sizeof(float) * num_atoms);
    memcpy(mol->atom.velocity.z, atom_vel_z, sizeof(float) * num_atoms);

    memcpy(mol->atom.radius, atom_radius, sizeof(float) * num_atoms);

    memcpy(mol->atom.element, atom_element, sizeof(Element) * num_atoms);
    memcpy(mol->atom.label, atom_label, sizeof(Label) * num_atoms);
    memcpy(mol->atom.res_idx, atom_res_idx, sizeof(ResIdx) * num_atoms);

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
