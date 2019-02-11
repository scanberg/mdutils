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

    DynamicArray<vec3> positions;
    DynamicArray<vec3> velocities;
    DynamicArray<Label> labels;
    DynamicArray<Element> elements;
    DynamicArray<ResIdx> residue_indices;
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
        vec3 pos, vel;
        int atom_idx, res_idx;
        char atom_name[8] = {};
        char res_name[8] = {};

        line = extract_line(gro_string);  // line becomes zero terminated upon assignment
        auto result = sscanf(line, format, &res_idx, res_name, atom_name, &atom_idx, &pos.x, &pos.y, &pos.z, &vel.x, &vel.y, &vel.z);
        if (result > 0) {
            if (cur_res != res_idx) {
                cur_res = res_idx;
                res_count = (int)residues.count;
                CString res_name_trim = trim(CString(res_name));
                Residue res{};
                res.name = res_name_trim;
                res.id = res_idx;
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

            positions.push_back(pos);
            velocities.push_back(vel);
            labels.push_back(atom_name_trim);
            elements.push_back(elem);
            residue_indices.push_back((ResIdx)res_count);
        }
    }

    vec3 box{};
    line = extract_line(gro_string);
    sscanf(line, "%8f %8f %8f", &box.x, &box.y, &box.z);

    // Convert from nm to ångström
    for (auto& p : positions) {
        p *= 10.f;
    }
    for (auto& v : velocities) {
        v *= 10.f;
    }
    box *= 10.f;

    auto covalent_bonds = compute_covalent_bonds(residues, residue_indices, positions, elements);
    auto backbone_segments = compute_backbone_segments(residues, labels);
    auto backbone_sequences = compute_backbone_sequences(backbone_segments, residues);
    auto backbone_angles = compute_backbone_angles(positions, backbone_segments, backbone_sequences);
    auto chains = compute_chains(residues);
    auto donors = hydrogen_bond::compute_donors(elements, residue_indices, residues, covalent_bonds);
    auto acceptors = hydrogen_bond::compute_acceptors(elements);

    for (ChainIdx c = 0; c < chains.count; c++) {
        for (auto i = chains[c].res_idx.beg; i < chains[c].res_idx.end; i++) {
            residues[i].chain_idx = c;
        }
    }

    init_molecule_structure(mol, num_atoms, (int32)covalent_bonds.size(), (int32)residues.size(), (int32)chains.size(), (int32)backbone_segments.size(), (int32)backbone_sequences.size(),
                            (int32)donors.size(), (int32)acceptors.size());

    // Copy data into molecule
    memcpy(mol->atom.positions, positions.ptr, positions.size_in_bytes());
    memcpy(mol->atom.elements, elements.ptr, elements.size_in_bytes());
    memcpy(mol->atom.labels, labels.ptr, labels.size_in_bytes());
    memcpy(mol->atom.residue_indices, residue_indices.ptr, residue_indices.size_in_bytes());

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
