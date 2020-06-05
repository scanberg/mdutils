#include "gro_utils.h"
#include <core/string_utils.h>
#include <core/log.h>
#include <mol/element.h>
#include <mol/element_utils.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>

#include <stdio.h>

namespace gro {

struct LineFormat {
    int pos_width = -1;
};

bool load_molecule_from_file(MoleculeStructure* mol, CStringView filename) {
    StringView txt = allocate_and_read_textfile(filename);
    defer { FREE(txt.cstr()); };
    if (!txt) {
        LOG_ERROR("Could not read file: '%.*s'.", filename.length(), filename);
        return false;
    }
    return load_molecule_from_string(mol, txt);
}

inline LineFormat get_format(CStringView line) {
    LineFormat fmt{};

    // First float starts at offset 20, count
    if (line.length() > 20) {
        const char* c = line.beg() + 20;
        while (c != line.end() && *c != '\n' && *c == ' ') c++;
        while (c != line.end() && *c != '\n' && *c != ' ') c++;
        if (c != line.end()) {
            fmt.pos_width = (int)(c - (&line[20]));
        }
    }
    return fmt;
}

inline void extract_position_data(float* x, float* y, float* z, CStringView line, int width) {
    *x = str_to_float(line.substr(20 + 0 * width, width)) * 10.0f; // nm -> Å
    *y = str_to_float(line.substr(20 + 1 * width, width)) * 10.0f;
    *z = str_to_float(line.substr(20 + 2 * width, width)) * 10.0f;
}

inline void extract_box_ext_data(float* x, float* y, float* z, CStringView line) {
    *x = str_to_float(line.substr(0, 8));
    *y = str_to_float(line.substr(10, 8));
    *z = str_to_float(line.substr(20, 8));
}

bool load_molecule_from_string(MoleculeStructure* mol, CStringView gro_string) {
    CStringView header = extract_line(gro_string);
    CStringView length = extract_line(gro_string);
    (void)header;

    int num_atoms = to_int(length);
    if (num_atoms == 0) {
        return false;
    }

    const LineFormat format = get_format(peek_line(gro_string));
    if (format.pos_width == -1) {
        LOG_ERROR("Could not identify internal line format of gro file!");
        return false;
    }

    DynamicArray<AtomDescriptor> atoms;
    DynamicArray<ResidueDescriptor> residues;
    atoms.reserve(num_atoms);

    int cur_res = -1;
    CStringView line;

    for (int i = 0; i < num_atoms; ++i) {
        line = extract_line(gro_string);
        
        AtomDescriptor& atom = atoms.allocate_back();
        extract_position_data(&atom.x, &atom.y, &atom.z, line, format.pos_width);
        atom.name = trim(line.substr(10, 5));
        
        int res_id = to_int(line.substr(0, 5));
        if (cur_res != res_id) { 
            ResidueDescriptor& res = residues.allocate_back();
            res.name = trim(line.substr(5, 5));
            res.id = res_id;
            res.atom_range = {i, i};
            cur_res = res_id;
        }

        atoms.back().residue_index = (ResIdx)residues.size() - 1;
        residues.back().atom_range.end++;
    }

    MoleculeStructureDescriptor desc;
    desc.num_atoms = atoms.size();
    desc.atoms = atoms.data();
    desc.num_residues = residues.size();
    desc.residues = residues.data();

    return init_molecule_structure(mol, desc);
}

}  // namespace gro
