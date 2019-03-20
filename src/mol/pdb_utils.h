#pragma once

#include <core/types.h>
#include <core/string_types.h>
#include <mol/molecule_dynamic.h>

bool allocate_and_load_pdb_from_file(MoleculeDynamic* md, CString filename);
bool allocate_and_parse_pdb_from_string(MoleculeDynamic* md, CString string);

struct PdbInfo {
    int32 num_atoms = 0;
    int32 num_residues = 0;
    int32 num_chains = 0;
    int32 num_frames = 0;
    DynamicArray<Range<int64>> frame_byte_ranges = {};
};

bool extract_pdb_info(PdbInfo* info, CString pdb_string);
