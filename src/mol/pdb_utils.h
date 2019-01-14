#pragma once

#include <mol/molecule_dynamic.h>
#include <core/common.h>
#include <core/string_utils.h>

enum PdbLoadParams {
	PDB_READ_ATOM = BIT(0),
	PDB_READ_HETATM = BIT(1),
	PDB_TREAT_MODELS_AS_FRAMES = BIT(2),
	PDB_DEFAULT = 0xFFFFFFFF
};

bool allocate_and_load_pdb_from_file(MoleculeDynamic* md, CString filename, PdbLoadParams params = PDB_DEFAULT);
bool allocate_and_parse_pdb_from_string(MoleculeDynamic* md, CString string, PdbLoadParams params = PDB_DEFAULT);

struct Result {
    int32 num_frames = 0;
};

Result load_pdb(CString filename);