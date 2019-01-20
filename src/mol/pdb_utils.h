#pragma once

#include <mol/molecule_dynamic.h>
#include <core/common.h>
#include <core/string_utils.h>

bool allocate_and_load_pdb_from_file(MoleculeDynamic* md, CString filename);
bool allocate_and_parse_pdb_from_string(MoleculeDynamic* md, CString string);