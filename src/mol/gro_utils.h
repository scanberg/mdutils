#pragma once

#include <mol/molecule_structure.h>
#include <core/string_utils.h>

namespace gro {
bool load_molecule_from_file(MoleculeStructure* mol, CString filename);
bool load_molecule_from_string(MoleculeStructure* mol, CString string);
}
