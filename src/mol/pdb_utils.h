#pragma once

#include <core/types.h>
#include <core/string_types.h>
#include <mol/molecule_dynamic.h>

namespace pdb {

constexpr uint32 PDB_FILE_TAG = 0x50001;

bool load_molecule_from_file(MoleculeStructure* mol, CString filename);
bool load_molecule_from_string(MoleculeStructure* mol, CString string);

bool load_dynamic_from_file(MoleculeDynamic* md, CString filename);
bool load_dynamic_from_string(MoleculeDynamic* md, CString string);

struct Info {
    int32 num_atoms = 0;
    int32 num_residues = 0;
    int32 num_chains = 0;
    int32 num_frames = 0;
    DynamicArray<Range<int64>> frame_byte_ranges = {};
};

bool extract_info(Info* info, CString pdb_string);

bool init_dynamic_from_file(MoleculeDynamic* md, CString filename);

bool read_next_trajectory_frame(MoleculeTrajectory* traj);
bool close_file_handle(MoleculeTrajectory* traj);

}