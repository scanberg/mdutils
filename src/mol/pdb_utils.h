#pragma once

#include <core/types.h>
#include <core/string_types.h>
#include <mol/molecule_dynamic.h>

namespace pdb {

constexpr uint32 PDB_FILE_TAG = 0x50001;

struct MoleculeInfo {
	int32 num_atoms = 0;
	int32 num_residues = 0;
	int32 num_chains = 0;
};

// Loads molecule data
bool load_molecule_from_file(MoleculeStructure* mol, CString filename);
bool load_molecule_from_string(MoleculeStructure* mol, CString string);

// Loads entire trajectory
bool load_trajectory_from_file(MoleculeTrajectory* traj, CString filename);
bool load_trajectory_from_string(MoleculeTrajectory* traj, CString string);

// Initializes a trajectory from file, but does not load frames (for async operations)
bool init_trajectory_from_file(MoleculeTrajectory* traj, CString filename);
bool read_next_trajectory_frame(MoleculeTrajectory* traj);
bool close_file_handle(MoleculeTrajectory* traj);

// Extract molecule info from a pdb string
bool extract_molecule_info(MoleculeInfo* info, CString pdb_string);

}