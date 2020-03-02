#pragma once

#include <core/types.h>
#include <core/string_types.h>
#include <mol/molecule_dynamic.h>
#include <mol/trajectory_utils.h>

namespace pdb {

struct MoleculeInfo {
	i32 num_atoms = 0;
	i32 num_residues = 0;
	i32 num_chains = 0;
};

// Loads molecule data
bool load_molecule_from_file(MoleculeStructure* mol, CStringView filename);
bool load_molecule_from_string(MoleculeStructure* mol, CStringView string);

// Loads entire trajectory
bool load_trajectory_from_file(MoleculeTrajectory* traj, CStringView filename);
bool load_trajectory_from_string(MoleculeTrajectory* traj, CStringView string);

// Extract molecule info from a pdb string
bool extract_molecule_info(MoleculeInfo* info, CStringView pdb_string);

// --- Core Trajectory Functionality ---
bool read_trajectory_num_frames(i32* num_frames, CStringView filename);
bool read_trajectory_simulation_box(mat3* sim_box, CStringView filename);

// Reads byte offset and length of frames within pdb file
bool read_trajectory_frame_bytes(FrameBytes* frame_bytes, CStringView filename);

// Extracts trajectory frame data from a raw-chunk of pdb_data
bool extract_trajectory_frame(TrajectoryFrame* frame, i32 num_atoms, Array<u8> data);
}