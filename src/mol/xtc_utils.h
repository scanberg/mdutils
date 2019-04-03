#pragma once

#include <mol/molecule_trajectory.h>

namespace xtc {

constexpr uint32 XTC_FILE_TAG = 0x50002;

bool init_trajectory_from_file(MoleculeTrajectory* traj, int32 mol_atom_count, CString filename);
bool close_file_handle(MoleculeTrajectory* traj);

bool read_next_trajectory_frame(MoleculeTrajectory* traj);
}

