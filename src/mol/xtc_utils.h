#pragma once

#include <mol/molecule_trajectory.h>

namespace xtc {

constexpr u32 XTC_FILE_TAG = 0x50002;

bool init_trajectory_from_file(MoleculeTrajectory* traj, i32 num_atoms, CStringView filename);
bool close_file_handle(MoleculeTrajectory* traj);

bool read_next_trajectory_frame(MoleculeTrajectory* traj);
bool read_trajectory_frame_from_memory(TrajectoryFrame* frame, i32 num_atoms, void* ptr, int64_t num_bytes);
}