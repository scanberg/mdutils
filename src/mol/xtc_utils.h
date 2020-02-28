#pragma once

#include <mol/molecule_trajectory.h>

namespace xtc {

constexpr u32 XTC_FILE_TAG = 0x50002;

bool init_trajectory_from_file(MoleculeTrajectory* traj, i32 num_atoms, CStringView filename);
bool close_file_handle(MoleculeTrajectory* traj);

bool read_next_trajectory_frame(MoleculeTrajectory* traj);
bool read_trajectory_frame_from_memory(TrajectoryFrame* frame, i32 num_atoms, void* ptr, int64_t num_bytes);

i32 read_num_frames(CStringView filename);
DynamicArray<i64> read_frame_offsets(CStringView filename);
bool decompress_trajectory_frame(TrajectoryFrame* frame, i32 num_atoms, Array<u8> raw_data);

}