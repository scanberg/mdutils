#pragma once

#include <mol/molecule_trajectory.h>
#include <mol/trajectory_utils.h>

namespace xtc {

// Helper functions
bool read_trajectory_frames_from_file(Array<TrajectoryFrame> frames, Array<const i64> frame_file_offsets, i32 num_atoms, CStringView filename);

// --- Core functionality ---
bool read_trajectory_num_frames(i32* num_frames, CStringView filename);
bool read_trajectory_frame_bytes(FrameBytes* frame_bytes, CStringView filename);
bool decompress_trajectory_frame(TrajectoryFrame* frame, i32 num_atoms, Array<u8> raw_data);

}