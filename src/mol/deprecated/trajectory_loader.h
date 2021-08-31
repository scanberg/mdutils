#pragma once

#include <stdint.h>

struct Trajectory;
struct TrajectoryFrame;

struct TrajectoryLoader {
    bool (*initialize_trajectory)(Trajectory* traj, int32_t num_atoms, int32_t num_frames) = 0;
    bool (*load_frame_data)(TrajectoryFrame* frame, void* ptr, int64_t num_bytes) = 0;
    bool (*free_trajectory)(Trajectory* traj) = 0;
};


TrajectoryLoader* get_loader(uint64_t id);