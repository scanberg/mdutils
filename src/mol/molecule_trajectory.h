#pragma once

#include <core/types.h>
#include <core/array_types.h>
#include <core/vector_types.h>
#include <core/common.h>

enum class SimulationType {
    Undefined,
    NVT,
    NPT,
    NVE
};

struct TrajectoryFrame {
    i32 index = 0;
    f32 time = 0;
    mat3 box = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    struct {
        float* x = nullptr;
        float* y = nullptr;
        float* z = nullptr;
    } atom_position;
};

struct MoleculeTrajectory {
    i32 num_atoms = 0;
    i32 num_frames = 0;
    f32 total_simulation_time = 0;
    SimulationType simulation_type = SimulationType::Undefined;

    // @NOTE: The frame_buffer may not contain all frames in trajectory.
    // If the trajectory is large, frame_buffer will be used as a cache towards the trajectory streamed from disk.
    // @NOTE: This is not implemented at the moment.
    Array<TrajectoryFrame> frame_buffer{};

    // This is the position data of the full trajectory
    struct {
        float* x = nullptr;
        float* y = nullptr;
        float* z = nullptr;
    } position_data;

    // These are the offsets for each frame inside the file on disk.
    //Array<i64> frame_offsets{};

    operator bool() const { return num_atoms > 0 && frame_buffer.count > 0; }
};

// Allocates memory and initializes trajectory
bool init_trajectory(MoleculeTrajectory* traj, i32 num_atoms, i32 num_frames, f32 time_between_frames = 1.0f, const mat3& sim_box = {});

// Frees memory allocated by trajectory
void free_trajectory(MoleculeTrajectory* traj);