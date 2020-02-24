#include "molecule_trajectory.h"
#include <core/log.h>
#include <core/hash.h>
#include <mol/trajectory_utils.h>

#define ALIGNMENT 64

bool init_trajectory(MoleculeTrajectory* traj, i32 num_atoms, i32 num_frames, f32 time_between_frames, const mat3& sim_box) {
    ASSERT(traj);

    traj->num_atoms = num_atoms;
    traj->num_frames = num_frames;
    traj->total_simulation_time = 0;
    traj->simulation_type = SimulationType::Undefined;
	traj->file = {};

    const i64 pos_mem_size = (num_frames * num_atoms * sizeof(float) + ALIGNMENT) * 3;
    void* pos_mem = ALIGNED_MALLOC(pos_mem_size, ALIGNMENT);

    const i64 frame_mem_size = num_frames * sizeof(TrajectoryFrame);
    void* frame_mem = MALLOC(frame_mem_size);

    if (!pos_mem) {
        LOG_ERROR("Could not allocate memory for trajectory positions");
        return false;
    }

    if (!frame_mem) {
        LOG_ERROR("Could not allocate memory for trajectory frames");
        return false;
    }

    float* pos_data_x = (float*)pos_mem;
    float* pos_data_y = (float*)get_next_aligned_adress(pos_data_x + num_frames * num_atoms, ALIGNMENT);
    float* pos_data_z = (float*)get_next_aligned_adress(pos_data_y + num_frames * num_atoms, ALIGNMENT);

    ASSERT(IS_ALIGNED(pos_data_x, ALIGNMENT));
    ASSERT(IS_ALIGNED(pos_data_y, ALIGNMENT));
    ASSERT(IS_ALIGNED(pos_data_z, ALIGNMENT));

    traj->frame_offsets = {};
    traj->position_data.x = pos_data_x;
    traj->position_data.y = pos_data_y;
    traj->position_data.z = pos_data_z;
    traj->frame_buffer = {(TrajectoryFrame*)frame_mem, num_frames};

    for (i32 i = 0; i < num_frames; i++) {
        traj->frame_buffer[i].index = i;
        traj->frame_buffer[i].time = i * time_between_frames;
        traj->frame_buffer[i].box = sim_box;
        traj->frame_buffer[i].atom_position.x = traj->position_data.x + i * num_atoms;
        traj->frame_buffer[i].atom_position.y = traj->position_data.y + i * num_atoms;
        traj->frame_buffer[i].atom_position.z = traj->position_data.z + i * num_atoms;
    }

    return true;
}

void free_trajectory(MoleculeTrajectory* traj) {
    ASSERT(traj);

	close_file_handle(traj);
    traj->file.path = "";
    if (traj->frame_offsets.ptr) FREE(traj->frame_offsets.ptr);
    if (traj->position_data.x) ALIGNED_FREE(traj->position_data.x);
    if (traj->frame_buffer.ptr) FREE(traj->frame_buffer.ptr);

    *traj = {};
}

