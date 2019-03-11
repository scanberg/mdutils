#include "molecule_trajectory.h"

#define ALIGNMENT 64

bool init_trajectory(MoleculeTrajectory* traj, int32 num_atoms, int32 num_frames) {
    ASSERT(traj);

    traj->num_atoms = num_atoms;
    traj->num_frames = num_frames;
    traj->total_simulation_time = 0;
    traj->simulation_type = MoleculeTrajectory::NVT;
    traj->path_to_file = {};
    traj->file_handle = nullptr;

    float* pos_data = (float*)ALIGNED_MALLOC((num_frames * num_atoms * sizeof(float) + ALIGNMENT) * 3, ALIGNMENT);
    float* pos_data_x = pos_data;
    float* pos_data_y = (float*)get_next_aligned_adress(pos_data_x + num_frames * num_atoms, ALIGNMENT);
    float* pos_data_z = (float*)get_next_aligned_adress(pos_data_y + num_frames * num_atoms, ALIGNMENT);

    traj->frame_offsets = {};
    traj->position_data = {pos_data_x, pos_data_y, pos_data_z};
    traj->frame_buffer = {(TrajectoryFrame*)CALLOC(num_frames, sizeof(TrajectoryFrame)), num_frames};

    return true;
}

void free_trajectory(MoleculeTrajectory* traj) {
    ASSERT(traj);

    free_string(&traj->path_to_file);
    if (traj->frame_offsets.ptr) FREE(traj->frame_offsets.ptr);
    if (traj->position_data.x) ALIGNED_FREE(traj->position_data.x);
    if (traj->frame_buffer.ptr) FREE(traj->frame_buffer.ptr);

    *traj = {};
}
