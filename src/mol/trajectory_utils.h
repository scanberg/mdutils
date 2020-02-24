#pragma once

#include <mol/molecule_trajectory.h>

inline bool all_trajectory_frames_read(const MoleculeTrajectory& traj) { return (traj.num_frames == (i32)traj.frame_offsets.count); }

bool read_next_trajectory_frame(MoleculeTrajectory* traj);

bool close_file_handle(MoleculeTrajectory* traj);

inline TrajectoryFrame& get_trajectory_frame(MoleculeTrajectory& traj, int frame_index) {
    ASSERT(-1 < frame_index && frame_index < traj.num_frames);
    return traj.frame_buffer.ptr[frame_index];
}

inline const TrajectoryFrame& get_trajectory_frame(const MoleculeTrajectory& traj, int frame_index) {
    ASSERT(-1 < frame_index && frame_index < traj.num_frames);
    return traj.frame_buffer.ptr[frame_index];
}

inline Array<float> get_trajectory_position_x(MoleculeTrajectory& traj, int frame_index) {
    ASSERT(-1 < frame_index && frame_index < traj.num_frames);
    return {traj.frame_buffer.ptr[frame_index].atom_position.x, traj.num_atoms};
}
inline Array<const float> get_trajectory_position_x(const MoleculeTrajectory& traj, int frame_index) {
    ASSERT(-1 < frame_index && frame_index < traj.num_frames);
    return {traj.frame_buffer.ptr[frame_index].atom_position.x, traj.num_atoms};
}

inline Array<float> get_trajectory_position_y(MoleculeTrajectory& traj, int frame_index) {
    ASSERT(-1 < frame_index && frame_index < traj.num_frames);
    return {traj.frame_buffer.ptr[frame_index].atom_position.y, traj.num_atoms};
}
inline Array<const float> get_trajectory_position_y(const MoleculeTrajectory& traj, int frame_index) {
    ASSERT(-1 < frame_index && frame_index < traj.num_frames);
    return {traj.frame_buffer.ptr[frame_index].atom_position.y, traj.num_atoms};
}

inline Array<float> get_trajectory_position_z(MoleculeTrajectory& traj, int frame_index) {
    ASSERT(-1 < frame_index && frame_index < traj.num_frames);
    return {traj.frame_buffer.ptr[frame_index].atom_position.z, traj.num_atoms};
}
inline Array<const float> get_trajectory_position_z(const MoleculeTrajectory& traj, int frame_index) {
    ASSERT(-1 < frame_index && frame_index < traj.num_frames);
    return {traj.frame_buffer.ptr[frame_index].atom_position.z, traj.num_atoms};
}