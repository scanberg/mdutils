#pragma once

#include <mol/molecule_trajectory.h>

// Reads the header info of a trajectory and allocates space for it
// Only supports XTC currently
bool load_and_allocate_trajectory(MoleculeTrajectory* traj, CString path);

// Reads the actual trajectory position information... Necessary?
bool read_trajectory_data(MoleculeTrajectory* traj);

bool read_next_trajectory_frame(MoleculeTrajectory* traj);
bool all_trajectory_frames_read(MoleculeTrajectory* traj);
bool close_file_handle(MoleculeTrajectory* traj);

//void copy_trajectory_frame(TrajectoryFrame* dst, const MoleculeTrajectory& src_traj, int frame_index);
//void copy_trajectory_positions(Array<vec3> dst_array, const MoleculeTrajectory& traj, int frame_index);
//void read_trajectory_box_vectors(vec3 box_vectors[3], const MoleculeTrajectory& traj, int frame_index);

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
	return { traj.frame_buffer.ptr[frame_index].atom_position.x, traj.num_atoms };
}
inline Array<const float> get_trajectory_position_x(const MoleculeTrajectory& traj, int frame_index) {
    ASSERT(-1 < frame_index && frame_index < traj.num_frames);
	return { traj.frame_buffer.ptr[frame_index].atom_position.x, traj.num_atoms };
}

inline Array<float> get_trajectory_position_y(MoleculeTrajectory& traj, int frame_index) {
	ASSERT(-1 < frame_index && frame_index < traj.num_frames);
	return { traj.frame_buffer.ptr[frame_index].atom_position.y, traj.num_atoms };
}
inline Array<const float> get_trajectory_position_y(const MoleculeTrajectory& traj, int frame_index) {
	ASSERT(-1 < frame_index && frame_index < traj.num_frames);
	return { traj.frame_buffer.ptr[frame_index].atom_position.y, traj.num_atoms };
}

inline Array<float> get_trajectory_position_z(MoleculeTrajectory& traj, int frame_index) {
	ASSERT(-1 < frame_index && frame_index < traj.num_frames);
	return { traj.frame_buffer.ptr[frame_index].atom_position.z, traj.num_atoms };
}
inline Array<const float> get_trajectory_position_z(const MoleculeTrajectory& traj, int frame_index) {
	ASSERT(-1 < frame_index && frame_index < traj.num_frames);
	return { traj.frame_buffer.ptr[frame_index].atom_position.z, traj.num_atoms };
}