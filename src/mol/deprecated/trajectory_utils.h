#pragma once

#include <core/string_types.h>
#include <mol/molecule_trajectory.h>

struct FrameBytes {
    u64 offset : 40; // max: ~1.1TB
    u64 extent : 24; // max: 16MB
};

constexpr u64 INVALID_UID = 0;

u64 generate_UID(CStringView filename);

// Reads the number of frames and unique ID of trajectory
bool read_trajectory_cache_header(u64* UID, i64* num_frames, CStringView cache_filename);

// Read trajectory Frame Byte Cache and the unique ID
bool read_trajectory_cache(FrameBytes* frame_bytes, CStringView cache_filename);

// Write trajectory Frame Byte Cache with a unique ID (fingerprint)
bool write_trajectory_cache(u64 UID, const FrameBytes* frame_bytes, i64 num_frames, CStringView cache_filename);


inline TrajectoryFrame& get_trajectory_frame(MoleculeTrajectory& traj, int frame_index) {
    ASSERT(-1 < frame_index && frame_index < traj.num_frames);
    return traj.frame_buffer.ptr[frame_index];
}

inline const TrajectoryFrame& get_trajectory_frame(const MoleculeTrajectory& traj, int frame_index) {
    ASSERT(-1 < frame_index && frame_index < traj.num_frames);
    return traj.frame_buffer.ptr[frame_index];
}

inline soa_vec3 get_trajectory_positions(MoleculeTrajectory& traj, int frame_index) {
    ASSERT(0 <= frame_index && frame_index < traj.num_frames);
    return traj.frame_buffer.ptr[frame_index].atom_position;
}

inline const soa_vec3 get_trajectory_positions(const MoleculeTrajectory& traj, int frame_index) {
    ASSERT(0 <= frame_index && frame_index < traj.num_frames);
    return traj.frame_buffer.ptr[frame_index].atom_position;
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