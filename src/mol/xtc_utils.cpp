#include "xtc_utils.h"
#include <core/string_utils.h>
#include <core/log.h>

#include <stdio.h>
#include <xdrfile_xtc.h>

namespace xtc {
bool init_trajectory_from_file(MoleculeTrajectory* traj, int32 mol_atom_count, CString filename) {
    ASSERT(traj);
    free_trajectory(traj);

    CString directory = get_directory(filename);
    CString file = get_file_without_extension(filename);
    StringBuffer<512> cache_file = directory;
    cache_file += "/";
    cache_file += file;
    cache_file += ".cache";

    XDRFILE* file_handle = xdrfile_open(filename.cstr(), "r");
    if (!file_handle) {
        LOG_ERROR("Could not open file: %s", filename);
        return false;
    }

    int32 num_atoms = 0;
    int32 num_frames = 0;
    int64* offsets = nullptr;
    if (read_xtc_natoms(filename.cstr(), &num_atoms) != exdrOK) {
        LOG_ERROR("Could not extract number of atoms in trajectory");
		xdrfile_close(file_handle);
        return false;
    }

	if (num_atoms != mol_atom_count) {
		LOG_ERROR("Trajectory atom count did not match molecule atom count");
		xdrfile_close(file_handle);
		return false;
	}

    FILE* offset_cache_handle = fopen(cache_file.cstr(), "rb");
    if (offset_cache_handle) {
        fseek(offset_cache_handle, 0, SEEK_END);
        int64 byte_size = ftell(offset_cache_handle);
        offsets = (int64*)malloc(byte_size);
		ASSERT(offsets != 0);
        num_frames = (int32)(byte_size / sizeof(int64));
        fread(offsets, sizeof(int64), num_frames, offset_cache_handle);
        fclose(offset_cache_handle);
    } else {
        if (read_xtc_frame_offsets(filename.cstr(), &num_frames, &offsets) != exdrOK) {
            LOG_ERROR("Could not read frame offsets in trajectory");
			xdrfile_close(file_handle);
            return false;
        }
        FILE* write_cache_handle = fopen(cache_file.cstr(), "wb");
        if (write_cache_handle) {
            fwrite(offsets, sizeof(int64), num_frames, write_cache_handle);
            fclose(write_cache_handle);
        }
    }

    if (!offsets) {
        return false;
    }

    init_trajectory(traj, num_atoms, num_frames);

    traj->num_atoms = num_atoms;
    traj->num_frames = 0;
    traj->total_simulation_time = 0;
    traj->simulation_type = MoleculeTrajectory::NVT;
    traj->file.path = allocate_string(filename);
    traj->file.handle = file_handle;
	traj->file.tag = XTC_FILE_TAG;
    traj->frame_offsets = {offsets, num_frames};

    return true;
}

bool read_next_trajectory_frame(MoleculeTrajectory* traj) {
    ASSERT(traj);
    if (!traj->file.handle) return false;
    auto num_frames = traj->frame_offsets.count;
    if (traj->num_frames == num_frames) return false;

    // Next index to be loaded
    int i = traj->num_frames;

    int step;
    float precision;
    float time;
    float matrix[3][3];
    float* pos_buf = (float*)TMP_MALLOC(traj->num_atoms * 3 * sizeof(float));
    defer { TMP_FREE(pos_buf); };

    read_xtc((XDRFILE*)traj->file.handle, traj->num_atoms, &step, &time, matrix, (float(*)[3])pos_buf, &precision);

    TrajectoryFrame* frame = traj->frame_buffer.ptr + i;
    for (int j = 0; j < traj->num_atoms; j++) {
        frame->atom_position.x[j] = 10.f * pos_buf[j * 3 + 0];
        frame->atom_position.y[j] = 10.f * pos_buf[j * 3 + 1];
        frame->atom_position.z[j] = 10.f * pos_buf[j * 3 + 2];
    }
    frame->box = mat3(matrix[0][0], matrix[0][1], matrix[0][2], matrix[1][0], matrix[1][1], matrix[1][2], matrix[2][0], matrix[2][1], matrix[2][2]) * 10.f;

    traj->num_frames++;
    return true;
}

bool close_file_handle(MoleculeTrajectory* traj) {
    ASSERT(traj);
	if (traj->file.tag != xtc::XTC_FILE_TAG) {
		LOG_ERROR("Wrong file tag for closing file handle... Expected XTC_FILE_TAG");
		return false;
	}

    if (traj->file.handle) {
        xdrfile_close((XDRFILE*)traj->file.handle);
        traj->file.handle = nullptr;
        return true;
    }
    return false;
}

}  // namespace xtc