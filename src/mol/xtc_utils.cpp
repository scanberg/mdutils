#include "xtc_utils.h"
#include <core/string_utils.h>
#include <core/log.h>
#include <core/file.h>

#include <stdio.h>
#include <xdrfile_xtc.h>

namespace xtc {
bool init_trajectory_from_file(MoleculeTrajectory* traj, i32 mol_atom_count, CStringView filename) {
    ASSERT(traj);
    free_trajectory(traj);

    CStringView directory = get_directory(filename);
    CStringView file = get_file_without_extension(filename);
    StringBuffer<512> cache_file = directory;
    cache_file += "/";
    cache_file += file;
    cache_file += ".cache";

    i64 xtc_file_size = 0;
    {
        FILE* f = fopen(filename, "rb");
        if (f) {
            fseeki64(f, 0, SEEK_END);
            xtc_file_size = ftelli64(f);
            fclose(f);
        }
    }

    XDRFILE* xtc_file_handle = xdrfile_open(filename.cstr(), "r");
    if (!xtc_file_handle) {
        LOG_ERROR("Could not open file: %s", filename);
        return false;
    }

    i32 num_atoms = 0;
    i32 num_frames = 0;
    i64* offsets = nullptr;

    if (read_xtc_natoms(filename.cstr(), &num_atoms) != exdrOK) {
        LOG_ERROR("Could not extract number of atoms in trajectory");
        xdrfile_close(xtc_file_handle);
        return false;
    }

    if (num_atoms != mol_atom_count) {
        LOG_ERROR("Trajectory atom count did not match molecule atom count");
        xdrfile_close(xtc_file_handle);
        return false;
    }

    // Try to read offsets from cache-file
    FILE* cache_handle = fopen(cache_file, "rb");
    if (cache_handle) {
        fseeki64(cache_handle, 0, SEEK_END);
        const i64 file_size = ftelli64(cache_handle);
        rewind(cache_handle);
        i64 expected_xtc_file_size;
        fread(&expected_xtc_file_size, sizeof(i64), 1, cache_handle);  // First i64 is the expected size of the xtc

        if (expected_xtc_file_size == xtc_file_size) {
            const i64 offset_size = file_size - sizeof(i64);  // minus expected size
            offsets = (i64*)malloc(offset_size);
            ASSERT(offsets);
            num_frames = (i32)(offset_size / sizeof(i64));
            fread(offsets, sizeof(i64), num_frames, cache_handle);
        }
        fclose(cache_handle);
    }

    // If no cache-file was found or was not fitting, then re-read all offsets and re-create offset-cache
    if (!offsets) {
        if (read_xtc_header(filename.cstr(), &num_atoms, &num_frames, &offsets) != exdrOK) {
            LOG_ERROR("Could not read frame offsets in trajectory");
            xdrfile_close(xtc_file_handle);
            return false;
        }
        FILE* write_cache_handle = fopen(cache_file, "wb");
        if (write_cache_handle) {
            fwrite(&xtc_file_size, sizeof(i64), 1, write_cache_handle);
            fwrite(offsets, sizeof(i64), num_frames, write_cache_handle);
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
    traj->simulation_type = SimulationType::NVT;
    traj->file.path = filename;
    traj->file.handle = xtc_file_handle;
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

    int natoms, step;
    float precision;
    float time;
    float matrix[3][3];
    float* pos_buf = (float*)TMP_MALLOC(traj->num_atoms * 3 * sizeof(float));
    defer { TMP_FREE(pos_buf); };

    read_xtc((XDRFILE*)traj->file.handle, &natoms, &step, &time, matrix, (float(*)[3])pos_buf, &precision);

    TrajectoryFrame* frame = traj->frame_buffer.ptr + i;
    for (int j = 0; j < traj->num_atoms; j++) {
        frame->atom_position.x[j] = 10.f * pos_buf[j * 3 + 0];
        frame->atom_position.y[j] = 10.f * pos_buf[j * 3 + 1];
        frame->atom_position.z[j] = 10.f * pos_buf[j * 3 + 2];
    }
    frame->box =
        mat3(matrix[0][0], matrix[0][1], matrix[0][2], matrix[1][0], matrix[1][1], matrix[1][2], matrix[2][0], matrix[2][1], matrix[2][2]) * 10.f;

    traj->num_frames++;
    return true;
}

bool read_trajectory_frame_from_memory(TrajectoryFrame* frame, i32 num_atoms, void* ptr, int64_t num_bytes) {
    ASSERT(frame);
    ASSERT(ptr);

    XDRFILE* file = xdrfile_mem(ptr, num_bytes, "r");
    defer { xdrfile_close(file); };
    if (!file) {
        return false;
    }

    int natoms, step;
    float precision;
    float time;
    float matrix[3][3];
    float* pos_buf = (float*)TMP_MALLOC(num_atoms * 3 * sizeof(float));
    defer { TMP_FREE(pos_buf); };

    read_xtc(file, &natoms, &step, &time, matrix, (float(*)[3])pos_buf, &precision);

    for (i32 i = 0; i < num_atoms; i++) {
        frame->atom_position.x[i] = 10.f * pos_buf[i * 3 + 0];
        frame->atom_position.y[i] = 10.f * pos_buf[i * 3 + 1];
        frame->atom_position.z[i] = 10.f * pos_buf[i * 3 + 2];
    }
    frame->box =
        mat3(matrix[0][0], matrix[0][1], matrix[0][2], matrix[1][0], matrix[1][1], matrix[1][2], matrix[2][0], matrix[2][1], matrix[2][2]) * 10.f;

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

DynamicArray<i64> read_frame_offsets(CStringView filename) {
    i64 file_size = 0;
    {
        FILE* f = fopen(filename, "rb");
        if (f) {
            fseeki64(f, 0, SEEK_END);
            file_size = ftelli64(f);
            fclose(f);
        }
    }

    if (file_size == 0) {
        LOG_ERROR("Could not locate file %.*s", filename.size(), filename.beg());
        return {};
    }

    StringBuffer<512> cache_file = get_directory(filename);
    cache_file += "/";
    cache_file += get_file_without_extension(filename);
    cache_file += ".cache";

    DynamicArray<i64> offsets = {};

    // Try to read offsets from cache-file
    FILE* cache_handle = fopen(cache_file, "rb");
    if (cache_handle) {
        fseeki64(cache_handle, 0, SEEK_END);
        const i64 cache_file_size = ftelli64(cache_handle);
        rewind(cache_handle);
        i64 expected_file_size;
        fread(&expected_file_size, sizeof(i64), 1, cache_handle);  // First entry is the expected size of the trajectory file

        if (expected_file_size == file_size) {
            const i64 num_offsets = cache_file_size / sizeof(i64) - 1;  // Minus first entry
            offsets.resize(num_offsets);
            fread(offsets.data(), sizeof(i64), num_offsets, cache_handle);
        }
        fclose(cache_handle);
    }

    if (offsets.empty()) {
        i64* tmp_offsets;
        i32 num_atoms, num_frames;
        if (read_xtc_header(filename.cstr(), &num_atoms, &num_frames, &tmp_offsets) != exdrOK) {
            LOG_ERROR("Could not read frame offsets in trajectory");
            return {};
        }

        FILE* write_cache_handle = fopen(cache_file, "wb");
        if (write_cache_handle) {
            fwrite(&file_size, sizeof(i64), 1, write_cache_handle);
            fwrite(offsets.data(), sizeof(i64), offsets.size(), write_cache_handle);
            fclose(write_cache_handle);
        }

        offsets.resize(num_frames);
        memcpy(offsets.data(), tmp_offsets, num_frames * sizeof(i64));
    }

    offsets.push_back(file_size);

    return offsets;
}

i32 read_num_frames(CStringView filename) { return (i32)read_frame_offsets(filename).size(); }

bool decompress_trajectory_frame(TrajectoryFrame* frame, i32 num_atoms, Array<u8> raw_data) {
    ASSERT(frame);
    XDRFILE* file = xdrfile_mem(raw_data.beg(), raw_data.size_in_bytes(), "r");
    if (!file) {
        LOG_ERROR("Could not open XDR-stream to memory location");
        return false;
    }

    float* pos_buf = (float*)TMP_MALLOC(num_atoms * 3 * sizeof(float));
    defer { TMP_FREE(pos_buf); };

    int natoms;
    float precision;
    read_xtc(file, &natoms, &frame->index, &frame->time, (float(&)[3][3])frame->box, (float(*)[3])pos_buf, &precision);

    // nm -> Ångström
    for (i32 i = 0; i < num_atoms; i++) {
        frame->atom_position.x[i] = 10.f * pos_buf[i * 3 + 0];
        frame->atom_position.y[i] = 10.f * pos_buf[i * 3 + 1];
        frame->atom_position.z[i] = 10.f * pos_buf[i * 3 + 2];
    }

    frame->box *= 10.0f;

    return true;
}

}  // namespace xtc