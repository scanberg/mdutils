#include "xtc_utils.h"
#include <core/string_utils.h>
#include <core/log.h>
#include <core/file.h>

#include <stdio.h>
#include <xdrfile_xtc.h>

namespace xtc {
/*
bool read_trajectory_frames_from_file(Array<TrajectoryFrame> frames, i32 num_atoms, CStringView filename, i64 byte_offset) {
    // Not implemented...
    ASSERT(false);
    return false;
}
*/

static bool generate_cache(CStringView filename) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        LOG_ERROR("Could not open file '.*s", filename.length(), filename.beg());
        return false;
    }
    fseeki64(file, 0, SEEK_END);
    const u64 file_size = (u64)ftelli64(file);
    fclose(file);

    StringBuffer<512> zfilename = filename;  // Make sure it is zero terminated
    StringBuffer<512> cache_file = get_directory(filename);
    cache_file += "/";
    cache_file += get_file_without_extension(filename);
    cache_file += ".cache";

    i64* tmp_offsets;
    i32 num_atoms, num_frames;
    if (read_xtc_header(zfilename.cstr(), &num_atoms, &num_frames, &tmp_offsets) != exdrOK) {
        LOG_ERROR("Could not read frame offsets in trajectory");
        return false;
    }
    defer { free(tmp_offsets); };

    FrameBytes* frame_bytes = (FrameBytes*)TMP_MALLOC(num_frames * sizeof(FrameBytes));
    defer { TMP_FREE(frame_bytes); };

    for (i32 i = 0; i < num_frames - 1; i++) {
        frame_bytes[i].offset = tmp_offsets[i];
        frame_bytes[i].extent = tmp_offsets[i + 1] - tmp_offsets[i];
    }
    frame_bytes[num_frames - 1].offset = tmp_offsets[num_frames - 1];
    frame_bytes[num_frames - 1].extent = file_size - tmp_offsets[num_frames - 1];

    u64 UID = generate_UID(filename);
    return write_trajectory_cache(UID, frame_bytes, num_frames, cache_file);
}

bool read_trajectory_num_frames(i32* num_frames, CStringView filename) {
    u64 UID;
    if ((UID = generate_UID(filename)) != INVALID_UID) {
        StringBuffer<512> cache_file = get_directory(filename);
        cache_file += "/";
        cache_file += get_file_without_extension(filename);
        cache_file += ".cache";

        u64 c_UID;
        i64 c_num_frames;
        if (read_trajectory_cache_header(&c_UID, &c_num_frames, cache_file) && (UID == c_UID)) {
            // Cache is valid
            *num_frames = (i32)c_num_frames;
            return true;
        } else {
            // Cache is invalid
            // Regenerate data
            if (generate_cache(filename)) {
                if (read_trajectory_cache_header(&c_UID, &c_num_frames, cache_file) && (UID == c_UID)) {
                    *num_frames = (i32)c_num_frames;
                    return true;
                }
            }
        }
    }
    return false;
}

bool read_trajectory_frame_bytes(FrameBytes* frame_bytes, CStringView filename) {
    u64 UID;
    if ((UID = generate_UID(filename)) != INVALID_UID) {
        StringBuffer<512> cache_file = get_directory(filename);
        cache_file += "/";
        cache_file += get_file_without_extension(filename);
        cache_file += ".cache";

        u64 c_UID;
        i64 c_num_frames;
        if (read_trajectory_cache_header(&c_UID, &c_num_frames, cache_file) && (UID == c_UID)) {
            // Cache is valid
            return read_trajectory_cache(frame_bytes, cache_file);
        } else {
            // Cache is invalid
            // Regenerate data
            if (generate_cache(filename)) {
                if (read_trajectory_cache_header(&c_UID, &c_num_frames, cache_file) && (UID == c_UID)) {
                    return read_trajectory_cache(frame_bytes, cache_file);
                }
            }
        }
    }
    return false;
}

bool decompress_trajectory_frame(TrajectoryFrame* frame, i32 num_atoms, Array<u8> raw_data) {
    ASSERT(frame);
    XDRFILE* file = xdrfile_mem(raw_data.beg(), raw_data.size_in_bytes(), "r");
    defer { xdrfile_close(file); };

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