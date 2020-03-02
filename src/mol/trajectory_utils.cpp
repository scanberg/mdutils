#include "trajectory_utils.h"
#include <core/file.h>

u64 generate_UID(CStringView filename) {
    // @NOTE: Ideally, we want to use some type of hash generated from the entire file.
    // For now we just use the filesize mixed up with some versioning

    FILE* file = fopen(filename, "rb");
    if (file) {
        fseeki64(file, 0, SEEK_END);
        const u64 file_size = (u64)ftelli64(file);
        rewind(file);
        fclose(file);

        const u64 version = 3;
        return version ^ file_size;
    }

    return INVALID_UID;
}

// Reads the number of frames and unique ID of trajectory
bool read_trajectory_cache_header(u64* UID, i64* num_frames, CStringView cache_filename) {
    ASSERT(UID);
    ASSERT(num_frames);

    FILE* file = fopen(cache_filename, "rb");
    if (file) {
        fseeki64(file, 0, SEEK_END);
        const i64 file_size = ftelli64(file);
        rewind(file);
        fread(UID, sizeof(u64), 1, file);  // First entry is the UID
        fclose(file);
        *num_frames = file_size / sizeof(i64) - 1;  // Minus UID
        return true;
    }
    return false;
}

// Read trajectory Frame Byte Cache and the unique ID
bool read_trajectory_cache(FrameBytes* frame_bytes, CStringView cache_filename) {
    ASSERT(sizeof(FrameBytes) == 8);
    FILE* file = fopen(cache_filename, "rb");
    if (file) {
        fseeki64(file, 0, SEEK_END);
        const i64 file_size = ftelli64(file);
        rewind(file);
        const i64 num_frames = file_size / sizeof(i64) - 1;  // Minus UID
        fseeki64(file, sizeof(u64), SEEK_SET); // Skip UID
        fread(frame_bytes, sizeof(FrameBytes), num_frames, file);
        fclose(file);
        return true;
    }
    return false;
}

// Write trajectory Frame Byte Cache with a unique ID (fingerprint)
bool write_trajectory_cache(u64 UID, const FrameBytes* frame_bytes, i64 num_frames, CStringView cache_filename) {
    FILE* file = fopen(cache_filename, "wb");
    if (file) {
        fwrite(&UID, sizeof(u64), 1, file);
        fwrite(frame_bytes, sizeof(FrameBytes), num_frames, file);
        fclose(file);
        return true;
    }
    return false;
}
