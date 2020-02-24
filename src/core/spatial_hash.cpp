#include "spatial_hash.h"
#include <core/common.h>
#include <core/math_utils.h>

#include <atomic>

namespace spatialhash {

void compute_frame(Frame* frame, const float* pos_x, const float* pos_y, const float* pos_z, i64 count, const vec3& cell_ext) {
    ASSERT(frame);
    if (count == 0) {
        *frame = {};
        return;
    };

    vec3 min_box(FLT_MAX);
    vec3 max_box(-FLT_MAX);
    for (i64 i = 0; i < count; i++) {
        const vec3 p = {pos_x[i], pos_y[i], pos_z[i]};
        min_box = math::min(min_box, p);
        max_box = math::max(max_box, p);
    }
    min_box -= 1.f;
    max_box += 1.f;
    compute_frame(frame, pos_x, pos_y, pos_z, count, cell_ext, min_box, max_box);
}

Frame compute_frame(const float* pos_x, const float* pos_y, const float* pos_z, i64 count, const vec3& cell_ext) {
    Frame frame;
    compute_frame(&frame, pos_x, pos_y, pos_z, count, cell_ext);
    return frame;
}

Frame compute_frame(const float* pos_x, const float* pos_y, const float* pos_z, i64 count, const vec3& cell_ext, const vec3& min_box, const vec3& max_box) {
    Frame frame;
    compute_frame(&frame, pos_x, pos_y, pos_z, count, cell_ext, min_box, max_box);
    return frame;
}

void compute_frame(Frame* frame, const float* pos_x, const float* pos_y, const float* pos_z, i64 count, const vec3& cell_ext, const vec3& min_box, const vec3& max_box) {
    ASSERT(frame);
    if (count == 0) return;

    frame->min_box = min_box;
    frame->max_box = max_box;
    frame->cell_count = math::max(ivec3(1), ivec3((max_box - min_box) / cell_ext));
    frame->cell_ext = (max_box - min_box) / (vec3)frame->cell_count;
    frame->cells.resize(frame->cell_count.x * frame->cell_count.y * frame->cell_count.z);
    frame->entries.resize(count);
    memset(frame->cells.data(), 0, frame->cells.size() * sizeof(Cell));

    u32* l_idx = (u32*)TMP_MALLOC(count * sizeof(u32));
    u32* g_idx = (u32*)TMP_MALLOC(count * sizeof(u32));
    defer {
        TMP_FREE(l_idx);
        TMP_FREE(g_idx);
    };

    for (int i = 0; i < count; i++) {
        const vec3 p = {pos_x[i], pos_y[i], pos_z[i]};
        int cell_idx = compute_cell_idx(*frame, p);
        l_idx[i] = frame->cells[cell_idx].count++;
        g_idx[i] = cell_idx;
    }

    for (int i = 1; i < frame->cells.size(); i++) {
        frame->cells[i].offset = frame->cells[i - 1].offset + frame->cells[i - 1].count;
    }

    for (int i = 0; i < frame->entries.size(); i++) {
        int dst = frame->cells[g_idx[i]].offset + l_idx[i];
        frame->entries[dst].position = {pos_x[i], pos_y[i], pos_z[i]};
        frame->entries[dst].index = i;
    }
}

}  // namespace spatialhash
