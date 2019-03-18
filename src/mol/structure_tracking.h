#pragma once

#include <core/types.h>
#include <core/array_types.h>
#include <core/vector_types.h>
#include <core/hash.h>

struct MoleculeDynamic;

enum TransformFlag_ {
	TransformFlag_Translate = 1,
	TransformFlag_Rotate = 2,
	TransformFlag_All = 0xFFFFFFFF
};

typedef int TransformFlags;

struct Transform {
    mat3 rotation = {};
    vec3 translation = {};
};

namespace structure_tracking {

typedef uint32 ID;
inline ID get_id(CString str) { return hash::crc32(str); }

void initialize();
void shutdown();

bool create_structure(ID structure_id);
bool remove_structure(ID structure_id);

void clear_structures();

mat3 compute_rotation(const float* RESTRICT x0, const float* RESTRICT y0, const float* RESTRICT z0,
					  const float* RESTRICT x1, const float* RESTRICT y1, const float* RESTRICT z1,
					  const float* RESTRICT mass, int64 count, const vec3& com0, const vec3& com1);

bool compute_trajectory_transform_data(ID structure_id, Array<const bool> atom_mask, const MoleculeDynamic& dynamic, int32 target_frame_idx = 0);

const Transform& get_transform_to_target_frame(ID structure_id, int32 source_frame);

void apply_transform(float* RESTRICT x, float* RESTRICT y, float* RESTRICT z, int64 count, const Transform& t, TransformFlags flags = TransformFlag_All);

}  // namespace structure_tracking
