#pragma once

#include <core/types.h>
#include <core/array_types.h>
#include <core/vector_types.h>
#include <core/hash.h>

struct MoleculeDynamic;

struct Transform {
    mat3 rotation = {};
	vec3 com = {};
};

namespace structure_tracking {
typedef uint32 ID;

void initialize();
void shutdown();

ID create_structure();
bool remove_structure(ID structure_id);

void clear_structures();

mat3 compute_rotation(const float* RESTRICT x0, const float* RESTRICT y0, const float* RESTRICT z0,
					  const float* RESTRICT x1, const float* RESTRICT y1, const float* RESTRICT z1,
					  const float* RESTRICT mass, int64 count, const vec3& com0, const vec3& com1);

bool compute_trajectory_transform_data(ID structure_id, Array<const bool> atom_mask, const MoleculeDynamic& dynamic, int32 target_frame_idx = 0);

const Transform& get_transform_to_target_frame(ID structure_id, int32 source_frame);

}  // namespace structure_tracking
