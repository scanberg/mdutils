#pragma once

#include <core/types.h>
#include <core/array_types.h>
#include <core/vector_types.h>
#include <core/hash.h>

struct MoleculeStructure;
struct MoleculeTrajectory;

namespace structure_tracking {

typedef uint32 ID;
inline ID get_id(CString str) { return hash::crc32(str); }

void initialize();
void shutdown();

void compute_tracking_data(ID structure_id, Array<const bool> atom_mask, const MoleculeStructure& mol, const MoleculeTrajectory& traj, int32 reference_frame_idx = 0);

void transform_coordinates_to_reference(float* RESTRICT x, float* RESTRICT y, float* RESTRICT z, ID structure_id);

}  // namespace structure_tracking
