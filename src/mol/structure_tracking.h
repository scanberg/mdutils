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

bool create_structure(ID structure_id);
bool remove_structure(ID structure_id);

void clear_structures();

bool compute_tracking_data(ID structure_id, Array<const bool> atom_mask, const MoleculeStructure& mol, const MoleculeTrajectory& traj, int32 reference_frame_idx = 0, float rbf_radial_cutoff = 10.f);
bool transform_coordinates_to_reference(float* RESTRICT x, float* RESTRICT y, float* RESTRICT z, ID structure_id, int32 frame_idx);

}  // namespace structure_tracking
