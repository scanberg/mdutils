#pragma once

#include <core/bitfield.h>
#include <core/array_types.h>
#include <core/math_utils.h>
#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/molecule_dynamic.h>
#include <mol/aminoacid.h>

struct BackboneAnglesTrajectory {
    int num_segments = 0;
    int num_frames = 0;
    ArrayView<BackboneAngle> angle_data{};
};

struct AABB {
    vec3 min = {};
    vec3 max = {};
};

struct EigenFrame {
    mat3 vectors;
    vec3 values;
};

inline ArrayView<BackboneAngle> get_backbone_angles(BackboneAnglesTrajectory& backbone_angle_traj, int frame_index) {
    if (backbone_angle_traj.angle_data.count == 0 || backbone_angle_traj.num_segments == 0) return {};
    ASSERT(frame_index < backbone_angle_traj.angle_data.count / backbone_angle_traj.num_segments);
    return ArrayView<BackboneAngle>(&backbone_angle_traj.angle_data[frame_index * backbone_angle_traj.num_segments], backbone_angle_traj.num_segments);
}

inline ArrayView<BackboneAngle> get_backbone_angles(BackboneAnglesTrajectory& backbone_angle_traj, int frame_offset, int frame_count) {
    if (backbone_angle_traj.angle_data.count == 0 || backbone_angle_traj.num_segments == 0) return {};
#ifdef DEBUG
    int32 num_frames = (int32)backbone_angle_traj.angle_data.count / backbone_angle_traj.num_segments;
    ASSERT(frame_offset < num_frames);
    ASSERT(frame_offset + frame_count <= num_frames);
#endif
    return backbone_angle_traj.angle_data.subarray(frame_offset * backbone_angle_traj.num_segments, frame_count * backbone_angle_traj.num_segments);
}

inline int32 get_backbone_angles_trajectory_current_frame_count(const BackboneAnglesTrajectory& backbone_angle_traj) {
    if (backbone_angle_traj.angle_data.count == 0 || backbone_angle_traj.num_segments == 0) return 0;
    return (int32)backbone_angle_traj.angle_data.count / backbone_angle_traj.num_segments;
}

inline ArrayView<BackboneAngle> get_backbone_angles(BackboneAnglesTrajectory& backbone_angle_traj, int frame_index, Chain chain) {
    return get_backbone_angles(backbone_angle_traj, frame_index).subarray(chain.res_range);
}

DynamicArray<BackboneSequence> compute_backbone_sequences(ArrayView<const BackboneSegment> segments, ArrayView<const Residue> residues);
DynamicArray<BackboneSegment> compute_backbone_segments(ArrayView<const Residue> residues, ArrayView<const Label> atom_labels);
// DynamicArray<SplineSegment> compute_spline(Array<const vec3> atom_pos, Array<const uint32> colors, Array<const BackboneSegment> backbone, int32 num_subdivisions = 1, float tension = 0.5f);

// Computes the dihedral angles within the backbone:
// phi   = dihedral( C[i-1], N[i],  CA[i],  C[i])
// psi   = dihedral( N[i],  CA[i],   C[i],  N[i+1])
// As seen here https://en.wikipedia.org/wiki/Ramachandran_plot.
DynamicArray<BackboneAngle> compute_backbone_angles(ArrayView<const BackboneSegment> backbone_segments, const float* pos_x, const float* pos_y, const float* pos_z);
void compute_backbone_angles(ArrayView<BackboneAngle> dst, ArrayView<const BackboneSegment> backbone_segments, const float* pos_x, const float* pos_y, const float* pos_z);

DynamicArray<BackboneAngle> compute_backbone_angles(ArrayView<const BackboneSegment> segments, ArrayView<const BackboneSequence> sequences, const float* pos_x, const float* pos_y, const float* pos_z);
void compute_backbone_angles(ArrayView<BackboneAngle> dst, ArrayView<const BackboneSegment> segments, ArrayView<const BackboneSequence> sequences, const float* pos_x, const float* pos_y, const float* pos_z);

void init_backbone_angles_trajectory(BackboneAnglesTrajectory* data, const MoleculeDynamic& dynamic);
void free_backbone_angles_trajectory(BackboneAnglesTrajectory* data);
void compute_backbone_angles_trajectory(BackboneAnglesTrajectory* bb_angle_traj, const MoleculeDynamic& dynamic);

void translate(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, int64 count, const vec3& translation);

// Transforms points as homogeneous vectors[x,y,z,w*] with supplied transformation matrix (NO perspective division is done)
// W-component is supplied by user
void transform_ref(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, int64 count, const mat4& transformation, float w_comp = 1.0f);

// Transforms points as homogeneous vectors[x,y,z,w*] with supplied transformation matrix (NO perspective division is done)
// W-component is supplied by user
void transform(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, int64 count, const mat4& transformation, float w_comp = 1.0f);
void transform(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z, float* RESTRICT in_x, float* RESTRICT in_y, float* RESTRICT in_z, int64 count, const mat4& transformation,
               float w_comp = 1.0f);

// Transforms points as homogeneous vectors[x,y,z,1] with supplied transformation matrix and applies 'perspective' division
void homogeneous_transform(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, int64 count, const mat4& transformation);

// Computes the minimun spanning Axis aligned bounding box which contains all supplied points [x,y,z,(r)adius)]
// Writes the results to out variables as float[3] = {x,y,z}
AABB compute_aabb(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, int64 count);
AABB compute_aabb(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, const float* in_r, int64 count);

vec3 compute_com(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, int64 count);
vec3 compute_com(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, const float* RESTRICT in_mass, int64 count);
vec3 compute_com(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, const Element* RESTRICT element, int64 count);

vec3 compute_com_periodic(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, const float* RESTRICT in_mass, int64 count, const mat3& box);
vec3 compute_com_periodic_ref(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, const float* RESTRICT in_mass, int64 count, const mat3& box);

mat3 compute_covariance_matrix(const float* RESTRICT x, const float* RESTRICT y, const float* RESTRICT z, const float* RESTRICT mass, int64 count, const vec3& com);

EigenFrame compute_eigen_frame(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, const float* RESTRICT in_mass, int64 count);

// clang-format off
void linear_interpolation(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
						  const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
						  const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
						  int64 count, float t);

void linear_interpolation_pbc(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
							  const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
							  const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
							  int64 count, float t, const mat3& sim_box);

void cubic_interpolation(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
						 const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
						 const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
						 const float* RESTRICT in_x2, const float* RESTRICT in_y2, const float* RESTRICT in_z2,
						 const float* RESTRICT in_x3, const float* RESTRICT in_y3, const float* RESTRICT in_z3,
						 int64 count, float t);

void cubic_interpolation_pbc(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
							 const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
							 const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
                             const float* RESTRICT in_x2, const float* RESTRICT in_y2, const float* RESTRICT in_z2,
							 const float* RESTRICT in_x3, const float* RESTRICT in_y3, const float* RESTRICT in_z3,
							 int64 count, float t, const mat3& sim_box);
// clang-format on

inline vec3 apply_pbc(const vec3& pos, const mat3& sim_box) {
    const vec3 ext = sim_box * vec3(1.0f);
    return math::fract(pos / ext) * ext;
}

void apply_pbc(float* RESTRICT x, float* RESTRICT y, float* RESTRICT z, const float* RESTRICT mass, ArrayView<const Sequence> sequences, const mat3& sim_box);

//void recenter_trajectory_on_residue(MoleculeDynamic* dynamic, ResIdx target_residue);

void recenter_trajectory(MoleculeDynamic* dynamic, Bitfield atom_mask);

// This computes heuristical covalent bonds in a hierarchical way (first internal, then external per residue) and stores the indices to the bonds
// within the residues. Only adjacent residues can form external covalent bonds in this function.
DynamicArray<Bond> compute_covalent_bonds(ArrayView<Residue> residues, const float* pos_x, const float* pos_y, const float* pos_z, const Element* element, int64 count);

// This is computes heuristical covalent bonds between any atoms without hierarchical constraints.
DynamicArray<Bond> compute_covalent_bonds(const float* pos_x, const float* pos_y, const float* pos_z, const Element* element, int64 count);

bool has_covalent_bond(const Residue& res_a, const Residue& res_b);
bool valid_segment(const BackboneSegment& seg);

DynamicArray<Sequence> compute_sequences(ArrayView<const Residue> residue);

DynamicArray<float> compute_atom_radii(ArrayView<const Element> elements);
void compute_atom_radii(float* out_radii, const Element* element, int64 count);

DynamicArray<float> compute_atom_masses(ArrayView<const Element> elements);
void compute_atom_masses(float* out_mass, const Element* element, int64 count);

bool is_amino_acid(const Residue& res);
bool is_dna(const Residue& res);

DynamicArray<Label> get_unique_residue_types(const MoleculeStructure& mol);
DynamicArray<ResIdx> get_residues_by_name(const MoleculeStructure& mol, CStringView name);

DynamicArray<AtomRange> find_equivalent_structures(const MoleculeStructure& mol, AtomRange ref);
DynamicArray<int> find_equivalent_structures(const MoleculeStructure& mol, Bitfield ref_mask, int ref_offset);