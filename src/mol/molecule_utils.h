#pragma once

#include <core/bitfield.h>
#include <core/array_types.h>
#include <core/math_utils.h>
#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/molecule_dynamic.h>
#include <mol/aminoacid.h>

struct AABB {
    vec3 min = {};
    vec3 max = {};
    vec3 ext() const { return max - min; }
};

struct EigenFrame {
    vec3 vectors[3];
    float values[3];
};

// Translates a set of points
void translate(soa_vec3 in_out, i64 count, const vec3& translation);
void translate(soa_vec3 out, const soa_vec3 in, i64 count, const vec3& translation);

// Transforms points as homogeneous vectors[x,y,z,w*] with supplied transformation matrix (NO 'perspective' division is done)
// W-component is supplied by user
void transform(soa_vec3 in_out, i64 count, const mat4& transformation, float w_comp = 1.0f);
void transform(soa_vec3 out, const soa_vec3 in, i64 count, const mat4& transformation, float w_comp = 1.0f);

// Transforms points as homogeneous vectors[x,y,z,w*] with supplied transformation matrix and applies 'perspective' division
// W-component is supplied by user
void homogeneous_transform(soa_vec3 in_out, i64 count, const mat4& transformation, float w_comp = 1.0f);
void homogeneous_transform(soa_vec3 out, const soa_vec3 in, i64 count, const mat4& transformation, float w_comp = 1.0f);

// Computes the minimun spanning Axis aligned bounding box which contains all supplied points [x,y,z,(r)adius)]
AABB compute_aabb(const soa_vec3 in_pos, i64 count);
AABB compute_aabb(const soa_vec3 in_pos, const float in_radius[], i64 count);

vec3 compute_com(const float in_x[], const float in_y[], const float in_z[], i64 count);
inline vec3 compute_com(const soa_vec3& in_pos, i64 count) { return compute_com(in_pos.x, in_pos.y, in_pos.z, count); }

vec3 compute_com(const float in_x[], const float in_y[], const float in_z[], const float in_mass[], i64 count);
inline vec3 compute_com(const soa_vec3& in_pos, const float in_mass[], i64 count) { return compute_com(in_pos.x, in_pos.y, in_pos.z, in_mass, count); }

vec3 compute_com_periodic(const soa_vec3 in_position, const float in_mass[], i64 count, const mat3& box);
vec3 compute_com_periodic_ref(const soa_vec3 in_position, const float in_mass[], i64 count, const mat3& box);

mat3 compute_covariance_matrix(const soa_vec3 in_position, const float in_mass[], i64 count, const vec3& com);

EigenFrame compute_eigen_frame(const soa_vec3 in_position, const float in_mass[], i64 count);

void linear_interpolation(soa_vec3 out_position, const soa_vec3 in_pos[2], i64 count, float t);
void linear_interpolation_pbc(soa_vec3 out_position, const soa_vec3 in_pos[2], i64 count, float t, const mat3& sim_box);

void cubic_interpolation(soa_vec3 out_position, const soa_vec3 in_pos[4], i64 count, float t);
void cubic_interpolation_pbc(soa_vec3 out_position, const soa_vec3 in_pos[4], i64 count, float t, const mat3& sim_box);

inline vec3 apply_pbc(const vec3& pos, const mat3& sim_box) {
    const vec3 ext = sim_box * vec3(1.0f);
    vec3 p = math::fract(pos / ext);
    if (p.x < 0.0f) p.x += 1.0f;
    if (p.y < 0.0f) p.y += 1.0f;
    if (p.z < 0.0f) p.z += 1.0f;
    return p * ext;
}

inline vec3 apply_pbc(const vec3& pos) {
    vec3 p = math::fract(pos);
    if (p.x < 0.0f) p.x += 1.0f;
    if (p.y < 0.0f) p.y += 1.0f;
    if (p.z < 0.0f) p.z += 1.0f;
    return p;
}

// Makes sure that all atomic positions fall inside of the simulation box, otherwise they are periodically transformed to end up within the box
void apply_pbc(soa_vec3 in_out_position, i64 count, const mat3& sim_box);

// Ranges correspond to e.g. either chains or residues which shall remain 'uncut' across the periodic boundary.
void apply_pbc(soa_vec3 out_position, const AtomRange range[], i64 num_ranges, const mat3& sim_box);

// Recenters a trajectory given a range of atoms which define a reference structure to compute a center of mass from
void recenter_trajectory(MoleculeDynamic* dynamic, AtomRange range);

inline bool valid_segment(const BackboneSegment& seg) {
    return (seg.ca_idx != seg.c_idx) && (seg.ca_idx != seg.n_idx) && (seg.ca_idx != seg.o_idx)
        && (seg.c_idx != seg.n_idx) && (seg.c_idx != seg.o_idx) && (seg.n_idx != seg.o_idx);
}

// Computes the dihedral angles within the backbone:
// phi   = dihedral( C[i-1], N[i],  CA[i],  C[i])
// psi   = dihedral( N[i],  CA[i],   C[i],  N[i+1])
// As explained here https://en.wikipedia.org/wiki/Ramachandran_plot.
void compute_backbone_angles(BackboneAngle out_angle[], const soa_vec3 in_pos, const BackboneSegment in_segments[], i64 num_segments);

void compute_atom_radius(float out_radii[], const Element in_element[], i64 count);
void compute_atom_mass(float out_mass[], const Element in_element[], i64 count);

bool is_amino_acid(const Label& residue_label);
bool is_dna(const Label& residue_label);

DynamicArray<Label> get_unique_residue_types(const MoleculeStructure& mol);
DynamicArray<ResIdx> get_residues_by_name(const MoleculeStructure& mol, CStringView name);

DynamicArray<AtomRange> find_equivalent_structures(const MoleculeStructure& mol, AtomRange ref);
DynamicArray<int> find_equivalent_structures(const MoleculeStructure& mol, Bitfield ref_mask, int ref_offset);