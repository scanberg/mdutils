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

DynamicArray<BackboneSequence> compute_backbone_sequences(Array<const BackboneSegment> segments, Array<const Residue> residues);
DynamicArray<BackboneSegment> compute_backbone_segments(Array<const Residue> residues, Array<const Label> atom_labels);

// Computes the dihedral angles within the backbone:
// phi   = dihedral( C[i-1], N[i],  CA[i],  C[i])
// psi   = dihedral( N[i],  CA[i],   C[i],  N[i+1])
// As explained here https://en.wikipedia.org/wiki/Ramachandran_plot.

void compute_backbone_angles(BackboneAngle* dst, const BackboneSegment* backbone_segments, const soa_vec3 in_position, i64 num_segments);

void translate(soa_vec3 in_out, i64 count, const vec3& translation);

// Transforms points as homogeneous vectors[x,y,z,w*] with supplied transformation matrix (NO perspective division is done)
// W-component is supplied by user
void transform_ref(soa_vec3 in_out, i64 count, const mat4& transformation, float w_comp = 1.0f);

// Transforms points as homogeneous vectors[x,y,z,w*] with supplied transformation matrix (NO perspective division is done)
// W-component is supplied by user
void transform(soa_vec3 in_out, i64 count, const mat4& transformation, float w_comp = 1.0f);
void transform(soa_vec3 out, const soa_vec3 in, i64 count, const mat4& transformation, float w_comp = 1.0f);

// Transforms points as homogeneous vectors[x,y,z,1] with supplied transformation matrix and applies 'perspective' division
void homogeneous_transform(soa_vec3 in_out, i64 count, const mat4& transformation);

// Computes the minimun spanning Axis aligned bounding box which contains all supplied points [x,y,z,(r)adius)]
AABB compute_aabb(const soa_vec3 in_position, i64 count);
AABB compute_aabb(const soa_vec3 in_position, const float in_radius[], i64 count);

vec3 compute_com(const soa_vec3 in_position, i64 count);
vec3 compute_com(const soa_vec3 in_position, const float in_mass[], i64 count);
vec3 compute_com(const soa_vec3 in_position, const Element element[], i64 count);

vec3 compute_com_periodic(const soa_vec3 in_position, const float in_mass[], i64 count, const mat3& box);
vec3 compute_com_periodic_ref(const soa_vec3 in_position, const float in_mass[], i64 count, const mat3& box);

mat3 compute_covariance_matrix(const soa_vec3 in_position, const float in_mass[], i64 count, const vec3& com);

EigenFrame compute_eigen_frame(const soa_vec3 in_position, const float in_mass[], i64 count);

// clang-format off
void linear_interpolation(soa_vec3 out_position, const soa_vec3 in_position[2], i64 count, float t);
void linear_interpolation_pbc(soa_vec3 out_position, const soa_vec3 in_position[2], i64 count, float t, const mat3& sim_box);

void cubic_interpolation(soa_vec3 out_position, const soa_vec3 in_position[4], i64 count, float t);
void cubic_interpolation_pbc(soa_vec3 out_position, const soa_vec3 in_position[4], i64 count, float t, const mat3& sim_box);
// clang-format on

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
void apply_pbc(soa_vec3 in_out_position, const float in_mass[], i64 count, const mat3& sim_box);

// Ranges correspond to e.g. either chains or residues which shall remain 'uncut' across the periodic boundary.
void apply_pbc(soa_vec3 out_position, const float in_mass[], const AtomRange range[], i64 num_ranges, const mat3& sim_box);

// Recenters a trajectory given a mask which define a set of atoms to compute a center of mass from
void recenter_trajectory(MoleculeDynamic* dynamic, Bitfield atom_mask);

// This is computes heuristical covalent bonds in a hierarchical fashion
DynamicArray<Bond> compute_covalent_bonds(const MoleculeStructure& mol);

bool has_covalent_bond(const BondRange& res_bond_range_a, const BondRange& res_bond_range_b);
bool valid_segment(const BackboneSegment& seg);

//DynamicArray<ResRange> compute_sequences(Array<const Residue> residue);

//DynamicArray<float> compute_atom_radii(Array<const Element> elements);
void compute_atom_radii(float out_radii[], const Element in_element[], i64 count);

//DynamicArray<float> compute_atom_masses(Array<const Element> elements);
void compute_atom_masses(float out_mass[], const Element in_element[], i64 count);

bool is_amino_acid(const Label& residue_label);
bool is_dna(const Label& residue_label);

DynamicArray<Label> get_unique_residue_types(const MoleculeStructure& mol);
DynamicArray<ResIdx> get_residues_by_name(const MoleculeStructure& mol, CStringView name);

DynamicArray<AtomRange> find_equivalent_structures(const MoleculeStructure& mol, AtomRange ref);
DynamicArray<int> find_equivalent_structures(const MoleculeStructure& mol, Bitfield ref_mask, int ref_offset);