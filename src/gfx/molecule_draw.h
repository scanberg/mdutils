#pragma once

#include <core/types.h>
#include <core/gl.h>
#include <core/array_types.h>
#include <core/vector_types.h>
#include <mol/molecule_structure.h>
#include <mol/molecule_utils.h>

namespace draw {
void initialize();
void shutdown();

void draw_instanced_quads(int num_instances);

void draw_vdw(GLuint atom_position_radius_buffer, GLuint atom_color_buffer, int32 atom_count, const mat4& view_mat, const mat4& proj_mat,
              float radius_scale = 1.f);

void draw_licorice(GLuint atom_position_buffer, GLuint atom_color_buffer, GLuint bond_buffer, int32 bond_count, const mat4& view_mat,
                   const mat4& proj_mat, float radius_scale = 1.f);

void draw_ribbons(Array<const BackboneSegment> backbone_segments, Array<const Chain> chains, Array<const vec3> atom_positions,
                  Array<const uint32> atom_colors, const mat4& view_mat, const mat4& proj_mat, int num_subdivisions = 8, float tension = 0.5f,
                  float width_scale = 1.f, float thickness_scale = 1.f);

void compute_backbone_control_points(GLuint dst_buffer, GLuint atom_position_buffer, GLuint backbone_index_buffer, int num_backbone_indices);

void compute_backbone_spline(GLuint dst_buffer, GLuint control_point_buffer, GLuint control_point_index_buffer, int num_control_point_indices,
                             int num_subdivisions = 8, float tension = 0.5);

//void compute_backbone_spline2(GLuint dst_buffer, GLuint atom_position_buffer, GLuint backbone_index_buffer, int num_backbone_indices,
//                             int num_subdivisions = 8, float tension = 0.5);

// DEBUG
void draw_backbone(Array<const BackboneSegment> backbone, Array<const vec3> atom_positions, const mat4& view_mat, const mat4& proj_mat);
void draw_spline(Array<const SplineSegment> spline, const mat4& view_mat, const mat4& proj_mat);

void draw_spline(GLuint spline_buffer, GLuint spline_index_buffer, int32 num_spline_indices, const mat4& view_proj_mat, uint32 color = 0xFF00FF00);
void draw_support_vectors(GLuint spline_buffer, GLuint spline_index_buffer, int32 num_spline_indices, const mat4& view_proj_mat,
                          uint32 color = 0xFF00FF00);

void draw_ribbons(GLuint spline_buffer, GLuint spline_index_buffer, int32 num_spline_indices, const mat4& view_mat, const mat4& proj_mat, uint32 color);

}  // namespace draw
