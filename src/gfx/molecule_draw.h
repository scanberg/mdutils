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

// This the layout for control points that will be captured by the transform feedback into control_point_buffer and spline_buffer
// This is packed for 32-bit alignment, but transform feedback only outputs full 32-bit types, so the real output is packed into uint32s.
struct ControlPoint {
    float control_point[3];
    short support_vector[3];
    short tangent_vector[3];
    short backbone_angles[2];
    uint32 atom_index;
};

constexpr int cp_size = sizeof(ControlPoint);

void draw_vdw(GLuint atom_position_radius_buffer, GLuint atom_color_buffer, int32 atom_count, const mat4& view_mat, const mat4& proj_mat, float radius_scale = 1.f);
void draw_licorice(GLuint atom_position_buffer, GLuint atom_color_buffer, GLuint bond_buffer, int32 bond_count, const mat4& view_mat, const mat4& proj_mat, float radius_scale = 1.f);
void draw_ribbons(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, int32 num_spline_indices, const mat4& view_mat, const mat4& proj_mat);
void draw_cartoon(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, int32 num_spline_indices, GLuint ramachandran_tex, const mat4& view_mat, const mat4& proj_mat);
void draw_spline(GLuint spline_buffer, GLuint spline_index_buffer, int32 num_spline_indices, const mat4& view_proj_mat, uint32 s_color = 0xFF00FF00, uint32 v_color = 0xFF0000FF,
                 uint32 t_color = 0xFFFF0000);

void compute_backbone_control_points(GLuint dst_buffer, GLuint atom_position_buffer, GLuint backbone_index_buffer, int num_backbone_indices);
void compute_backbone_spline(GLuint dst_buffer, GLuint control_point_buffer, GLuint control_point_index_buffer, int num_control_point_indices, int num_subdivisions = 8, float tension = 0.5);

}  // namespace draw
