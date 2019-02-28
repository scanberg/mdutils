#pragma once

#include <core/types.h>
#include <core/gl.h>
#include <core/array_types.h>
#include <core/vector_types.h>
#include <mol/molecule_structure.h>
#include <mol/molecule_utils.h>
#include <gfx/view_param.h>

namespace draw {
void initialize();
void shutdown();

// These are just dummy structs to show the expected layout of bufferdata

// This the layout for control points that will be captured by the transform feedback into control_point_buffer and spline_buffer
// This is packed for 32-bit alignment because transform feedback only outputs full 32-bit types, so the real output is packed into uint32s.
struct ControlPoint {
    float position[3];
    int16 support_vector[3];
    int16 tangent_vector[3];
    uint8 classification[4];  // classification probabilities: Coil, Sheet, Helix etc.
    uint32 atom_index;
};

struct AtomPosition {
    float x, y, z;
};

struct AtomRadius {
    float radius;
};

struct AtomColor {
    uint8 r, g, b, a;
};

struct AtomVelocity {
    float x, y, z;
};

struct AtomMask {
    uint8 mask;
};

struct Bond {
    uint32 atom_idx[2];
};

void draw_vdw(GLuint atom_position_buffer, GLuint atom_radius_buffer, GLuint atom_color_buffer, GLuint atom_velocity_buffer, int32 atom_count, const ViewParam& view_param, float radius_scale = 1.f);
void draw_licorice(GLuint atom_position_buffer, GLuint atom_color_buffer, GLuint atom_velocity_buffer, GLuint bond_buffer, int32 bond_count, const ViewParam& view_param, float radius_scale = 1.f);
void draw_ribbons(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, GLuint atom_velocity_buffer, int32 num_spline_indices, const ViewParam& view_param);
void draw_cartoon(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, int32 num_spline_indices, const ViewParam& view_param);
void draw_spline(GLuint spline_buffer, GLuint spline_index_buffer, int32 num_spline_indices, const ViewParam& view_param, uint32 s_color = 0xFF00FF00, uint32 v_color = 0xFF0000FF,
                 uint32 t_color = 0xFFFF0000);

namespace lean_and_mean {
void draw_vdw(GLuint atom_position_buffer, GLuint atom_radius_buffer, GLuint atom_color_buffer, GLuint atom_mask_buffer, int32 atom_count, const ViewParam& view_param, float radius_scale = 1.f,
              vec4 color = vec4(1, 1, 1, 1), uint32 mask = 0xFFFFFFFFU);
void draw_licorice(GLuint atom_position_buffer, GLuint atom_color_buffer, GLuint atom_mask_buffer, GLuint bond_buffer, int32 bond_count, const ViewParam& view_param, float radius_scale = 1.f,
                   vec4 color = vec4(1, 1, 1, 1), uint32 mask = 0xFFFFFFFFU);
void draw_ribbons(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, GLuint atom_mask_buffer, int32 num_spline_indices, const ViewParam& view_param, float scale = 1.f,
                  vec4 color = vec4(1, 1, 1, 1), uint32 mask = 0xFFFFFFFFU);

// void draw_ribbons(GLuint spline_buffer, GLuint spline_index_buffer, int32 num_spline_indices, const ViewParam& view_param);
// void draw_cartoon(GLuint spline_buffer, GLuint spline_index_buffer, GLuint atom_color_buffer, int32 num_spline_indices, const ViewParam& view_param);
}  // namespace lean_and_mean

void compute_backbone_control_points(GLuint dst_buffer, GLuint atom_position_buffer, GLuint backbone_index_buffer, int num_backbone_indices, GLuint ramachandran_tex);
void compute_backbone_spline(GLuint dst_buffer, GLuint control_point_buffer, GLuint control_point_index_buffer, int num_control_point_indices, float tension = 0.5);

}  // namespace draw
