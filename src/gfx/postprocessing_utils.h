#pragma once

#include <core/gl.h>
#include <core/types.h>
#include <core/vector_types.h>

namespace postprocessing {

void initialize(int width, int height);
void shutdown();

typedef int Tonemapping;
enum Tonemapping_ { Tonemapping_Passthrough, Tonemapping_ExposureGamma, Tonemapping_Filmic };

void shade_deferred(GLuint depth_tex, GLuint color_tex, GLuint normal_tex, const mat4& inv_proj_matrix);
void highlight_selection(GLuint atom_idx_tex, GLuint selection_buffer);
void apply_dof(GLuint depth_tex, GLuint color_tex, const mat4& proj_matrix, float focus_point, float focus_scale);
void apply_ssao(GLuint depth_tex, GLuint normal_tex, const mat4& proj_mat, float intensity = 1.5f, float radius = 3.f, float bias = 0.1f);
void apply_tonemapping(GLuint color_tex, Tonemapping tonemapping, float exposure = 1.0f, float gamma = 2.2f);
void blit_texture(GLuint tex);

}  // namespace postprocessing
