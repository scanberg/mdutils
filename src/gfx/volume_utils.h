#pragma once

#include <core/gl.h>
#include <core/types.h>
#include <core/string_types.h>
#include <core/volume.h>

namespace volume {

void initialize();
void shutdown();

void create_volume_texture(GLuint* texture, const ivec3& dim);
void free_volume_texture(GLuint texture);
void create_tf_texture(GLuint* texture, int* width, CString path_to_file);
void set_volume_texture_data(GLuint texture, ivec3 dim, void* data);
mat4 compute_model_to_world_matrix(const vec3& min_world_aabb, const vec3& max_world_aabb);
mat4 compute_texture_to_model_matrix(const ivec3& dim);

void save_volume_to_file(const Volume& volume, CString path_to_file);

/*
    Renders a volumetric texture using OpenGL.
    - volume_texture: An OpenGL 3D texture containing the data.
    - tf_texture:     An OpenGL 1D texture containing the transfer function
    - depth_texture: An OpenGL 2D texture containing the depth data in the frame (for stopping ray traversal).
    - texture_matrix: Matrix containing texture to model transformation of the volume.
    - model_matrix: Matrix containing model to world transformation of the volume, which is assumed to occupy a unit cube [0,1] in its model-space.
    - view_matrix: Matrix containing world to view transformation of the camera.
    - proj_matrix: Matrix containing view to clip transformation of the camera.
    - color: The color of the voxels within the volume
    - density_scale: global scaling of densities
    - alpha_scale:   global alpha scaling of the transfer function
*/
void render_volume_texture(GLuint volume_texture, GLuint tf_texture, GLuint depth_texture, const mat4& texture_matrix, const mat4& model_matrix, const mat4& view_matrix,
                           const mat4& proj_matrix, float density_scale = 1.0f, float alpha_scale = 1.0f);

}  // namespace volume
