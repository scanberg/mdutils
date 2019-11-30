#pragma once

#include <core/camera.h>
#include <core/vector_types.h>

/*
struct FpsControllerState {
    struct {
        bool forward_button = false;
        bool backward_button = false;
        bool left_button = false;
        bool right_button = false;
        vec2 mouse_vel_ndc = {0, 0};
        float delta_time = 0.f;
    } input;

    struct {
        float move_speed = 1.f;
        float rotation_speed = 1.f;
    } params;
};
*/

struct TrackballControllerState {
    struct {
        bool rotate_button = false;
        bool pan_button = false;
        bool dolly_button = false;
        float dolly_delta = 0;
        vec2 mouse_coord_prev = vec2(0, 0);
        vec2 mouse_coord_curr = vec2(0, 0);
        vec2 screen_size = vec2(0, 0);
    } input;

    struct {
        float pan_scale = 0.5f;
        float pan_exponent = 1.f;
        float dolly_drag_scale = 0.01f;
        float dolly_drag_exponent = 1.1f;
        float dolly_delta_scale = 0.1f;
        float dolly_delta_exponent = 1.1f;
        float min_distance = 1.f;
        float max_distance = 1000.f;
    } params;

    // Distance to the rotational focus point
    float distance = 14.f;
};

enum TrackballFlags_ { TrackballFlags_RotateReturnsTrue = 0x1, TrackballFlags_PanReturnsTrue = 0x2, TrackballFlags_DollyReturnsTrue = 0x4, TrackballFlags_AnyInteractionReturnsTrue = 0xFFFFFFFF };
typedef uint32_t TrackballFlags;

mat3 look_at(const vec3& look_from, const vec3& look_at, const vec3& look_up = {0, 1, 0});

void camera_trackball(Camera* camera, vec2 prev_ndc, vec2 curr_ndc);
void camera_move(Camera* camera, vec3 vec);

mat4 compute_world_to_view_matrix(const Camera& camera);
mat4 compute_view_to_world_matrix(const Camera& camera);
mat4 compute_perspective_projection_matrix(const Camera& camera, int width, int height);
mat4 compute_perspective_projection_matrix(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y);
mat4 compute_orthographic_projection_matrix(const Camera& camera, int width, int height);
mat4 compute_orthographic_projection_matrix(const Camera& camera, int width, int height, float texel_offset_x, float texel_offset_y);


bool camera_controller_trackball(vec3* position, quat* orientation, TrackballControllerState* trackball_state, TrackballFlags flags = 0xFFFFFFFF);

//void camera_controller_fps(Camera* camera, const FpsControllerState& state);