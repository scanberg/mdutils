#version 150 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_proj_mat;
uniform mat4 u_curr_view_to_prev_clip_mat;
uniform float u_exposure = 1.0;

in GS_FS {
    flat vec4 color;
    flat vec4 view_sphere;
    flat vec4 picking_color;
    smooth vec4 view_coord;
    flat vec4 view_velocity;
} in_frag;

#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif
layout(location = 0) out vec4 out_color_alpha;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_ss_vel;
layout(location = 3) out vec4 out_picking_color;

// https://aras-p.info/texts/CompactNormalStorage.html
vec4 encode_normal (vec3 n) {
    float p = sqrt(n.z*8+8);
    return vec4(n.xy/p + 0.5,0,0);
}

void main() {
    vec3 center = in_frag.view_sphere.xyz;
    float radius = in_frag.view_sphere.w;
    vec3 view_dir = -normalize(in_frag.view_coord.xyz);

    vec3 m = -center;
    vec3 d = -view_dir;
    float r = radius;
    float b = dot(m, d);
    float c = dot(m, m) - r*r;
    float discr = b*b - c;
    if (discr < 0.0) discard;
    float t = -b -sqrt(discr);

    vec3 view_coord = d * t;
    vec3 view_normal = (view_coord - center) / radius;
    vec4 clip_coord = u_proj_mat * vec4(view_coord, 1);

    vec3 prev_view_coord = view_coord - in_frag.view_velocity.xyz;
    vec4 prev_clip_coord = u_curr_view_to_prev_clip_mat * vec4(prev_view_coord, 1);

    vec2 curr_uv = clip_coord.xy / clip_coord.w * 0.5 + 0.5;
    vec2 prev_uv = prev_clip_coord.xy / prev_clip_coord.w * 0.5 + 0.5;
    vec2 ss_vel = curr_uv - prev_uv;

    //gl_FragDepth = (-u_proj_mat[2][2] - u_proj_mat[3][2] / view_coord.z) * 0.5 + 0.5;
    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;
    out_color_alpha = in_frag.color;
    out_normal = encode_normal(view_normal);
    out_ss_vel = vec4(ss_vel * 100, 0, 0);
    out_picking_color = in_frag.picking_color;
}