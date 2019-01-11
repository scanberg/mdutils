#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_view_mat;
uniform mat4 u_proj_mat;
uniform mat4 u_inv_proj_mat;

uniform float u_radius_scale = 1.0;

layout (location = 0) in vec4 in_pos_rad;
layout (location = 1) in vec4 in_color;

out VS_GS {
    flat vec4 view_sphere;
    flat vec4 color;
    flat vec4 picking_color;
    flat vec2 axis_a;
    flat vec2 axis_b;
    flat vec2 center;
    flat float inv_aspect_ratio;
    flat float z;
} out_geom;

vec4 pack_u32(uint data) {
    return vec4(
        (data & uint(0x000000FF)) >> 0,
        (data & uint(0x0000FF00)) >> 8,
        (data & uint(0x00FF0000)) >> 16,
        (data & uint(0xFF000000)) >> 24) / 255.0;
}

// From Inigo Quilez!
void proj_sphere(in vec4 sphere, 
                 in float fle,
                 out vec2 axis_a,
                 out vec2 axis_b,
                 out vec2 center) {
    vec3  o = sphere.xyz;
    float r2 = sphere.w*sphere.w;
    float z2 = o.z*o.z; 
    float l2 = dot(o,o);
    
    // axis
    axis_a = fle*sqrt(-r2*(r2-l2)/((l2-z2)*(r2-z2)*(r2-z2)))*vec2( o.x,o.y);
    axis_b = fle*sqrt(-r2*(r2-l2)/((l2-z2)*(r2-z2)*(r2-l2)))*vec2(-o.y,o.x);
    center = -fle*o.z*o.xy/(z2-r2);
}

void main() {
    vec4 pos_rad = in_pos_rad;
    vec4 color = in_color;
    vec4 picking_color = pack_u32(uint(gl_VertexID));

    vec3 pos = pos_rad.xyz;
    float rad = pos_rad.w * u_radius_scale;
    vec4 view_coord = u_view_mat * vec4(pos, 1.0);
    vec4 view_sphere = vec4(view_coord.xyz, rad);

    // Focal length
    float fle = u_proj_mat[1][1];
    // 1.0 / aspect_ratio
    float inv_ar = u_proj_mat[0][0] / u_proj_mat[1][1];
    // Compute view depth (with bias of radius) 
    float z = -u_proj_mat[2][2] - u_proj_mat[3][2] / (view_coord.z + rad);

    vec2 axis_a;
    vec2 axis_b;
    vec2 center;
    proj_sphere(view_sphere, fle, axis_a, axis_b, center);

    out_geom.view_sphere = view_sphere;
    out_geom.color = color;
    out_geom.picking_color = picking_color;
    out_geom.axis_a = axis_a;
    out_geom.axis_b = axis_b;
    out_geom.center = center;
    out_geom.inv_aspect_ratio = inv_ar;
    out_geom.z = z;

    gl_Position = view_coord;
}