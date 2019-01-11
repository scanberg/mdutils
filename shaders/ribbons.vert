#version 330 core
uniform samplerBuffer u_atom_colors;

layout (location = 0) in vec3 v_control_point;
layout (location = 1) in vec3 v_support_vector;
layout (location = 2) in vec3 v_support_tangent;
layout (location = 3) in uint v_atom_index;

out Vertex {
    vec4 control_point;
    vec4 support_vector;
    vec4 support_tangent;
    vec4 color;
    vec4 picking_color;
    uint cap;
} out_vert;

vec4 pack_u32(uint data) {
    return vec4(
        (data & uint(0x000000FF)) >> 0,
        (data & uint(0x0000FF00)) >> 8,
        (data & uint(0x00FF0000)) >> 16,
        (data & uint(0xFF000000)) >> 24) / 255.0;
}

void main() {
    out_vert.control_point = vec4(v_control_point, 1);
    out_vert.support_vector = vec4(v_support_vector, 0);
    out_vert.support_tangent = vec4(v_support_tangent, 0);
    out_vert.color = texelFetch(u_atom_colors, int(v_atom_index));
    out_vert.picking_color = pack_u32(v_atom_index & uint(0x7FFFFFFF));
    //out_vert.cap = v_atom_index & uint(0x80000000);
    out_vert.cap = 1U;
}