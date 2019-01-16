#version 330 core
uniform samplerBuffer u_atom_color_tex;
uniform sampler2D u_ramachandran_tex;

layout (location = 0) in vec3 v_control_point;
layout (location = 1) in vec3 v_support_vector;
layout (location = 2) in vec3 v_support_tangent;
layout (location = 3) in vec2 v_backbone_angles;
layout (location = 4) in uint v_atom_index;

out Vertex {
    vec4 control_point;
    vec4 support_vector;
    vec4 support_tangent;
    vec4 color;
    vec4 picking_color;
    vec3 weights;
} out_vert;

vec4 pack_u32(uint data) {
    return vec4(
        (data & uint(0x000000FF)) >> 0,
        (data & uint(0x0000FF00)) >> 8,
        (data & uint(0x00FF0000)) >> 16,
        (data & uint(0xFF000000)) >> 24) / 255.0;
}

void main() {
    // [-1, 1] -> [0,1]
    vec2 tc = vec2(0,1) + vec2(1,-1)*(v_backbone_angles * 0.5 + 0.5);
    vec3 rc = texture(u_ramachandran_tex, tc).rgb;

    out_vert.control_point = vec4(v_control_point, 1);
    out_vert.support_vector = vec4(v_support_vector, 0);
    out_vert.support_tangent = vec4(v_support_tangent, 0);
    out_vert.color = texelFetch(u_atom_color_tex, int(v_atom_index));
    out_vert.color = vec4(rc, 1);
    out_vert.picking_color = pack_u32(v_atom_index);
    float sheet_w = rc.b;
    float helix_w = rc.r;
    float tube_w = clamp(1.0 - sheet_w - helix_w, 0, 1);
    out_vert.weights = vec3(tube_w, helix_w, sheet_w);
    //out_vert.weights = vec3(0.5,0.5,0);
}