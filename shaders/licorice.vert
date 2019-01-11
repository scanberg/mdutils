#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_view_mat;

layout(location = 0) in vec3 v_position;
layout(location = 1) in vec4 v_color;

out Vertex {
    flat vec4 color;
    flat uint picking_id;
} out_vert;

void main() {
    gl_Position = u_view_mat * vec4(v_position, 1.0);
    out_vert.color = v_color;
    out_vert.picking_id = uint(gl_VertexID);
}