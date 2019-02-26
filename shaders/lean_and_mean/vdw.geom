#version 150 core

uniform mat4 u_inv_proj_mat;

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

in VS_GS {
    flat vec2 axis_a;
    flat vec2 axis_b;
    flat vec2 center;
    flat float z;
    flat float radius;
} in_vert[];

out GS_FS {
    smooth vec2 uv;
} out_frag;

void emit_vertex(vec2 uv) {
    vec2 axis_a = in_vert[0].axis_a;
    vec2 axis_b = in_vert[0].axis_b;
    vec2 center = in_vert[0].center;
    float z = in_vert[0].z;

    vec2 xy = (center + axis_a * uv.x + axis_b * uv.y);
    vec4 pc = vec4(xy, z, 1);

    out_frag.uv = uv;
    gl_Position = pc;
    EmitVertex();
}

void main()
{
    if (in_vert[0].radius == 0) {
        EndPrimitive();
        return;
    }

    emit_vertex(vec2(-1,-1));
    emit_vertex(vec2( 1,-1));
    emit_vertex(vec2(-1, 1));
    emit_vertex(vec2( 1, 1));

    EndPrimitive();
}