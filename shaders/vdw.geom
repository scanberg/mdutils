#version 150 core

uniform mat4 u_inv_proj_mat;

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

in VS_GS {
    flat vec4 view_sphere;
    flat vec4 color;
    flat vec4 picking_color;
    flat vec2 axis_a;
    flat vec2 axis_b;
    flat vec2 center;
    flat float inv_aspect_ratio;
    flat float z;
} in_vert[];

out GS_FS {
    flat vec4 color;
    flat vec4 view_sphere;
    flat vec4 picking_color;
    smooth vec4 view_coord;
} out_frag;

void emit_vertex(vec2 uv) {
    vec2 axis_a = in_vert[0].axis_a;
    vec2 axis_b = in_vert[0].axis_b;
    vec2 center = in_vert[0].center;
    float inv_aspect_ratio = in_vert[0].inv_aspect_ratio;
    float z = in_vert[0].z;

    vec2 xy = (center + axis_a * uv.x + axis_b * uv.y) * vec2(inv_aspect_ratio, 1.0);
    vec4 pc = vec4(xy, z, 1);
    vec4 vc = u_inv_proj_mat * pc;

    out_frag.view_coord = vc / vc.w;
    gl_Position = pc;
    EmitVertex();
}

void main()
{
    if (in_vert[0].color.a == 0 || in_vert[0].view_sphere.w == 0) {
        EndPrimitive();
        return;
    }

    out_frag.color = in_vert[0].color;
    out_frag.view_sphere = in_vert[0].view_sphere;
    out_frag.picking_color = in_vert[0].picking_color;

    emit_vertex(vec2(-1,-1));
    emit_vertex(vec2( 1,-1));
    emit_vertex(vec2(-1, 1));
    emit_vertex(vec2( 1, 1));

    EndPrimitive();
}