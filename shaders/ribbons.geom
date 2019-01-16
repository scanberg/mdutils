#version 330 core

/*
    This shader can be optimized alot furhter.
    Especially for Intel integrated GPUs, where we want to minimize the workload in the geometry shader
    Move matrix multiplications and transformations to vertex shader
    Output 'Ready' corner points for the box that is the body for the ribbon-segment
    Output 'Ready' normals for the edges of the boxes
    In geometry shader: Pick the right indices and construct geometry
*/

uniform mat4 u_normal_mat;
uniform mat4 u_view_proj_mat;
uniform vec2 u_scale = vec2(1.0, 0.1);

layout(lines) in;
layout(triangle_strip, max_vertices = 24) out;

in Vertex {
    vec4 control_point;
    vec4 support_vector;
    vec4 support_tangent;
    vec4 color;
    vec4 picking_color;
} in_vert[];

out Fragment {
    smooth vec4 color;
    smooth vec4 view_normal;
    flat vec4 picking_color;
} out_frag;

void main() {
    vec4 p[2];
    vec4 x[2];
    vec4 y[2];
    vec4 z[2];
    
    p[0] = in_vert[0].control_point;
    p[1] = in_vert[1].control_point;
    x[0] = in_vert[0].support_vector;
    x[1] = in_vert[1].support_vector * sign(dot(in_vert[0].support_vector, in_vert[1].support_vector));
    z[0] = in_vert[0].support_tangent;
    z[1] = in_vert[1].support_tangent;
    y[0] = vec4(cross(z[0].xyz, x[0].xyz), 0); // To maintain right-handedness
    y[1] = vec4(cross(z[1].xyz, x[1].xyz), 0);

#if 0
    // This is to possible fix tesselation direction so it follows the curve of the segment
    float flip_sign = sign(dot(x[1].xyz, y[0].xyz));
    x[0].xyz *= flip_sign;
    x[1].xyz *= flip_sign;
    y[0].xyz *= flip_sign;
    y[1].xyz *= flip_sign;
    z[0].xyz *= flip_sign;
    z[1].xyz *= flip_sign;
#endif

    mat4 M[2];
    M[0] = mat4(x[0], y[0], z[0], p[0]);
    M[1] = mat4(x[1], y[1], z[1], p[1]);

    mat4 N[2];
    N[0] = u_normal_mat * M[0];
    N[1] = u_normal_mat * M[1];

    vec4 v0[4];
    v0[0] = u_view_proj_mat * M[0] * vec4(vec2(-1,-1) * u_scale, 0, 1);
    v0[1] = u_view_proj_mat * M[0] * vec4(vec2( 1,-1) * u_scale, 0, 1);
    v0[2] = u_view_proj_mat * M[0] * vec4(vec2(-1, 1) * u_scale, 0, 1);
    v0[3] = u_view_proj_mat * M[0] * vec4(vec2( 1, 1) * u_scale, 0, 1);

    vec4 v1[4];
    v1[0] = u_view_proj_mat * M[1] * vec4(vec2(-1,-1) * u_scale, 0, 1);
    v1[1] = u_view_proj_mat * M[1] * vec4(vec2( 1,-1) * u_scale, 0, 1);
    v1[2] = u_view_proj_mat * M[1] * vec4(vec2(-1, 1) * u_scale, 0, 1);
    v1[3] = u_view_proj_mat * M[1] * vec4(vec2( 1, 1) * u_scale, 0, 1);
 
    // BOTTOM
    out_frag.color = in_vert[0].color;
    out_frag.picking_color = in_vert[0].picking_color;
    out_frag.view_normal = N[0] * vec4( 0, -1, 0, 0);
    gl_Position = v0[0]; EmitVertex();
    gl_Position = v0[1]; EmitVertex();
    
    out_frag.color = in_vert[1].color;
    out_frag.picking_color = in_vert[1].picking_color;
    out_frag.view_normal = N[1] * vec4( 0, -1, 0, 0);
    gl_Position = v1[0]; EmitVertex();
    gl_Position = v1[1]; EmitVertex();
    EndPrimitive();

    // TOP
    out_frag.color = in_vert[0].color;
    out_frag.picking_color = in_vert[0].picking_color;
    out_frag.view_normal = N[0] * vec4( 0, 1, 0, 0);
    gl_Position = v0[3]; EmitVertex();
    gl_Position = v0[2]; EmitVertex();

    out_frag.color = in_vert[1].color;
    out_frag.picking_color = in_vert[1].picking_color;
    out_frag.view_normal = N[1] * vec4( 0, 1, 0, 0);
    gl_Position = v1[3]; EmitVertex();
    gl_Position = v1[2]; EmitVertex();
    EndPrimitive();

    // LEFT
    out_frag.color = in_vert[0].color;
    out_frag.picking_color = in_vert[0].picking_color;
    out_frag.view_normal = N[0] * vec4(-1, 0, 0, 0);
    gl_Position = v0[2]; EmitVertex();
    gl_Position = v0[0]; EmitVertex();

    out_frag.color = in_vert[1].color;
    out_frag.picking_color = in_vert[1].picking_color;
    out_frag.view_normal = N[1] * vec4(-1, 0, 0, 0);
    gl_Position = v1[2]; EmitVertex();
    gl_Position = v1[0]; EmitVertex();
    EndPrimitive();

    // RIGHT
    out_frag.color = in_vert[0].color;
    out_frag.picking_color = in_vert[0].picking_color;
    out_frag.view_normal = N[0] * vec4( 1, 0, 0, 0);
    gl_Position = v0[1]; EmitVertex();
    gl_Position = v0[3]; EmitVertex();

    out_frag.color = in_vert[1].color;
    out_frag.picking_color = in_vert[1].picking_color;
    out_frag.view_normal = N[1] * vec4( 1, 0, 0, 0);
    gl_Position = v1[1]; EmitVertex();
    gl_Position = v1[3]; EmitVertex();
    EndPrimitive();

    // FRONT
    out_frag.color = in_vert[1].color;
    out_frag.picking_color = in_vert[1].picking_color;
    out_frag.view_normal = N[1] * vec4( 0, 0, 1, 0);
    gl_Position = v1[0]; EmitVertex();
    gl_Position = v1[1]; EmitVertex();
    gl_Position = v1[2]; EmitVertex();
    gl_Position = v1[3]; EmitVertex();
    EndPrimitive();

    // BACK
    out_frag.color = in_vert[0].color;
    out_frag.picking_color = in_vert[0].picking_color;
    out_frag.view_normal = N[0] * vec4( 0, 0, -1, 0);
    gl_Position = v0[1]; EmitVertex();
    gl_Position = v0[0]; EmitVertex();
    gl_Position = v0[3]; EmitVertex();
    gl_Position = v0[2]; EmitVertex();
    EndPrimitive();
}