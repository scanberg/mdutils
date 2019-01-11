#version 330 core
//#define DEBUG_MODE 1
uniform mat4 u_view_mat;
uniform mat4 u_view_proj_mat;
uniform vec2 u_scale = vec2(1.0, 0.1);
layout(lines) in;
#ifdef DEBUG_MODE
layout(line_strip, max_vertices = 12) out;
#else
layout(triangle_strip, max_vertices = 24) out;
#endif

in Vertex {
    vec4 control_point;
    vec4 support_vector;
    vec4 support_tangent;
    vec4 color;
    vec4 picking_color;
    uint cap;
} in_vert[];

out Fragment {
    smooth vec4 color;
    smooth vec3 view_normal;
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

    mat4 m[2];
    m[0] = mat4(x[0], y[0], z[0], p[0]);
    m[1] = mat4(x[1], y[1], z[1], p[1]);

#ifdef DEBUG_MODE
    out_frag.color = vec4(1,0,0,1);
    gl_Position = m[0] * vec4(0,0,0,1); EmitVertex();
    gl_Position = m[0] * vec4(1,0,0,1); EmitVertex();
    EndPrimitive();
    gl_Position = m[1] * vec4(0,0,0,1); EmitVertex();
    gl_Position = m[1] * vec4(1,0,0,1); EmitVertex();
    EndPrimitive();
    out_frag.color = vec4(0,1,0,1);
    gl_Position = m[0] * vec4(0,0,0,1); EmitVertex();
    gl_Position = m[0] * vec4(0,1,0,1); EmitVertex();
    EndPrimitive();
    gl_Position = m[1] * vec4(0,0,0,1); EmitVertex();
    gl_Position = m[1] * vec4(0,1,0,1); EmitVertex();
    EndPrimitive();
    out_frag.color = vec4(0,0,1,1);
    gl_Position = m[0] * vec4(0,0,0,1); EmitVertex();
    gl_Position = m[0] * vec4(0,0,1,1); EmitVertex();
    EndPrimitive();
    gl_Position = m[1] * vec4(0,0,0,1); EmitVertex();
    gl_Position = m[1] * vec4(0,0,1,1); EmitVertex();
    EndPrimitive();

#else   
    mat3 normal_mat = inverse(transpose(mat3(u_view_mat)));

    // BOTTOM
    out_frag.color = in_vert[0].color;
    out_frag.picking_color = in_vert[0].picking_color;
    out_frag.view_normal = normal_mat * vec3(m[0] * vec4( 0, -1, 0, 0));
    gl_Position = u_view_proj_mat * m[0] * vec4(vec2(-1,-1) * u_scale, 0, 1); EmitVertex();
    gl_Position = u_view_proj_mat * m[0] * vec4(vec2( 1,-1) * u_scale, 0, 1); EmitVertex();
    
    out_frag.color = in_vert[1].color;
    out_frag.picking_color = in_vert[1].picking_color;
    out_frag.view_normal = normal_mat * vec3(m[1] * vec4( 0, -1, 0, 0));
    gl_Position = u_view_proj_mat * m[1] * vec4(vec2(-1,-1) * u_scale, 0, 1); EmitVertex();
    gl_Position = u_view_proj_mat * m[1] * vec4(vec2( 1,-1) * u_scale, 0, 1); EmitVertex();
    EndPrimitive();

    // TOP
    out_frag.color = in_vert[0].color;
    out_frag.picking_color = in_vert[0].picking_color;
    out_frag.view_normal = normal_mat * vec3(m[0] * vec4( 0, 1, 0, 0));
    gl_Position = u_view_proj_mat * m[0] * vec4(vec2( 1, 1) * u_scale, 0, 1); EmitVertex();
    gl_Position = u_view_proj_mat * m[0] * vec4(vec2(-1, 1) * u_scale, 0, 1); EmitVertex();

    out_frag.color = in_vert[1].color;
    out_frag.picking_color = in_vert[1].picking_color;
    out_frag.view_normal = normal_mat * vec3(m[1] * vec4( 0, 1, 0, 0));
    gl_Position = u_view_proj_mat * m[1] * vec4(vec2( 1, 1) * u_scale, 0, 1); EmitVertex();
    gl_Position = u_view_proj_mat * m[1] * vec4(vec2(-1, 1) * u_scale, 0, 1); EmitVertex();
    EndPrimitive();

    // LEFT
    out_frag.color = in_vert[0].color;
    out_frag.picking_color = in_vert[0].picking_color;
    out_frag.view_normal = normal_mat * vec3(m[0] * vec4(-1, 0, 0, 0));
    gl_Position = u_view_proj_mat * m[0] * vec4(vec2(-1, 1) * u_scale, 0, 1); EmitVertex();
    gl_Position = u_view_proj_mat * m[0] * vec4(vec2(-1,-1) * u_scale, 0, 1); EmitVertex();

    out_frag.color = in_vert[1].color;
    out_frag.picking_color = in_vert[1].picking_color;
    out_frag.view_normal = normal_mat * vec3(m[1] * vec4(-1, 0, 0, 0));
    gl_Position = u_view_proj_mat * m[1] * vec4(vec2(-1, 1) * u_scale, 0, 1); EmitVertex();
    gl_Position = u_view_proj_mat * m[1] * vec4(vec2(-1,-1) * u_scale, 0, 1); EmitVertex();
    EndPrimitive();

    // RIGHT
    out_frag.color = in_vert[0].color;
    out_frag.picking_color = in_vert[0].picking_color;
    out_frag.view_normal = normal_mat * vec3(m[0] * vec4( 1, 0, 0, 0));
    gl_Position = u_view_proj_mat * m[0] * vec4(vec2( 1,-1) * u_scale, 0, 1); EmitVertex();
    gl_Position = u_view_proj_mat * m[0] * vec4(vec2( 1, 1) * u_scale, 0, 1); EmitVertex();

    out_frag.color = in_vert[1].color;
    out_frag.picking_color = in_vert[1].picking_color;
    out_frag.view_normal = normal_mat * vec3(m[1] * vec4( 1, 0, 0, 0));
    gl_Position = u_view_proj_mat * m[1] * vec4(vec2( 1,-1) * u_scale, 0, 1); EmitVertex();
    gl_Position = u_view_proj_mat * m[1] * vec4(vec2( 1, 1) * u_scale, 0, 1); EmitVertex();
    EndPrimitive();

    // FRONT
    if (in_vert[0].cap != 0U) {
        out_frag.color = in_vert[1].color;
        out_frag.picking_color = in_vert[1].picking_color;
        out_frag.view_normal = normal_mat * vec3(m[1] * vec4( 0, 0, 1, 0));
        gl_Position = u_view_proj_mat * m[1] * vec4(vec2(-1,-1) * u_scale, 0, 1); EmitVertex();
        gl_Position = u_view_proj_mat * m[1] * vec4(vec2( 1,-1) * u_scale, 0, 1); EmitVertex();
        gl_Position = u_view_proj_mat * m[1] * vec4(vec2(-1, 1) * u_scale, 0, 1); EmitVertex();
        gl_Position = u_view_proj_mat * m[1] * vec4(vec2( 1, 1) * u_scale, 0, 1); EmitVertex();
        EndPrimitive();
    }

    // BACK
    if (in_vert[1].cap != 0U) {
        out_frag.color = in_vert[0].color;
        out_frag.picking_color = in_vert[0].picking_color;
        out_frag.view_normal = normal_mat * vec3(m[0] * vec4( 0, 0, -1, 0));
        gl_Position = u_view_proj_mat * m[0] * vec4(vec2( 1,-1) * u_scale, 0, 1); EmitVertex();
        gl_Position = u_view_proj_mat * m[0] * vec4(vec2(-1,-1) * u_scale, 0, 1); EmitVertex();
        gl_Position = u_view_proj_mat * m[0] * vec4(vec2( 1, 1) * u_scale, 0, 1); EmitVertex();
        gl_Position = u_view_proj_mat * m[0] * vec4(vec2(-1, 1) * u_scale, 0, 1); EmitVertex();
        EndPrimitive();
    }

#endif
}