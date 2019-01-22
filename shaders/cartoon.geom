#version 330 core

uniform mat4 u_normal_mat;
uniform mat4 u_view_proj_mat;
uniform vec2 u_scale = vec2(1.0, 1.0);

layout(lines) in;
layout(triangle_strip, max_vertices = 26) out;

in Vertex {
    vec4 control_point;
    vec4 support_vector;
    vec4 support_tangent;
    vec4 color;
    vec4 picking_color;
    vec4 weights;
} in_vert[];

out Fragment {
    smooth vec4 color;
    smooth vec4 view_normal;
    flat vec4 picking_color;
} out_frag;

void emit_vertex(in vec4 pos, in vec4 normal, in int idx) {
    out_frag.color = in_vert[idx].color;
    out_frag.picking_color = in_vert[idx].picking_color;
    out_frag.view_normal = u_normal_mat * normal;
    gl_Position = pos;
    EmitVertex();
}

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

    const float TWO_PI = 2.0 * 3.14159265;

    mat4 M[2];
    M[0] = mat4(x[0], y[0], z[0], p[0]);
    M[1] = mat4(x[1], y[1], z[1], p[1]);

    mat4 N[2];
    N[0] = u_normal_mat * M[0];
    N[1] = u_normal_mat * M[1];

    vec2 coil[12];
    vec2 sheet[12];
    vec2 helix[12];

    for (int i = 0; i < 12; i++) {
        float t = float(i) / 12.0 * TWO_PI;
        coil[i] = vec2(cos(t), sin(t)) * u_scale * 0.2;
    }

    vec2 sheet_scale = u_scale * vec2(1, 0.2);
    sheet[0] = vec2(1,-1) * sheet_scale;
    sheet[1] = vec2(1,-1) * sheet_scale;
    sheet[2] = vec2(1,-1) * sheet_scale;

    sheet[3] = vec2(1,1) * sheet_scale;
    sheet[4] = vec2(1,1) * sheet_scale;
    sheet[5] = vec2(1,1) * sheet_scale;

    sheet[6] = vec2(-1,1) * sheet_scale;
    sheet[7] = vec2(-1,1) * sheet_scale;
    sheet[8] = vec2(-1,1) * sheet_scale;

    sheet[9]  = vec2(-1,-1) * sheet_scale;
    sheet[10] = vec2(-1,-1) * sheet_scale;
    sheet[11] = vec2(-1,-1) * sheet_scale;

    for (int i = 0; i < 12; i++) {
        float t = float(i) / 12.0 * TWO_PI;
        helix[i] = vec2(cos(t), sin(t)) * vec2(1, 0.2) * u_scale;
    }

    vec4 w0 = in_vert[0].weights;
    vec4 w1 = in_vert[1].weights;

    vec2 v0[12];
    vec2 n0[12];
    for (int i = 0; i < 12; i++) {
        v0[i] = w0.x * coil[i] + w0.y * sheet[i] + w0.z * helix[i];
    }

    vec2 v1[12];
    vec2 n1[12];
    for (int i = 0; i < 12; i++) {
        v1[i] = w1.x * coil[i] + w1.y * sheet[i] + w1.z * helix[i];
    }

    vec4 pv0[12];
    vec4 nv0[12];
    for (int i = 0; i < 12; i++) {
        pv0[i] = u_view_proj_mat * M[0] * vec4(v0[i], 0, 1);
        nv0[i] = normalize(M[0] * vec4(v0[i], 0, 0));
    }

    vec4 pv1[12];
    vec4 nv1[12];
    for (int i = 0; i < 12; i++) {
        pv1[i] = u_view_proj_mat * M[1] * vec4(v1[i], 0, 1);
        nv1[i] = normalize(M[1] * vec4(v1[i], 0, 0));
    }
 
    for (int i = 0; i < 12; i++) {
        emit_vertex(pv1[i], nv1[i], 1);
        emit_vertex(pv0[i], nv0[i], 0);
    }
    emit_vertex(pv1[0], nv1[0], 1);
    emit_vertex(pv0[0], nv0[0], 0);
}