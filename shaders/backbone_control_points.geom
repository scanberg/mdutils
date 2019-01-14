#version 330 core

layout (triangles_adjacency) in;
layout (points, max_vertices = 1) out;

in uint atom_index[];

out vec3 out_control_point;
out vec3 out_support_vector;
out vec3 out_tangent_vector;
out vec2 out_backbone_angles;
out uint out_atom_index;

float dihedral_angle(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3) {
    vec3 b1 = p1 - p0;
    vec3 b2 = p2 - p1;
    vec3 b3 = p3 - p2;
    vec3 c1 = cross(b1, b2);
    vec3 c2 = cross(b2, b3);
    return atan(dot(cross(c1, c2), normalize(b2)), dot(c1, c2));
}

void main() {
    vec3 ca  = gl_in[0].gl_Position.xyz; // Ca[i]
    vec3 c   = gl_in[1].gl_Position.xyz; // C[i]
    vec3 o   = gl_in[2].gl_Position.xyz; // O[i]
    vec3 n   = gl_in[3].gl_Position.xyz; // N[i]
    vec3 c_p = gl_in[4].gl_Position.xyz; // C[i-1]
    vec3 n_n = gl_in[5].gl_Position.xyz; // N[i+1]

    vec3 p = ca;
    vec3 v = normalize(o - c);
    vec3 t = vec3(0);   // Placeholder, tangent is computed analytically in spline shader

    float phi = dihedral_angle(c_p, n, ca, c);
    float psi = dihedral_angle(n, ca, c, n_n);

    out_control_point = p;
    out_support_vector = v;
    out_tangent_vector = t;
    out_backbone_angles = vec2(phi, psi);
    out_atom_index = atom_index[0];
    EmitVertex();
    EndPrimitive();
}