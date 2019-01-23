#version 150 core

uniform sampler2D u_texture_atom_idx;
uniform samplerBuffer u_buffer_selection;

in vec2 tc;
out vec4 out_frag;

uint unpackUnorm4x8(in vec4 v) {
    uvec4 uv = uvec4(v * 255);
    return uv.x | (uv.y << 8) | (uv.z << 16) | (uv.w << 24);
}

float selected(ivec2 coord) {
    vec4 c = texelFetch(u_texture_atom_idx, ivec2(gl_FragCoord.xy) + coord, 0);
    uint atom_idx = unpackUnorm4x8(c);
    return texelFetch(u_buffer_selection, int(atom_idx)).x;
}

void main() {
    float c = selected(ivec2(0,0));
    float xn = selected(ivec2(-1,0));
    float xp = selected(ivec2(+1,0));
    float yn = selected(ivec2(0,-1));
    float yp = selected(ivec2(0,+1));

    float t = max(0, xn + xp + yn + yp - 4*c);
    vec3 color = mix(vec3(0), vec3(5), t);
    out_frag = vec4(color, 0);
}