#version 150 core

uniform sampler2D u_texture_atom_idx;
uniform usamplerBuffer u_buffer_selection;

in vec2 tc;
out vec4 out_frag;

uint unpackUnorm4x8(in vec4 v) {
    uvec4 uv = uvec4(v * 255.f);
    return uv.x | (uv.y << 8) | (uv.z << 16) | (uv.w << 24);
}

float selected(ivec2 coord) {
    vec4 c = texelFetch(u_texture_atom_idx, ivec2(gl_FragCoord.xy) + coord, 0);
    uint atom_idx = unpackUnorm4x8(c);
    return texelFetch(u_buffer_selection, int(atom_idx)).x;
}

void main() {
    float c = selected(ivec2(0,0));
    float xn1 = selected(ivec2(-1,0));
    //float xn2 = selected(ivec2(-2,0));
    float xp1 = selected(ivec2(+1,0));
    //float xp2 = selected(ivec2(+2,0));
    float yn1 = selected(ivec2(0,-1));
    //float yn2 = selected(ivec2(0,-2));
    float yp1 = selected(ivec2(0,+1));
    //float yp2 = selected(ivec2(0,+2));

    float line_t = max(0, -xn1 -xp1 -yn1 -yp1 + 4*c);
    float fill_t = c;

    const vec3 line_color = vec3(10,10,0);
    const vec3 fill_color = vec3(2,2,2);

    vec3 color = mix(vec3(0), line_color, line_t) + mix(vec3(0), fill_color, fill_t);
    out_frag = vec4(color, 0);
}