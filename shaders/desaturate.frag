#version 150 core

uniform sampler2D u_texture_atom_color;
uniform sampler2D u_texture_atom_idx;
uniform usamplerBuffer u_buffer_selection;
uniform int u_selecting;

// http://lolengine.net/blog/2013/07/27/rgb-to-hsv-in-glsl.
vec3 rgb2hsv(vec3 c) {
    vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
    vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));

    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10;
    return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

vec3 hsv2rgb(vec3 c) {
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

in vec2 tc;
out vec4 out_frag;

uint unpackUnorm4x8(in vec4 v) {
    uvec4 uv = uvec4(v * 255.f);
    return uv.x | (uv.y << 8) | (uv.z << 16) | (uv.w << 24);
}

uint fetch_value(ivec2 coord) {
    vec4 c = texelFetch(u_texture_atom_idx, ivec2(gl_FragCoord.xy) + coord, 0);
    uint atom_idx = unpackUnorm4x8(c);
    return texelFetch(u_buffer_selection, int(atom_idx)).x;
}

void main() {
    uint c  = fetch_value(ivec2(0,0));

    float highlight_t = float((c & 1U) != 0U);
    float selection_t = float((c & 2U) != 0U);

    vec3 atom_color = texelFetch(u_texture_atom_color, ivec2(gl_FragCoord.xy), 0).rgb;
    vec3 color = atom_color;

    if (highlight_t > 0) {
        //vec3 hsv = rgb2hsv(atom_color);
        //hsv.z *= 1.5;
        //hsv.y *= 0.75;
        //hsv = mix(hsv, rgb2hsv(vec3(0, 0, 1)), 0.3);
        //hsv.y *= 1.5;
        //color = hsv2rgb(hsv); 
    } else if (u_selecting > 0) {
        vec3 hsv = rgb2hsv(atom_color);
        //hsv.y *= mix(0.4, 1.0, max(selection_t, highlight_t));
        hsv.y *= mix(0.1, 1.0, selection_t);
        color = hsv2rgb(hsv);
    }
    //color = atom_color;
    out_frag = vec4(color, 1);
}