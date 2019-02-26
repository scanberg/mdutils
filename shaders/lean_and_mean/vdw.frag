#version 150 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_explicit_attrib_location : enable

uniform vec4 u_color = vec4(1,1,1,1);

in GS_FS {
    smooth vec2 uv;
} in_frag;

layout(location = 0) out vec4 out_color;

void main() {
    if (dot(in_frag.uv, in_frag.uv) > 1.0) discard;
    out_color = u_color;
}