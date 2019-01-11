#version 150 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_proj_mat;
uniform float u_exposure = 1.0;
uniform float u_radius = 1.0;

in Fragment {
    flat vec4 color[2];
    flat vec4 picking_color[2];
    smooth vec3 view_pos;

    flat vec4  capsule_center_radius;
    flat vec4  capsule_axis_length;
} in_frag;

#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif
layout(location = 0) out vec4 out_color_alpha;
layout(location = 1) out vec4 out_f0_smoothness;
layout(location = 2) out vec4 out_normal;
layout(location = 3) out vec4 out_picking_id;

// Source from Ingo Quilez (https://www.shadertoy.com/view/Xt3SzX)
float intersect_capsule(in vec3 ro, in vec3 rd, in vec3 cc, in vec3 ca, float cr,
                      float ch, out vec3 normal, out int side)  // cc center, ca orientation axis, cr radius, ch height
{
    vec3 oc = ro - cc;
    ch *= 0.5;

    float card = dot(ca, rd);
    float caoc = dot(ca, oc);

    float a = 1.0 - card * card;
    float b = dot(oc, rd) - caoc * card;
    float c = dot(oc, oc) - caoc * caoc - cr * cr;
    float h = b * b - a * c;
    if (h < 0.0) return -1.0;
    float t = (-b - sqrt(h)) / a;

    float y = caoc + t * card;
    side = int(y > 0);

    // body
    if (abs(y) < ch) {
        normal = normalize(oc + t * rd - ca * y);
        return t;
    }

    // caps
    float sy = sign(y);
    oc = ro - (cc + sy * ca * ch);
    b = dot(rd, oc);
    c = dot(oc, oc) - cr * cr;
    h = b * b - c;
    if (h > 0.0) {
        t = -b - sqrt(h);
        normal = normalize(oc + rd * t);
        return t;
    }

    return -1.0;
}

vec4 encode_normal (vec3 n) {
    float p = sqrt(n.z*8+8);
    return vec4(n.xy/p + 0.5,0,0);
}

void main() {
    vec3 ro = vec3(0);
    vec3 rd = normalize(in_frag.view_pos);
    vec3 cc = in_frag.capsule_center_radius.xyz;
    float cr = in_frag.capsule_center_radius.w;
    vec3 ca = in_frag.capsule_axis_length.xyz;
    float ch = in_frag.capsule_axis_length.w;

    vec3 normal;
    int side;
    float t = intersect_capsule(ro, rd, cc, ca, cr, ch, normal, side);
    if (t < 0.0) {
        discard;
        return;
    }

    vec3 pos = ro + t*rd;
    vec4 color = in_frag.color[side];
    vec4 picking_color = in_frag.picking_color[side];

    vec4 coord = vec4(0, 0, pos.z, 1);
    coord = u_proj_mat * coord;
    coord = coord / coord.w;

    gl_FragDepth = coord.z * 0.5 + 0.5;
    out_color_alpha = color;
    out_f0_smoothness = vec4(0.04, 0.04, 0.04, 0.0);
    out_normal = encode_normal(normal);
    out_picking_id = picking_color;
}