#pragma once

#include <core/vector_types.h>
#include <core/math_utils.h>

// http://lolengine.net/blog/2013/07/27/rgb-to-hsv-in-glsl

inline vec3 rgb_to_hsv(vec3 c) {
    const vec4 K = vec4(0.0f, -1.0f / 3.0f, 2.0f / 3.0f, -1.0f);
    const vec4 p = mix(vec4(c.z, c.y, K.w, K.z), vec4(c.y, c.z, K.x, K.y), math::step(c.z, c.y));
    const vec4 q = mix(vec4(p.x, p.y, p.w, c.x), vec4(c.x, p.y, p.z, p.x), math::step(p.x, c.x));

    const float d = q.x - math::min(q.w, q.y);
    const float e = 1.0e-10f;
    return vec3(math::abs(q.z + (q.w - q.y) / (6.0f * d + e)), d / (q.x + e), q.x);
}

inline vec3 hsv_to_rgb(vec3 c) {
    vec4 K = vec4(1.0f, 2.0f / 3.0f, 1.0f / 3.0f, 3.0f);
    vec3 p = math::abs(math::fract(vec3(c.x) + vec3(K)) * 6.0f - vec3(K.w));
    return c.z * math::mix(vec3(K.x), math::clamp(p - vec3(K.x), 0.0f, 1.0f), c.y);
}

inline vec3 hcl_to_rgb(vec3 HCL) {
    constexpr float HCLgamma = 3;
    constexpr float HCLy0 = 100;
    constexpr float HCLmaxL = 0.530454533953517f;  // == exp(HCLgamma / HCLy0) - 0.5

    vec3 RGB = vec3(0);
    if (HCL.z != 0) {
        float H = HCL.x;
        float C = HCL.y;
        float L = HCL.z * HCLmaxL;
        float Q = math::exp((1 - C / (2 * L)) * (HCLgamma / HCLy0));
        float U = (2 * L - C) / (2 * Q - 1);
        float V = C / Q;
        float T = math::tan((H + math::min(math::fract(2 * H) / 4.f, math::fract(-2 * H) / 8.f)) * math::PI * 2);
        H *= 6;
        if (H <= 1) {
            RGB.r = 1;
            RGB.g = T / (1 + T);
        } else if (H <= 2) {
            RGB.r = (1 + T) / T;
            RGB.g = 1;
        } else if (H <= 3) {
            RGB.g = 1;
            RGB.b = 1 + T;
        } else if (H <= 4) {
            RGB.g = 1 / (1 + T);
            RGB.b = 1;
        } else if (H <= 5) {
            RGB.r = -1 / T;
            RGB.b = 1;
        } else {
            RGB.r = 1;
            RGB.b = -T;
        }
        RGB = RGB * V + U;
    }
    return RGB;
}

inline vec3 rgb_to_hcl(vec3 rgb) {
    constexpr float HCLgamma = 3;
    constexpr float HCLy0 = 100;
    constexpr float HCLmaxL = 0.530454533953517f;  // == exp(HCLgamma / HCLy0) - 0.5

    vec3 HCL;
    float H = 0;
    float U = math::min(rgb.x, math::min(rgb.y, rgb.z));
    float V = math::max(rgb.x, math::max(rgb.y, rgb.z));
    float Q = HCLgamma / HCLy0;
    HCL.y = V - U;
    if (HCL.y != 0) {
        H = math::atan(rgb.g - rgb.b, rgb.r - rgb.g) / math::PI;
        Q *= U / V;
    }
    Q = exp(Q);
    HCL.x = math::fract(H / 2.f - math::min(math::fract(H), math::fract(-H)) / 6.f);
    HCL.y *= Q;
    HCL.z = math::lerp(-U, V, Q) / (HCLmaxL * 2);
    return HCL;
}

// clang-format off
inline vec3 rgb_to_XYZ(vec3 rgb) {
    constexpr mat3 RGB_2_XYZ = {0.4124564, 0.3575761, 0.1804375,
                                0.2126729, 0.7151522, 0.0721750,
                                0.0193339, 0.1191920, 0.9503041};
    return RGB_2_XYZ * rgb;
}

inline vec3 XYZ_to_rgb(vec3 XYZ) {
    const mat3 XYZ_2_RGB = { 3.2404542, -1.5371385, -0.4985314,
                            -0.9692660,  1.8760108,  0.0415560,
                             0.0556434, -0.2040259,  1.0572252};
    return XYZ_2_RGB * XYZ;
}
// clang-format on

inline vec3 XYZ_to_Lab(vec3 XYZ) {
    const auto f = [](float t) {
        const float d = 6.f / 29.f;
        return t > d * d * d ? powf(t, 1.0f / 3.0f) : (t / (3.f * d * d) + 4.f / 29.f);
    };

    const float Xn = 0.950489f;  // reference white
    const float Yn = 1.0f;
    const float Zn = 0.825188f;
    const float fx = f(XYZ.x / Xn);
    const float fy = f(XYZ.y / Yn);
    const float fz = f(XYZ.z / Zn);
    const float L = 116.f * fy - 16.f;  // maximum L = 100
    const float a = 500.f * (fx - fy);
    const float b = 200.f * (fy - fz);

    return {L, a, b};
}

inline vec3 Lab_to_XYZ(vec3 Lab) {
    const auto f = [](float t) {
        const float d = 6.f / 29.f;
        return t > d ? t * t * t : 3.0f * d * d * (t - 4.f / 29.f);
    };

    const float Xn = 0.950489f;  // reference white
    const float Yn = 1.0f;
    const float Zn = 0.825188f;
    const float X = Xn * f((Lab.x + 16.f) / 116.f + Lab.y / 500.f);
    const float Y = Yn * f((Lab.x + 16.f) / 116.f);
    const float Z = Zn * f((Lab.x + 16.f) / 116.f - Lab.z / 200.f);

    return {X, Y, Z};
}

inline vec3 rgb_to_Lab(vec3 rgb) { return XYZ_to_Lab(rgb_to_XYZ(rgb)); }
inline vec3 Lab_to_rgb(vec3 Lab) { return XYZ_to_rgb(Lab_to_XYZ(Lab)); }

inline vec3 hcl_to_rgb(float h, float c, float l) { return hcl_to_rgb({h, c, l}); }
inline vec3 rgb_to_hcl(float r, float g, float b) { return rgb_to_hcl({r, g, b}); }

inline vec4 convert_color(u32 color) { return glm::unpackUnorm4x8(color); }
inline u32 convert_color(vec4 color) { return glm::packUnorm4x8(color); }

}  // namespace math
