#pragma once

#define NOMINMAX

#include <core/types.h>
#include <core/vector_types.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/random.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/spline.hpp>
#include <glm/gtx/quaternion.hpp>
#include <stdlib.h>

namespace math {
constexpr float PI = 3.14159265358979323846264338327950288f;
constexpr float SQRT_TWO = 1.41421356237309504880f;
constexpr float EPSILON = 1.192092896e-07f;
constexpr float FLOAT_MAX = 3.402823466e+38f;

// @Note: The only reason here for using templates is to support vectors as well...
template <typename T>
constexpr T rad_to_deg(const T& rad) {
    return rad * (180.0f / PI);
}

template <typename T>
constexpr T deg_to_rad(const T& deg) {
    return deg * (PI / 180.0f);
}

// Core
using glm::abs;
using glm::ceil;
using glm::clamp;
using glm::exp;
using glm::floor;
using glm::fract;
using glm::inversesqrt;
using glm::max;
using glm::min;
using glm::mod;
using glm::modf;
using glm::pow;
using glm::round;
using glm::sign;
using glm::sqrt;
using glm::step;

// Trigonometry
using glm::acos;
using glm::acosh;
using glm::asin;
using glm::asinh;
using glm::atan;
using glm::atanh;
using glm::cos;
using glm::radians;
using glm::sin;
using glm::tan;
using glm::tanh;

// Vector
using glm::cross;
using glm::distance;
using glm::distance2;
using glm::dot;
using glm::length;
using glm::length2;
using glm::normalize;

// Quaternion
using glm::angle;
using glm::axis;
using glm::conjugate;
using glm::slerp;
using glm::squad;
using glm::intermediate;

// Matrix
using glm::determinant;
using glm::inverse;
using glm::transpose;

// Cast
using glm::mat3_cast;
using glm::mat4_cast;
using glm::quat_cast;

template <typename T>
auto angle(const T& a, const T& b) {
    return acos(dot(normalize(a), normalize(b)));
}

template <typename T>
auto angle(const T& a, const T& b, const T& c) {
    return angle(a - b, c - b);
}

template <typename T>
auto euler_angles(const T& a) {
    return glm::eulerAngles(a);
}

inline float dihedral_angle(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3) {
    const vec3 b1 = p1 - p0;
    const vec3 b2 = p2 - p1;
    const vec3 b3 = p3 - p2;
    const vec3 c1 = math::cross(b1, b2);
    const vec3 c2 = math::cross(b2, b3);
    return glm::atan(glm::dot(glm::cross(c1, c2), glm::normalize(b2)), glm::dot(c1, c2));
}

inline float dihedral_angle(const vec3 p[4]) { return dihedral_angle(p[0], p[1], p[2], p[3]); }

// Interpolation
template <typename T, typename V>
T lerp(T const& a, T const& b, V t) {
    return ((V)1.0 - t) * a + t * b;
}

template <typename T, typename V>
T mix(T const& a, T const& b, V t) {
    return lerp(a, b, t);
}

inline quat nlerp(const quat& q0, const quat& q1, float s) { return normalize(lerp(q0, dot(q0, q1) < 0.0f ? -q1 : q1, s)); }

inline quat cubic_nlerp(const quat& q0, const quat& q1, const quat& q2, const quat& q3, float s) {
    // @NOTE: Avoid doing cubic interpolation if q1 and q2 are very similar.
    // This avoids stability issues which come from computing intermediate quaternions from very close quaternions
    if (dot(q1, q2) > 0.95f) {
        return nlerp(q1, q2, s);
    }

    const quat sq0 = dot(q0, q1) < 0.0f ? -q0 : q0;
    const quat sq2 = dot(q1, q2) < 0.0f ? -q2 : q2;
    const quat sq3 = dot(sq2, q3) < 0.0f ? -q3 : q3;

    const auto i1 = normalize(intermediate(sq0, q1, sq2));
    const auto i2 = normalize(intermediate(q1, sq2, sq3));

    const float s2 = 2.0f * s * (1.0f - s);
    return normalize(lerp(normalize(lerp(q1, sq2, s)), normalize(lerp(i1, i2, s)), s2));
}

/*
inline quat intermediate(const quat& prev, const quat& curr, const quat& next) {
    const quat inv = inverse(curr);
    const quat a = log(next * inv);
    const quat b = log(prev * inv);
    return exp((a + b) * -0.25f) * curr;
}
*/

inline quat cubic_slerp(const quat& q0, const quat& q1, const quat& q2, const quat& q3, float s) {
    // @NOTE: Avoid doing cubic interpolation if q1 and q2 are very similar.
    // This avoids stability issues which come from computing intermediate quaternions from very close quaternions
    if (dot(q1, q2) > 0.95f) {
        return nlerp(q1, q2, s);
    }

    const quat sq0 = dot(q0, q1) < 0.0f ? -q0 : q0;
    const quat sq2 = dot(q1, q2) < 0.0f ? -q2 : q2;
    const quat sq3 = dot(sq2, q3) < 0.0f ? -q3 : q3;

    const auto i1 = normalize(intermediate(sq0, q1, sq2));
    const auto i2 = normalize(intermediate(q1, sq2, sq3));

    return squad(q1, sq2, i1, i2, s);
}

template <typename T, typename V>
V catmull_rom(const T& v1, const T& v2, const T& v3, const T& v4, V s) {
    return glm::catmullRom(v1, v2, v3, v4, s);
}

template <typename T, typename V>
T cubic_spline(const T& p0, const T& p1, const T& p2, const T& p3, V s, V tension = (V)0.5) {
    T v0 = (p2 - p0) * tension;
    T v1 = (p3 - p1) * tension;
    V s2 = s * s;
    V s3 = s * s2;

    return ((V)2.0 * p1 - (V)2.0 * p2 + v0 + v1) * s3 + (-(V)3.0 * p1 + (V)3.0 * p2 - (V)2.0 * v0 - v1) * s2 + v0 * s + p1;
}

template <typename T, typename V>
T cubic_spline_tangent(const T& p0, const T& p1, const T& p2, const T& p3, V s, V tension = (V)0.5) {
    T v0 = (p2 - p0) * tension;
    T v1 = (p3 - p1) * tension;

    // f(t) = (2p1 - 2p2 + v0 + v1)s^3 + (-3p1 + 3p2 - 2v0 - v1)s^2 + v0s + p1;
    // df(t)/dt = (2p1 - 2p2 + v0 + v1)*3s^2 + (-3p1 + 3p2 - 2v0 - v1)*2s + v0;
    return ((V)2.0 * p1 - (V)2.0 * p2 + v0 + v1) * (V)3.0 * s * s + (-(V)3.0 * p1 + (V)3.0 * p2 - (V)2.0 * v0 - v1) * (V)2.0 * s + v0;
}

// Quaternion
inline quat angle_axis(float angle, const vec3& axis) { return glm::angleAxis(angle, axis); }

// Projection
inline vec3 unproject(const vec3& window_coords, const mat4& inv_view_proj_mat, const vec4& viewport) {
    vec4 tmp = vec4(window_coords, 1.f);
    tmp.x = (tmp.x - viewport[0]) / viewport[2];
    tmp.y = (tmp.y - viewport[1]) / viewport[3];
    tmp = tmp * 2.f - 1.f;

    vec4 obj = inv_view_proj_mat * tmp;
    obj /= obj.w;

    return vec3(obj);
}

// Barycentric
inline vec3 cartesian_to_barycentric(const vec2& a, const vec2& b, const vec2& c, const vec2& cartesian) {
    const vec2 v0 = b - a;
    const vec2 v1 = c - a;
    const vec2 v2 = cartesian - a;
    const float inv_denom = v0.x * v1.y - v1.x * v0.y;
    const float v = (v2.x * v1.y - v1.x * v2.y) * inv_denom;
    const float w = (v0.x * v2.y - v2.x * v0.y) * inv_denom;
    const float u = 1.0f - v - w;
    return {u, v, w};
}

template <typename T>
inline T barycentric_to_cartesian(const T& a, const T& b, const T& c, const vec3& barycentric) {
    return a * barycentric[0] + b * barycentric[1] + c * barycentric[2];
}

template <typename T>
inline vec3 cartesian_to_barycentric(const T& a, const T& b, const T& c, const T& cartesian) {
    const T v0 = b - a;
    const T v1 = c - a;
    const T v2 = cartesian - a;
    const float d00 = dot(v0, v0);
    const float d01 = dot(v0, v1);
    const float d11 = dot(v1, v1);
    const float d20 = dot(v2, v0);
    const float d21 = dot(v2, v1);
    const float inv_denom = d00 * d11 - d01 * d01;
    const float v = (d11 * d20 - d01 * d21) * inv_denom;
    const float w = (d00 * d21 - d01 * d20) * inv_denom;
    const float u = 1.0f - v - w;
    return {u, v, w};
}

// Random
inline float rnd() { return (float)rand() / (float)RAND_MAX; }
inline void set_rnd_seed(unsigned int seed) { srand(seed); }

inline float halton(int index, int base) {
    float f = 1;
    float r = 0;
    const float ifb = 1.f / base;
    while (index > 0) {
        f = f * ifb;
        r = r + f * fmodf((float)index, (float)base);
        index = (int)(index * ifb);
    }
    return r;
}

void generate_halton_sequence(float* dst, int count, int base);
void generate_halton_sequence(vec2* dst, int count, int base_x, int base_y);

// Color

// http://lolengine.net/blog/2013/07/27/rgb-to-hsv-in-glsl

inline vec3 rgb_to_hsv(vec3 c) {
    vec4 K = vec4(0.0f, -1.0f / 3.0f, 2.0f / 3.0f, -1.0f);
    vec4 p = mix(vec4(c.b, c.g, K.w, K.z), vec4(c.g, c.b, K.x, K.y), step(c.b, c.g));
    vec4 q = mix(vec4(p.x, p.y, p.w, c.r), vec4(c.r, p.y, p.z, p.x), step(p.x, c.r));

    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10f;
    return vec3(abs(q.z + (q.w - q.y) / (6.0f * d + e)), d / (q.x + e), q.x);
}

inline vec3 hsv_to_rgb(vec3 c) {
    vec4 K = vec4(1.0f, 2.0f / 3.0f, 1.0f / 3.0f, 3.0f);
    vec3 p = abs(fract(vec3(c.x) + vec3(K)) * 6.0f - vec3(K.w));
    return c.z * mix(vec3(K.x), clamp(p - vec3(K.x), 0.0f, 1.0f), c.y);
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
        float Q = exp((1 - C / (2 * L)) * (HCLgamma / HCLy0));
        float U = (2 * L - C) / (2 * Q - 1);
        float V = C / Q;
        float T = tan((H + min(fract(2 * H) / 4.f, fract(-2 * H) / 8.f)) * PI * 2);
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
    float U = min(rgb.r, min(rgb.g, rgb.b));
    float V = max(rgb.r, max(rgb.g, rgb.b));
    float Q = HCLgamma / HCLy0;
    HCL.y = V - U;
    if (HCL.y != 0) {
        H = atan(rgb.g - rgb.b, rgb.r - rgb.g) / PI;
        Q *= U / V;
    }
    Q = exp(Q);
    HCL.x = fract(H / 2.f - min(fract(H), fract(-H)) / 6.f);
    HCL.y *= Q;
    HCL.z = lerp(-U, V, Q) / (HCLmaxL * 2);
    return HCL;
}

inline vec3 rgb_to_XYZ(vec3 rgb) {
    const mat3 RGB_2_XYZ = {0.4124564, 0.3575761, 0.1804375, 0.2126729, 0.7151522, 0.0721750, 0.0193339, 0.1191920, 0.9503041};
    return RGB_2_XYZ * rgb;
}

inline vec3 XYZ_to_rgb(vec3 XYZ) {
    const mat3 XYZ_2_RGB = {3.2404542, -1.5371385, -0.4985314, -0.9692660, 1.8760108, 0.0415560, 0.0556434, -0.2040259, 1.0572252};
    return XYZ_2_RGB * XYZ;
}

inline vec3 XYZ_to_Lab(vec3 XYZ) {
    const auto f = [](float t) {
        const float d = 6.f / 29.f;
        return t > d * d * d ? powf(t, 0.33333333f) : (t / (3.f * d * d) + 4.f / 29.f);
    };

    const float Xn = 0.950489f;  // reference white
    const float Yn = 1.0f;
    const float Zn = 0.825188f;
    const float L = 116.f * f(XYZ.y / Yn) - 16.f;             // maximum L = 100
    const float a = 500.f * (f(XYZ.x / Xn) - f(XYZ.y / Yn));  // maximum
    const float b = 200.f * (f(XYZ.y / Yn) - f(XYZ.z / Zn));

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

inline vec4 convert_color(uint32 color) { return glm::unpackUnorm4x8(color); }
inline uint32 convert_color(vec4 color) { return glm::packUnorm4x8(color); }

}  // namespace math
