#pragma once

#include <core/types.h>

//#define GLM_FORCE_INTRINSICS
//#define GLM_FORCE_DEFAULT_ALIGNED_GENTYPES
//#define GLM_FORCE_SILENT_WARNINGS // silence warnings induced by glm
#define GLM_FORCE_XYZW_ONLY

#include <glm/detail/type_vec2.hpp>
#include <glm/detail/type_vec3.hpp>
#include <glm/detail/type_vec4.hpp>
#include <glm/detail/type_quat.hpp>
#include <glm/detail/type_mat2x2.hpp>
#include <glm/detail/type_mat3x3.hpp>
#include <glm/detail/type_mat4x4.hpp>

using quat = glm::tquat<float>;
using vec2 = glm::tvec2<float>;
using vec3 = glm::tvec3<float>;
using vec4 = glm::tvec4<float>;

using f32vec2 = glm::tvec2<float>;
using f32vec3 = glm::tvec3<float>;
using f32vec4 = glm::tvec4<float>;

using dvec2 = glm::tvec2<double>;
using dvec3 = glm::tvec3<double>;
using dvec4 = glm::tvec4<double>;

using f64vec2 = glm::tvec2<double>;
using f64vec3 = glm::tvec3<double>;
using f64vec4 = glm::tvec4<double>;

using ivec2 = glm::tvec2<int32_t>;
using ivec3 = glm::tvec3<int32_t>;
using ivec4 = glm::tvec4<int32_t>;

using i8vec2 = glm::tvec2<int8_t>;
using i8vec3 = glm::tvec3<int8_t>;
using i8vec4 = glm::tvec4<int8_t>;

using i16vec2 = glm::tvec2<int16_t>;
using i16vec3 = glm::tvec3<int16_t>;
using i16vec4 = glm::tvec4<int16_t>;

using i32vec2 = glm::tvec2<int32_t>;
using i32vec3 = glm::tvec3<int32_t>;
using i32vec4 = glm::tvec4<int32_t>;

using i64vec2 = glm::tvec2<int64_t>;
using i64vec3 = glm::tvec3<int64_t>;
using i64vec4 = glm::tvec4<int64_t>;

using uvec2 = glm::tvec2<uint32_t>;
using uvec3 = glm::tvec3<uint32_t>;
using uvec4 = glm::tvec4<uint32_t>;

using u8vec2 = glm::tvec2<uint8_t>;
using u8vec3 = glm::tvec3<uint8_t>;
using u8vec4 = glm::tvec4<uint8_t>;

using u16vec2 = glm::tvec2<uint16_t>;
using u16vec3 = glm::tvec3<uint16_t>;
using u16vec4 = glm::tvec4<uint16_t>;

using u32vec2 = glm::tvec2<uint32_t>;
using u32vec3 = glm::tvec3<uint32_t>;
using u32vec4 = glm::tvec4<uint32_t>;

using u64vec2 = glm::tvec2<uint64_t>;
using u64vec3 = glm::tvec3<uint64_t>;
using u64vec4 = glm::tvec4<uint64_t>;

using mat2 = glm::tmat2x2<float>;
using mat3 = glm::tmat3x3<float>;
using mat4 = glm::tmat4x4<float>;
using mat2x2 = glm::tmat2x2<float>;
using mat3x3 = glm::tmat3x3<float>;
using mat4x4 = glm::tmat4x4<float>;

using dmat2 = glm::tmat2x2<double>;
using dmat3 = glm::tmat3x3<double>;
using dmat4 = glm::tmat4x4<double>;
using dmat2x2 = glm::tmat2x2<double>;
using dmat3x3 = glm::tmat3x3<double>;
using dmat4x4 = glm::tmat4x4<double>;

struct soa_vec3 {
    float* __restrict x;
    float* __restrict y;
    float* __restrict z;
};

inline soa_vec3 operator+(const soa_vec3& in, i64 offset) {
    return {in.x + offset, in.y + offset, in.z + offset};
}

#ifndef HAS_VECTOR_TEMPLATE_INSTANTIATION
extern template struct glm::vec<2, float, glm::packed_highp>;
extern template struct glm::vec<3, float, glm::packed_highp>;
extern template struct glm::vec<4, float, glm::packed_highp>;

extern template struct glm::vec<2, double, glm::packed_highp>;
extern template struct glm::vec<3, double, glm::packed_highp>;
extern template struct glm::vec<4, double, glm::packed_highp>;

extern template struct glm::vec<2, int8_t, glm::packed_highp>;
extern template struct glm::vec<3, int8_t, glm::packed_highp>;
extern template struct glm::vec<4, int8_t, glm::packed_highp>;

extern template struct glm::vec<2, int16_t, glm::packed_highp>;
extern template struct glm::vec<3, int16_t, glm::packed_highp>;
extern template struct glm::vec<4, int16_t, glm::packed_highp>;

extern template struct glm::vec<2, int32_t, glm::packed_highp>;
extern template struct glm::vec<3, int32_t, glm::packed_highp>;
extern template struct glm::vec<4, int32_t, glm::packed_highp>;

extern template struct glm::vec<2, int64_t, glm::packed_highp>;
extern template struct glm::vec<3, int64_t, glm::packed_highp>;
extern template struct glm::vec<4, int64_t, glm::packed_highp>;

extern template struct glm::vec<2, uint8_t, glm::packed_highp>;
extern template struct glm::vec<3, uint8_t, glm::packed_highp>;
extern template struct glm::vec<4, uint8_t, glm::packed_highp>;

extern template struct glm::vec<2, uint16_t, glm::packed_highp>;
extern template struct glm::vec<3, uint16_t, glm::packed_highp>;
extern template struct glm::vec<4, uint16_t, glm::packed_highp>;

extern template struct glm::vec<2, uint32_t, glm::packed_highp>;
extern template struct glm::vec<3, uint32_t, glm::packed_highp>;
extern template struct glm::vec<4, uint32_t, glm::packed_highp>;

extern template struct glm::vec<2, uint64_t, glm::packed_highp>;
extern template struct glm::vec<3, uint64_t, glm::packed_highp>;
extern template struct glm::vec<4, uint64_t, glm::packed_highp>;

extern template struct glm::mat<2, 2, float, glm::packed_highp>;
extern template struct glm::mat<3, 3, float, glm::packed_highp>;
extern template struct glm::mat<4, 4, float, glm::packed_highp>;

extern template struct glm::mat<2, 2, double, glm::packed_highp>;
extern template struct glm::mat<3, 3, double, glm::packed_highp>;
extern template struct glm::mat<4, 4, double, glm::packed_highp>;
#endif