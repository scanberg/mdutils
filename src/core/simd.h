#pragma once

// https://stackoverflow.com/questions/11228855/header-files-for-x86-simd-intrinsics

#if defined(_MSC_VER)
/* Microsoft C/C++-compatible compiler */
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
/* GCC-compatible compiler, targeting x86/x86-64 */
#include <x86intrin.h>
#elif defined(__GNUC__) && defined(__ARM_NEON__)
/* GCC-compatible compiler, targeting ARM with NEON */
#include <arm_neon.h>
#elif defined(__GNUC__) && defined(__IWMMXT__)
/* GCC-compatible compiler, targeting ARM with WMMX */
#include <mmintrin.h>
#endif

#define INLINE inline

namespace simd {

typedef __m128 float4;
// typedef __m256 float8;
// typedef __m512 float16;

INLINE float4 set_float4(float v) { return _mm_set1_ps(v); }
INLINE float4 set_float4(float x, float y, float z, float w) { return _mm_set_ps(w, z, y, x); }

INLINE float4 zero_float4() { return _mm_setzero_ps(); }

INLINE bool all_zero(const float4 v) {
    const auto res = _mm_cmpeq_ps(v, _mm_setzero_ps());
    const auto mask = _mm_movemask_ps(res);
    return mask == 0x0000000F;
}

INLINE float4 load(const float* addr) { return _mm_load_ps(addr); }
INLINE void store(float* addr, float4 v) { _mm_store_ps(addr, v); }

INLINE float4 add(const float4 a, const float4 b) { return _mm_add_ps(a, b); }
INLINE float4 sub(const float4 a, const float4 b) { return _mm_sub_ps(a, b); }
INLINE float4 mul(const float4 a, const float4 b) { return _mm_mul_ps(a, b); }
INLINE float4 div(const float4 a, const float4 b) { return _mm_div_ps(a, b); }

INLINE float4 bit_and(const float4 a, const float4 b) { return _mm_and_ps(a, b); }
INLINE float4 bit_or(const float4 a, const float4 b) { return _mm_or_ps(a, b); }
INLINE float4 bit_xor(const float4 a, const float4 b) { return _mm_xor_ps(a, b); }

INLINE float4 cmp_gt(const float4 a, const float4 b) { return _mm_cmpgt_ps(a, b); }
INLINE float4 cmp_ge(const float4 a, const float4 b) { return _mm_cmpge_ps(a, b); }
INLINE float4 cmp_lt(const float4 a, const float4 b) { return _mm_cmplt_ps(a, b); }
INLINE float4 cmp_le(const float4 a, const float4 b) { return _mm_cmple_ps(a, b); }
INLINE float4 cmp_eq(const float4 a, const float4 b) { return _mm_cmpeq_ps(a, b); }

INLINE float4 abs(float4 a) { return _mm_and_ps(a, _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF))); }
INLINE float4 sign(float4 x) {
    const float4 zero = zero_float4();
    const float4 lz = cmp_lt(x, zero);
    const float4 gz = cmp_gt(x, zero);
    const float4 neg = bit_and(lz, set_float4(-1.0f));
    const float4 pos = bit_and(gz, set_float4(1.0f));
    const float4 res = bit_or(neg, pos);
    return res;
}

INLINE float4 step(const float4 edge, const float4 x) {
    const float4 cmp = cmp_ge(x, edge);
    const float4 res = bit_and(cmp, set_float4(1.f));
    return res;
}

INLINE __m128 lerp(const float4 a, const float4 b, float t) {
    const float4 one_minus_alpha = set_float4(1.0f - t);
    const float4 alpha = set_float4(t);
    const float4 res = add(mul(a, one_minus_alpha), mul(b, alpha));
    return res;
};

INLINE __m128 cubic_spline(const __m128 p0, const __m128 p1, const __m128 p2, const __m128 p3, float s, float tension = 0.5f) {
    const float4 s1 = set_float4(s);
    const float4 s2 = mul(s1, s1);
    const float4 s3 = mul(s2, s1);
    const float4 vt = set_float4(tension);
    const float4 v0 = mul(sub(p2, p0), vt);
    const float4 v1 = mul(sub(p3, p1), vt);
    const float4 x0 = add(mul(set_float4(2), sub(p1, p2)), add(v0, v1));
    const float4 x1 = sub(mul(set_float4(3), sub(p2, p1)), add(mul(set_float4(2), v0), v1));
    const float4 r0 = add(mul(x0, s3), mul(x1, s2));
    const float4 r1 = add(mul(v0, s1), p1);
    const float4 res = add(r0, r1);
    return res;
}

/*
INLINE float8 add(float8 a, float8 b) { return _mm256_add_ps(a, b); }
INLINE float8 sub(float8 a, float8 b) { return _mm256_sub_ps(a, b); }
INLINE float8 mul(float8 a, float8 b) { return _mm256_mul_ps(a, b); }
INLINE float8 div(float8 a, float8 b) { return _mm256_div_ps(a, b); }

INLINE float16 add(float16 a, float16 b) { return _mm512_add_ps(a, b); }
INLINE float16 sub(float16 a, float16 b) { return _mm512_sub_ps(a, b); }
INLINE float16 mul(float16 a, float16 b) { return _mm512_mul_ps(a, b); }
INLINE float16 div(float16 a, float16 b) { return _mm512_div_ps(a, b); }
*/

}  // namespace simd
