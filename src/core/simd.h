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
typedef __m256 float8;
typedef __m512 float16;

typedef __mmask16 mask16;

INLINE float4 set_128(float v) { return _mm_set1_ps(v); }
INLINE float4 set_128(float x, float y, float z, float w) { return _mm_set_ps(w, z, y, x); }

INLINE float4 zero_128() { return _mm_setzero_ps(); }

INLINE bool all_zero(const float4 v) {
    const auto res = _mm_cmpeq_ps(v, _mm_setzero_ps());
    const auto mask = _mm_movemask_ps(res);
    return mask == 0x0000000F;
}

INLINE float4 load_128(const float* addr) { return _mm_loadu_ps(addr); }
INLINE float4 load_aligned_128(const float* addr) {
    ASSERT(IS_ALIGNED(addr, 16));
    return _mm_load_ps(addr);
}

INLINE void store(float* addr, float4 v) { _mm_storeu_ps(addr, v); }
INLINE void store_aligned(float* addr, float4 v) {
    ASSERT(IS_ALIGNED(addr, 16));
    _mm_store_ps(addr, v);
}

INLINE float4 add(const float4 a, const float4 b) { return _mm_add_ps(a, b); }
INLINE float4 sub(const float4 a, const float4 b) { return _mm_sub_ps(a, b); }
INLINE float4 mul(const float4 a, const float4 b) { return _mm_mul_ps(a, b); }
INLINE float4 div(const float4 a, const float4 b) { return _mm_div_ps(a, b); }

INLINE float4 bit_and(const float4 a, const float4 b) { return _mm_and_ps(a, b); }
INLINE float4 bit_and_not(const float4 a, const float4 b) { return _mm_andnot_ps(a, b); }
INLINE float4 bit_or(const float4 a, const float4 b) { return _mm_or_ps(a, b); }
INLINE float4 bit_xor(const float4 a, const float4 b) { return _mm_xor_ps(a, b); }

INLINE float4 cmp_gt(const float4 a, const float4 b) { return _mm_cmpgt_ps(a, b); }
INLINE float4 cmp_ge(const float4 a, const float4 b) { return _mm_cmpge_ps(a, b); }
INLINE float4 cmp_lt(const float4 a, const float4 b) { return _mm_cmplt_ps(a, b); }
INLINE float4 cmp_le(const float4 a, const float4 b) { return _mm_cmple_ps(a, b); }
INLINE float4 cmp_eq(const float4 a, const float4 b) { return _mm_cmpeq_ps(a, b); }
INLINE float4 cmp_neq(const float4 a, const float4 b) { return _mm_cmpneq_ps(a, b); }

INLINE float4 abs(float4 a) { return _mm_and_ps(a, _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF))); }

// @NOTE: If 0.0f is given as input, it will be mapped to 1.0f
INLINE float4 sign(float4 x) {
	const float4 sgn = bit_and(x, set_128(-0.0f));
	const float4 res = bit_xor(sgn, set_128(1.0f));
	return res;
}

INLINE float4 min(const float4 a, const float4 b) {
	const float4 res = _mm_min_ps(a, b);
	return res;
}

INLINE float4 max(const float4 a, const float4 b) {
	const float4 res = _mm_max_ps(a, b);
	return res;
}

INLINE float horizontal_min(const float4 x) {
	const float4 min1 = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 0, 3, 2));
	const float4 min2 = _mm_min_ps(x, min1);
	const float4 min3 = _mm_shuffle_ps(min2, min2, _MM_SHUFFLE(0, 0, 0, 1));
	const float4 min4 = _mm_min_ps(min2, min3);
	float result = _mm_cvtss_f32(min4);
	return result;
}

INLINE float horizontal_max(const float4 x) {
	const float4 max1 = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 0, 3, 2));
	const float4 max2 = _mm_max_ps(x, max1);
	const float4 max3 = _mm_shuffle_ps(max2, max2, _MM_SHUFFLE(0, 0, 0, 1));
	const float4 max4 = _mm_max_ps(max2, max3);
	float result = _mm_cvtss_f32(max4);
	return result;
}

INLINE float4 step(const float4 edge, const float4 x) {
    const float4 cmp = cmp_ge(x, edge);
    const float4 res = bit_and(cmp, set_128(1.f));
    return res;
}

INLINE float4 lerp(const float4 a, const float4 b, float t) {
    const float4 one_minus_alpha = set_128(1.0f - t);
    const float4 alpha = set_128(t);
    const float4 res = add(mul(a, one_minus_alpha), mul(b, alpha));
    return res;
}

INLINE float4 mix(const float4 a, const float4 b, float t) {
	return lerp(a, b, t);
}

INLINE float4 cubic_spline(const __m128 p0, const __m128 p1, const __m128 p2, const __m128 p3, float s, float tension = 0.5f) {
    const float4 vt = set_128(tension);
    const float4 s1 = set_128(s);
    const float4 s2 = mul(s1, s1);
    const float4 s3 = mul(s2, s1);
    const float4 v0 = mul(sub(p2, p0), vt);
    const float4 v1 = mul(sub(p3, p1), vt);
    const float4 x0 = add(mul(set_128(2), sub(p1, p2)), add(v0, v1));
    const float4 x1 = sub(mul(set_128(3), sub(p2, p1)), add(mul(set_128(2), v0), v1));
    const float4 r0 = add(mul(x0, s3), mul(x1, s2));
    const float4 r1 = add(mul(v0, s1), p1);
    const float4 res = add(r0, r1);
    return res;
}

// 256-bit wide
#ifdef __AVX__

INLINE float8 set_256(float v) { return _mm256_set1_ps(v); }
INLINE float8 set_256(float x0, float y0, float z0, float w0, float x1, float y1, float z1, float w1) { return _mm256_set_ps(w1, z1, y1, x1, w0, z0, y0, x0); }

INLINE float8 zero_256() { return _mm256_setzero_ps(); }

/*
INLINE bool all_zero(const float8 v) {
	const auto res = _mm256_cmpeq_ps(v, _mm_setzero_ps());
	const auto mask = _mm_movemask_ps(res);
	return mask == 0x0000000F;
}
*/

INLINE float8 load_256(const float* addr) { return _mm256_loadu_ps(addr); }
INLINE float8 load_aligned_256(const float* addr) {
	ASSERT(IS_ALIGNED(addr, 16));
	return _mm256_load_ps(addr);
}

INLINE void store(float* addr, float8 v) { _mm256_storeu_ps(addr, v); }
INLINE void store_aligned(float* addr, float8 v) {
	ASSERT(IS_ALIGNED(addr, 16));
	_mm256_store_ps(addr, v);
}

INLINE float8 add(const float8 a, const float8 b) { return _mm256_add_ps(a, b); }
INLINE float8 sub(const float8 a, const float8 b) { return _mm256_sub_ps(a, b); }
INLINE float8 mul(const float8 a, const float8 b) { return _mm256_mul_ps(a, b); }
INLINE float8 div(const float8 a, const float8 b) { return _mm256_div_ps(a, b); }

INLINE float8 bit_and(const float8 a, const float8 b) { return _mm256_and_ps(a, b); }
INLINE float8 bit_and_not(const float8 a, const float8 b) { return _mm256_andnot_ps(a, b); }
INLINE float8 bit_or(const float8 a, const float8 b) { return _mm256_or_ps(a, b); }
INLINE float8 bit_xor(const float8 a, const float8 b) { return _mm256_xor_ps(a, b); }

INLINE float8 cmp_gt(const float8 a, const float8 b) { return _mm256_cmp_ps(a, b, _CMP_GT_OQ); }
INLINE float8 cmp_ge(const float8 a, const float8 b) { return _mm256_cmp_ps(a, b, _CMP_GE_OQ); }
INLINE float8 cmp_lt(const float8 a, const float8 b) { return _mm256_cmp_ps(a, b, _CMP_LT_OQ); }
INLINE float8 cmp_le(const float8 a, const float8 b) { return _mm256_cmp_ps(a, b, _CMP_LE_OQ); }
INLINE float8 cmp_eq(const float8 a, const float8 b) { return _mm256_cmp_ps(a, b, _CMP_EQ_OQ); }
INLINE float8 cmp_neq(const float8 a, const float8 b) { return _mm256_cmp_ps(a, b, _CMP_NEQ_OQ); }

INLINE float8 abs(float8 a) { return bit_and(a, _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF))); }

// @NOTE: If 0.0f is given as input, it will be mapped to 1.0f
INLINE float8 sign(float8 x) {
	const float8 sgn = bit_and(x, set_256(-0.0f));
	const float8 res = bit_xor(sgn, set_256(1.0f));
	return res;
}

INLINE float8 min(const float8 a, const float8 b) {
	const float8 res = _mm256_min_ps(a, b);
	return res;
}

INLINE float8 max(const float8 a, const float8 b) {
	const float8 res = _mm256_max_ps(a, b);
	return res;
}

INLINE float horizontal_min(const float8 x) {
	const float4 lo = _mm256_castps256_ps128(x);
	const float4 hi = _mm256_extractf128_ps(x, 0x1);
	const float4 min_val = min(lo, hi);
	const float result = horizontal_min(min_val);
	return result;
}

INLINE float horizontal_max(const float8 x) {
	const float4 lo = _mm256_castps256_ps128(x);
	const float4 hi = _mm256_extractf128_ps(x, 0x1);
	const float4 max_val = max(lo, hi);
	const float result = horizontal_max(max_val);
	return result;
}

INLINE float8 step(const float8 edge, const float8 x) {
	const float8 cmp = cmp_ge(x, edge);
	const float8 res = bit_and(cmp, set_256(1.f));
	return res;
}

INLINE float8 lerp(const float8 a, const float8 b, float t) {
	const float8 one_minus_alpha = set_256(1.0f - t);
	const float8 alpha = set_256(t);
	const float8 res = add(mul(a, one_minus_alpha), mul(b, alpha));
	return res;
}

INLINE float8 mix(const float8 a, const float8 b, float t) {
	return lerp(a, b, t);
}

INLINE float8 cubic_spline(const __m256 p0, const __m256 p1, const __m256 p2, const __m256 p3, float s, float tension = 0.5f) {
	const float8 vt = set_256(tension);
	const float8 s1 = set_256(s);
	const float8 s2 = mul(s1, s1);
	const float8 s3 = mul(s2, s1);
	const float8 v0 = mul(sub(p2, p0), vt);
	const float8 v1 = mul(sub(p3, p1), vt);
	const float8 x0 = add(mul(set_256(2), sub(p1, p2)), add(v0, v1));
	const float8 x1 = sub(mul(set_256(3), sub(p2, p1)), add(mul(set_256(2), v0), v1));
	const float8 r0 = add(mul(x0, s3), mul(x1, s2));
	const float8 r1 = add(mul(v0, s1), p1);
	const float8 res = add(r0, r1);
	return res;
}

#endif

// 512-bit wide

#if 0

INLINE float16 set_512(float v) { return _mm512_set1_ps(v); }
INLINE float16 set_512(float x0, float y0, float z0, float w0,
					   float x1, float y1, float z1, float w1,
					   float x2, float y2, float z2, float w2,
					   float x3, float y3, float z3, float w3) { return _mm512_set_ps(w3, z3, y3, x3, w2, z2, y2, x2, w1, z1, y1, x1, w0, z0, y0, x0); }

INLINE float16 zero_512() { return _mm512_setzero_ps(); }

/*
INLINE bool all_zero(const float16 v) {
	const auto res = _mm512_cmpeq_ps(v, _mm_setzero_ps());
	const auto mask = _mm_movemask_ps(res);
	return mask == 0x0000000F;
}
*/

INLINE float16 load512(const float* addr) { return _mm512_loadu_ps(addr); }
INLINE float16 load_aligned512(const float* addr) {
	ASSERT(IS_ALIGNED(addr, 16));
	return _mm512_load_ps(addr);
}

INLINE void store(float* addr, float16 v) { _mm512_storeu_ps(addr, v); }
INLINE void store_aligned(float* addr, float16 v) {
	ASSERT(IS_ALIGNED(addr, 16));
	_mm512_store_ps(addr, v);
}

INLINE float16 add(const float16 a, const float16 b) { return _mm512_add_ps(a, b); }
INLINE float16 sub(const float16 a, const float16 b) { return _mm512_sub_ps(a, b); }
INLINE float16 mul(const float16 a, const float16 b) { return _mm512_mul_ps(a, b); }
INLINE float16 div(const float16 a, const float16 b) { return _mm512_div_ps(a, b); }

INLINE float16 bit_and(const float16 a, const float16 b) { return _mm512_and_ps(a, b); }
INLINE float16 bit_or(const float16 a, const float16 b) { return _mm512_or_ps(a, b); }
INLINE float16 bit_xor(const float16 a, const float16 b) { return _mm512_xor_ps(a, b); }

INLINE mask16 cmp_gt(const float16 a, const float16 b) { return _mm512_cmp_ps_mask(a, b, _CMP_GT_OQ); }
INLINE mask16 cmp_ge(const float16 a, const float16 b) { return _mm512_cmp_ps_mask(a, b, _CMP_GE_OQ); }
INLINE mask16 cmp_lt(const float16 a, const float16 b) { return _mm512_cmp_ps_mask(a, b, _CMP_LT_OQ); }
INLINE mask16 cmp_le(const float16 a, const float16 b) { return _mm512_cmp_ps_mask(a, b, _CMP_LE_OQ); }
INLINE mask16 cmp_eq(const float16 a, const float16 b) { return _mm512_cmp_ps_mask(a, b, _CMP_EQ_OQ); }

INLINE float16 abs(float16 a) { return bit_and(a, _mm512_castsi512_ps(_mm512_set1_epi32(0x7FFFFFFF))); }
INLINE float16 sign(float16 x) {
	const mask16 lz = cmp_lt(x, zero_512());
	const float16 res = _mm512_mask_mov_ps(set_512(1.0f), lz, set_512(-1.0f));
	return res;
}

INLINE float16 step(const float16 edge, const float16 x) {
	const mask16 ge = cmp_ge(x, edge);
	const float16 res = _mm512_mask_mov_ps(zero_512(), ge, set_512(1.f));
	return res;
}

INLINE float16 lerp(const float16 a, const float16 b, float t) {
	const float16 one_minus_alpha = set_512(1.0f - t);
	const float16 alpha = set_512(t);
	const float16 res = add(mul(a, one_minus_alpha), mul(b, alpha));
	return res;
};

INLINE float16 cubic_spline(const __m512 p0, const __m512 p1, const __m512 p2, const __m512 p3, float s, float tension = 0.5f) {
	const float16 vt = set_512(tension);
	const float16 s1 = set_512(s);
	const float16 s2 = mul(s1, s1);
	const float16 s3 = mul(s2, s1);
	const float16 v0 = mul(sub(p2, p0), vt);
	const float16 v1 = mul(sub(p3, p1), vt);
	const float16 x0 = add(mul(set_512(2), sub(p1, p2)), add(v0, v1));
	const float16 x1 = sub(mul(set_512(3), sub(p2, p1)), add(mul(set_512(2), v0), v1));
	const float16 r0 = add(mul(x0, s3), mul(x1, s2));
	const float16 r1 = add(mul(v0, s1), p1);
	const float16 res = add(r0, r1);
	return res;
}

#endif

}  // namespace simd
