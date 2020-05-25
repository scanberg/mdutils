#pragma once

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

#ifdef __AVX__
#define SIMD_WIDTH 8
#define SIMD_TYPE_F __m256
#define SIMD_LOAD_F simd::load_f256
#define SIMD_SET_F simd::set_f256
#define SIMD_ZERO_F simd::zero_f256()
#define SIMD_TYPE_I __m256i
#define SUMD_LOAD_I simd::load_i256
#define SIMD_SET_I simd::set_i256
#define SIMD_ZERO_I simd::zero_i256
#else
#define SIMD_WIDTH 4
#define SIMD_TYPE_F __m128
#define SIMD_LOAD_F simd::load_f128
#define SIMD_SET_F simd::set_f128
#define SIMD_ZERO_F simd::zero_f128()
#define SIMD_TYPE_I __m128i
#define SUMD_LOAD_I simd::load_i128
#define SIMD_SET_I simd::set_i128
#define SIMD_ZERO_I simd::zero_i128
#endif
#define SIMD_STORE simd::store
#define SIMD_STORE_ALIGNED simd::store_aligned

#define INLINE inline

#ifdef min
#undef min
#endif

// @TODO: Add more functions for integer types

namespace simd {

typedef __m128 float128;
typedef __m128i int128;

#ifdef __AVX__
typedef __m256 float256;
typedef __m256i int256;
#endif
#if 0
typedef __m512 float512;
typedef __m512i int512;
typedef __mmask16 mask512;
#endif

INLINE float128 set_f128(float v) { return _mm_set1_ps(v); }
INLINE float128 set_f128(float x, float y, float z, float w) { return _mm_set_ps(w, z, y, x); }

INLINE int128 set_i128(int v) { return _mm_set1_epi32(v); }
INLINE int128 set_i128(int x, int y, int z, int w) { return _mm_set_epi32(w, z, y, x); }

INLINE float128 zero_f128() { return _mm_setzero_ps(); }
INLINE int128 zero_i128() { return _mm_setzero_si128(); }

INLINE bool all_zero(const float128 v) {
    const auto cmp = _mm_cmpeq_ps(v, _mm_setzero_ps());
    const auto mask = _mm_movemask_ps(cmp);
    return mask == 0x0000000F;
}

INLINE bool all_zero(const int128 v) {
	const auto cmp = _mm_cmpeq_epi32(v, _mm_setzero_si128());
	const auto mask = _mm_movemask_epi8(cmp);
	return mask == 0xFFFF;
}

INLINE float128 load_f128(const float* addr) { return _mm_loadu_ps(addr); }
INLINE float128 load_aligned_f128(const float* addr) {
    ASSERT(IS_ALIGNED(addr, 16));
    return _mm_load_ps(addr);
}

INLINE void store(float* addr, float128 v) { _mm_storeu_ps(addr, v); }
INLINE void store_aligned(float* addr, float128 v) {
    ASSERT(IS_ALIGNED(addr, 16));
    _mm_store_ps(addr, v);
}

INLINE float128 add(float128 a, float128 b) { return _mm_add_ps(a, b); }
INLINE float128 sub(float128 a, float128 b) { return _mm_sub_ps(a, b); }
INLINE float128 mul(float128 a, float128 b) { return _mm_mul_ps(a, b); }
INLINE float128 div(float128 a, float128 b) { return _mm_div_ps(a, b); }

INLINE int128 add(int128 a, int128 b) { return _mm_add_epi32(a, b); }
INLINE int128 sub(int128 a, int128 b) { return _mm_sub_epi32(a, b); }
INLINE int128 mul(int128 a, int128 b) { return _mm_mul_epi32(a, b); }
//INLINE int128 div(int128 a, int128 b) { return _mm_div_epi32(a, b); }

INLINE float128 bit_and(float128 a, float128 b) { return _mm_and_ps(a, b); }
INLINE float128 bit_and_not(float128 a, float128 b) { return _mm_andnot_ps(a, b); }
INLINE float128 bit_or(float128 a, float128 b) { return _mm_or_ps(a, b); }
INLINE float128 bit_xor(float128 a, float128 b) { return _mm_xor_ps(a, b); }

INLINE int128 bit_and(int128 a, int128 b) { return _mm_and_si128(a, b); }
INLINE int128 bit_and_not(int128 a, int128 b) { return _mm_andnot_si128(a, b); }
INLINE int128 bit_or(int128 a, int128 b) { return _mm_or_si128(a, b); }
INLINE int128 bit_xor(int128 a, int128 b) { return _mm_xor_si128(a, b); }

INLINE float128 cmp_gt(float128 a, float128 b) { return _mm_cmpgt_ps(a, b); }
INLINE float128 cmp_ge(float128 a, float128 b) { return _mm_cmpge_ps(a, b); }
INLINE float128 cmp_lt(float128 a, float128 b) { return _mm_cmplt_ps(a, b); }
INLINE float128 cmp_le(float128 a, float128 b) { return _mm_cmple_ps(a, b); }
INLINE float128 cmp_eq(float128 a, float128 b) { return _mm_cmpeq_ps(a, b); }
INLINE float128 cmp_neq(float128 a, float128 b) { return _mm_cmpneq_ps(a, b); }

INLINE float128 abs(float128 a) { return _mm_and_ps(a, _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF))); }
INLINE int128 abs(int128 a) { return _mm_abs_epi32(a); }

// @NOTE: If 0.0f is given as input, it will be mapped to 1.0f
INLINE float128 sign(float128 x) {
	const float128 sgn = bit_and(x, set_f128(-0.0f));
	const float128 res = bit_xor(sgn, set_f128(1.0f));
	return res;
}

INLINE float128 min(float128 a, float128 b) {
	const float128 res = _mm_min_ps(a, b);
	return res;
}

INLINE float128 max(float128 a, float128 b) {
	const float128 res = _mm_max_ps(a, b);
	return res;
}

INLINE float horizontal_min(float128 x) {
	const float128 min1 = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 0, 3, 2));
	const float128 min2 = _mm_min_ps(x, min1);
	const float128 min3 = _mm_shuffle_ps(min2, min2, _MM_SHUFFLE(0, 0, 0, 1));
	const float128 min4 = _mm_min_ps(min2, min3);
	float result = _mm_cvtss_f32(min4);
	return result;
}

INLINE float horizontal_max(float128 x) {
	const float128 max1 = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 0, 3, 2));
	const float128 max2 = _mm_max_ps(x, max1);
	const float128 max3 = _mm_shuffle_ps(max2, max2, _MM_SHUFFLE(0, 0, 0, 1));
	const float128 max4 = _mm_max_ps(max2, max3);
	float result = _mm_cvtss_f32(max4);
	return result;
}

// https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86
INLINE float horizontal_add(float128 x) {
	float128 shuf = _mm_shuffle_ps(x, x, _MM_SHUFFLE(2, 3, 0, 1));
	float128 sums = _mm_add_ps(x, shuf);
	shuf = _mm_movehl_ps(shuf, sums);
	sums = _mm_add_ss(sums, shuf);
	return _mm_cvtss_f32(sums);
}

INLINE float128 step(float128 edge, float128 x) {
    const float128 cmp = cmp_ge(x, edge);
    const float128 res = bit_and(cmp, set_f128(1.f));
    return res;
}

INLINE float128 lerp(float128 a, float128 b, float t) {
    const float128 one_minus_alpha = set_f128(1.0f - t);
    const float128 alpha = set_f128(t);
    const float128 res = add(mul(a, one_minus_alpha), mul(b, alpha));
    return res;
}

INLINE float128 mix(float128 a, float128 b, float t) {
	return lerp(a, b, t);
}

INLINE float128 cubic_spline(__m128 p0, __m128 p1, __m128 p2, __m128 p3, float s, float tension = 0.5f) {
    const float128 vt = set_f128(tension);
    const float128 s1 = set_f128(s);
    const float128 s2 = mul(s1, s1);
    const float128 s3 = mul(s2, s1);
    const float128 v0 = mul(sub(p2, p0), vt);
    const float128 v1 = mul(sub(p3, p1), vt);
    const float128 x0 = add(mul(set_f128(2), sub(p1, p2)), add(v0, v1));
    const float128 x1 = sub(mul(set_f128(3), sub(p2, p1)), add(mul(set_f128(2), v0), v1));
    const float128 r0 = add(mul(x0, s3), mul(x1, s2));
    const float128 r1 = add(mul(v0, s1), p1);
    const float128 res = add(r0, r1);
    return res;
}

// 256-bit wide
#ifdef __AVX__

INLINE float256 set_f256(float v) { return _mm256_set1_ps(v); }
INLINE float256 set_f256(float x0, float y0, float z0, float w0, float x1, float y1, float z1, float w1) { return _mm256_set_ps(w1, z1, y1, x1, w0, z0, y0, x0); }

INLINE float256 zero_f256() { return _mm256_setzero_ps(); }

/*
INLINE bool all_zero(float8 v) {
	const auto res = _mm256_cmpeq_ps(v, _mm_setzero_ps());
	const auto mask = _mm_movemask_ps(res);
	return mask == 0x0000000F;
}
*/

INLINE float256 load_f256(const float* addr) { return _mm256_loadu_ps(addr); }
INLINE float256 load_aligned_f256(const float* addr) {
	ASSERT(IS_ALIGNED(addr, 16));
	return _mm256_load_ps(addr);
}

INLINE void store(float* addr, float256 v) { _mm256_storeu_ps(addr, v); }
INLINE void store_aligned(float* addr, float256 v) {
	ASSERT(IS_ALIGNED(addr, 16));
	_mm256_store_ps(addr, v);
}

INLINE float256 add(float256 a, float256 b) { return _mm256_add_ps(a, b); }
INLINE float256 sub(float256 a, float256 b) { return _mm256_sub_ps(a, b); }
INLINE float256 mul(float256 a, float256 b) { return _mm256_mul_ps(a, b); }
INLINE float256 div(float256 a, float256 b) { return _mm256_div_ps(a, b); }

INLINE float256 bit_and(float256 a, float256 b) { return _mm256_and_ps(a, b); }
INLINE float256 bit_and_not(float256 a, float256 b) { return _mm256_andnot_ps(a, b); }
INLINE float256 bit_or(float256 a, float256 b) { return _mm256_or_ps(a, b); }
INLINE float256 bit_xor(float256 a, float256 b) { return _mm256_xor_ps(a, b); }

INLINE float256 cmp_gt(float256 a, float256 b) { return _mm256_cmp_ps(a, b, _CMP_GT_OQ); }
INLINE float256 cmp_ge(float256 a, float256 b) { return _mm256_cmp_ps(a, b, _CMP_GE_OQ); }
INLINE float256 cmp_lt(float256 a, float256 b) { return _mm256_cmp_ps(a, b, _CMP_LT_OQ); }
INLINE float256 cmp_le(float256 a, float256 b) { return _mm256_cmp_ps(a, b, _CMP_LE_OQ); }
INLINE float256 cmp_eq(float256 a, float256 b) { return _mm256_cmp_ps(a, b, _CMP_EQ_OQ); }
INLINE float256 cmp_neq(float256 a, float256 b) { return _mm256_cmp_ps(a, b, _CMP_NEQ_OQ); }

INLINE float256 abs(float256 a) { return bit_and(a, _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF))); }

// @NOTE: If 0.0f is given as input, it will be mapped to 1.0f
INLINE float256 sign(float256 x) {
	const float256 sgn = bit_and(x, set_f256(-0.0f));
	const float256 res = bit_xor(sgn, set_f256(1.0f));
	return res;
}

INLINE float256 min(float256 a, float256 b) {
	const float256 res = _mm256_min_ps(a, b);
	return res;
}

INLINE float256 max(float256 a, float256 b) {
	const float256 res = _mm256_max_ps(a, b);
	return res;
}

INLINE float horizontal_min(float256 x) {
	const float128 lo = _mm256_castps256_ps128(x);
	const float128 hi = _mm256_extractf128_ps(x, 0x1);
	const float128 min_val = min(lo, hi);
	const float result = horizontal_min(min_val);
	return result;
}

INLINE float horizontal_max(float256 x) {
	const float128 lo = _mm256_castps256_ps128(x);
	const float128 hi = _mm256_extractf128_ps(x, 0x1);
	const float128 max_val = max(lo, hi);
	const float result = horizontal_max(max_val);
	return result;
}

// https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86
INLINE float horizontal_add(float256 v) {
	float128 vlow = _mm256_castps256_ps128(v);
	const float128 vhigh = _mm256_extractf128_ps(v, 1); // high 128
	vlow = _mm_add_ps(vlow, vhigh);     // add the low 128
	float128 shuf = _mm_movehdup_ps(vlow);        // broadcast elements 3,1 to 2,0
	float128 sums = _mm_add_ps(vlow, shuf);
	shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
	sums = _mm_add_ss(sums, shuf);
	return _mm_cvtss_f32(sums);
}

INLINE float256 step(float256 edge, float256 x) {
	const float256 cmp = cmp_ge(x, edge);
	const float256 res = bit_and(cmp, set_f256(1.f));
	return res;
}

INLINE float256 lerp(float256 a, float256 b, float t) {
	const float256 one_minus_alpha = set_f256(1.0f - t);
	const float256 alpha = set_f256(t);
	const float256 res = add(mul(a, one_minus_alpha), mul(b, alpha));
	return res;
}

INLINE float256 mix(float256 a, float256 b, float t) {
	return lerp(a, b, t);
}

INLINE float256 cubic_spline(__m256 p0, __m256 p1, __m256 p2, __m256 p3, float s, float tension = 0.5f) {
	const float256 vt = set_f256(tension);
	const float256 s1 = set_f256(s);
	const float256 s2 = mul(s1, s1);
	const float256 s3 = mul(s2, s1);
	const float256 v0 = mul(sub(p2, p0), vt);
	const float256 v1 = mul(sub(p3, p1), vt);
	const float256 x0 = add(mul(set_f256(2), sub(p1, p2)), add(v0, v1));
	const float256 x1 = sub(mul(set_f256(3), sub(p2, p1)), add(mul(set_f256(2), v0), v1));
	const float256 r0 = add(mul(x0, s3), mul(x1, s2));
	const float256 r1 = add(mul(v0, s1), p1);
	const float256 res = add(r0, r1);
	return res;
}

#endif

}  // namespace simd
