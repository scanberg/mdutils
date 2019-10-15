#include "molecule_utils.h"

#include <core/common.h>
#include <core/simd.h>
#include <core/hash.h>
#include <core/log.h>
#include <core/spatial_hash.h>

#include <mol/trajectory_utils.h>

#include <svd3/svd3.h>
#include <ctype.h>

inline __m128 apply_pbc(const __m128 x, const __m128 box_ext) {
    const __m128 add = simd::bit_and(simd::cmp_lt(x, simd::zero_f128()), box_ext);
    const __m128 sub = simd::bit_and(simd::cmp_gt(x, box_ext), box_ext);
    const __m128 res = simd::bit_and(x, simd::sub(add, sub));
    return res;
}

template <typename T>
inline T de_periodize(T a, T b, T full_ext, T half_ext) {
    const T delta = b - a;
    const T signed_mask = math::sign(delta) * math::step(half_ext, math::abs(delta));
    const T res = b - full_ext * signed_mask;
    return res;
}

inline __m128 de_periodize(const __m128 a, const __m128 b, const __m128 full_ext, const __m128 half_ext) {
    const __m128 delta = simd::sub(b, a);
    const __m128 abs_delta = simd::abs(delta);
    const __m128 sign_delta = simd::sign(delta);
    const __m128 step_delta = simd::step(half_ext, abs_delta);
    const __m128 signed_mask = simd::mul(sign_delta, step_delta);
    const __m128 res = simd::sub(b, simd::mul(full_ext, signed_mask));
    return res;
}

#ifdef __AVX__
inline __m256 de_periodize(const __m256 a, const __m256 b, const __m256 full_ext, const __m256 half_ext) {
    const __m256 delta = simd::sub(b, a);
    const __m256 signed_mask = simd::mul(simd::sign(delta), simd::step(half_ext, simd::abs(delta)));
    const __m256 res = simd::sub(b, simd::mul(full_ext, signed_mask));
    return res;
}
#endif

#ifdef __AVX512__
inline __m512 de_periodize(const __m512 a, const __m512 b, const __m512 full_ext, const __m512 half_ext) {
    const __m512 delta = simd::sub(b, a);
    const __m512 signed_mask = simd::mul(simd::sign(delta), simd::step(half_ext, simd::abs(delta)));
    const __m512 res = simd::sub(b, simd::mul(full_ext, signed_mask));
    return res;
}
#endif

void translate(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, int64 count, const vec3& translation) {
    int64 i = 0;

    const int64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
    if (simd_count > 0) {
        SIMD_TYPE_F t_x = SIMD_SET_F(translation.x);
        SIMD_TYPE_F t_y = SIMD_SET_F(translation.y);
        SIMD_TYPE_F t_z = SIMD_SET_F(translation.z);

        for (; i < simd_count; i += SIMD_WIDTH) {
            SIMD_TYPE_F p_x = SIMD_LOAD_F(in_out_x + i);
            SIMD_TYPE_F p_y = SIMD_LOAD_F(in_out_y + i);
            SIMD_TYPE_F p_z = SIMD_LOAD_F(in_out_z + i);

            p_x = simd::add(p_x, t_x);
            p_y = simd::add(p_y, t_y);
            p_z = simd::add(p_z, t_z);

            SIMD_STORE(in_out_x + i, p_x);
            SIMD_STORE(in_out_y + i, p_y);
            SIMD_STORE(in_out_z + i, p_z);
        }
    }

    for (; i < count; i++) {
        in_out_x[i] += translation.x;
        in_out_y[i] += translation.y;
        in_out_z[i] += translation.z;
    }
}

void translate_ref(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, int64 count, const vec3& translation) {
    for (int64 i = 0; i < count; i++) {
        in_out_x[i] += translation.x;
        in_out_y[i] += translation.y;
        in_out_z[i] += translation.z;
    }
}

void transform_ref(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, int64 count, const mat4& transformation, float w_comp) {
    for (int64 i = 0; i < count; i++) {
        vec4 v = {in_out_x[i], in_out_y[i], in_out_z[i], w_comp};
        v = transformation * v;
        in_out_x[i] = v.x;
        in_out_y[i] = v.y;
        in_out_z[i] = v.z;
    }
}

void transform(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, int64 count, const mat4& transformation, float w_comp) {
    const SIMD_TYPE_F m11 = SIMD_SET_F(transformation[0][0]);
    const SIMD_TYPE_F m12 = SIMD_SET_F(transformation[0][1]);
    const SIMD_TYPE_F m13 = SIMD_SET_F(transformation[0][2]);

    const SIMD_TYPE_F m21 = SIMD_SET_F(transformation[1][0]);
    const SIMD_TYPE_F m22 = SIMD_SET_F(transformation[1][1]);
    const SIMD_TYPE_F m23 = SIMD_SET_F(transformation[1][2]);

    const SIMD_TYPE_F m31 = SIMD_SET_F(transformation[2][0]);
    const SIMD_TYPE_F m32 = SIMD_SET_F(transformation[2][1]);
    const SIMD_TYPE_F m33 = SIMD_SET_F(transformation[2][2]);

    const SIMD_TYPE_F m41 = SIMD_SET_F(transformation[3][0]);
    const SIMD_TYPE_F m42 = SIMD_SET_F(transformation[3][1]);
    const SIMD_TYPE_F m43 = SIMD_SET_F(transformation[3][2]);

    const SIMD_TYPE_F w = SIMD_SET_F(w_comp);

    int64 i = 0;
    const int64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
    for (; i < simd_count; i += SIMD_WIDTH) {
        const SIMD_TYPE_F x = SIMD_LOAD_F(in_out_x + i);
        const SIMD_TYPE_F y = SIMD_LOAD_F(in_out_y + i);
        const SIMD_TYPE_F z = SIMD_LOAD_F(in_out_z + i);

        const SIMD_TYPE_F m11x = simd::mul(m11, x);
        const SIMD_TYPE_F m21y = simd::mul(m21, y);
        const SIMD_TYPE_F m31z = simd::mul(m31, z);
        const SIMD_TYPE_F m41w = simd::mul(m41, w);

        const SIMD_TYPE_F m12x = simd::mul(m12, x);
        const SIMD_TYPE_F m22y = simd::mul(m22, y);
        const SIMD_TYPE_F m32z = simd::mul(m32, z);
        const SIMD_TYPE_F m42w = simd::mul(m42, w);

        const SIMD_TYPE_F m13x = simd::mul(m13, x);
        const SIMD_TYPE_F m23y = simd::mul(m23, y);
        const SIMD_TYPE_F m33z = simd::mul(m33, z);
        const SIMD_TYPE_F m43w = simd::mul(m43, w);

        const SIMD_TYPE_F res_x = simd::add(simd::add(m11x, m21y), simd::add(m31z, m41w));
        const SIMD_TYPE_F res_y = simd::add(simd::add(m12x, m22y), simd::add(m32z, m42w));
        const SIMD_TYPE_F res_z = simd::add(simd::add(m13x, m23y), simd::add(m33z, m43w));

        SIMD_STORE(in_out_x + i, res_x);
        SIMD_STORE(in_out_y + i, res_y);
        SIMD_STORE(in_out_z + i, res_z);
    }

    for (; i < count; i++) {
        const float x = in_out_x[i];
        const float y = in_out_y[i];
        const float z = in_out_z[i];

        in_out_x[i] = x * transformation[0][0] + y * transformation[1][0] + z * transformation[2][0] + w_comp * transformation[3][0];
        in_out_y[i] = x * transformation[0][1] + y * transformation[1][1] + z * transformation[2][1] + w_comp * transformation[3][1];
        in_out_z[i] = x * transformation[0][2] + y * transformation[1][2] + z * transformation[2][2] + w_comp * transformation[3][2];
    }
}

void homogeneous_transform(float* RESTRICT pos_x, float* RESTRICT pos_y, float* RESTRICT pos_z, int64 count, const mat4& transformation) {
    for (int64 i = 0; i < count; i++) {
        const vec4 p = transformation * vec4(pos_x[i], pos_y[i], pos_z[i], 1.0f);
        pos_x[i] = p.x / p.w;
        pos_y[i] = p.y / p.w;
        pos_z[i] = p.z / p.w;
    }
}

AABB compute_aabb(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, int64 count) {
    if (count == 0) {
        return {};
    }

    AABB aabb;

    int64 i = 0;
    if (count > SIMD_WIDTH) {  // @NOTE: There is probably some number where this makes most sense
        SIMD_TYPE_F min_x = SIMD_LOAD_F(in_x);
        SIMD_TYPE_F min_y = SIMD_LOAD_F(in_y);
        SIMD_TYPE_F min_z = SIMD_LOAD_F(in_z);

        SIMD_TYPE_F max_x = min_x;
        SIMD_TYPE_F max_y = min_y;
        SIMD_TYPE_F max_z = min_z;

        i += SIMD_WIDTH;
        const int64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
        for (; i < simd_count; i += SIMD_WIDTH) {
            const SIMD_TYPE_F x = SIMD_LOAD_F(in_x + i);
            const SIMD_TYPE_F y = SIMD_LOAD_F(in_y + i);
            const SIMD_TYPE_F z = SIMD_LOAD_F(in_z + i);

            min_x = simd::min(min_x, x);
            min_y = simd::min(min_y, y);
            min_z = simd::min(min_z, z);

            max_x = simd::max(max_x, x);
            max_y = simd::max(max_y, y);
            max_z = simd::max(max_z, z);
        }

        aabb.min = {simd::horizontal_min(min_x), simd::horizontal_min(min_y), simd::horizontal_min(min_z)};
        aabb.max = {simd::horizontal_max(max_x), simd::horizontal_max(max_y), simd::horizontal_max(max_z)};
    }

    for (; i < count; i++) {
        const vec3 p = {in_x[i], in_y[i], in_z[i]};
        aabb.min = math::min(aabb.min, p);
        aabb.max = math::max(aabb.max, p);
    }

    return aabb;
}

AABB compute_aabb(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, const float* in_r, int64 count) {
    if (count == 0) {
        return {};
    }

    AABB aabb;

    int64 i = 0;
    if (count > SIMD_WIDTH) {  // @NOTE: There is probably some number where this makes most sense
        SIMD_TYPE_F x = SIMD_LOAD_F(in_x);
        SIMD_TYPE_F y = SIMD_LOAD_F(in_y);
        SIMD_TYPE_F z = SIMD_LOAD_F(in_z);
        SIMD_TYPE_F r = SIMD_LOAD_F(in_r);

        SIMD_TYPE_F min_x = simd::sub(x, r);
        SIMD_TYPE_F min_y = simd::sub(y, r);
        SIMD_TYPE_F min_z = simd::sub(z, r);

        SIMD_TYPE_F max_x = simd::add(x, r);
        SIMD_TYPE_F max_y = simd::add(y, r);
        SIMD_TYPE_F max_z = simd::add(z, r);

        i += SIMD_WIDTH;
        const int64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
        for (; i < simd_count; i += SIMD_WIDTH) {
            x = SIMD_LOAD_F(in_x + i);
            y = SIMD_LOAD_F(in_y + i);
            z = SIMD_LOAD_F(in_z + i);
            r = SIMD_LOAD_F(in_r + i);

            min_x = simd::min(min_x, simd::sub(x, r));
            min_y = simd::min(min_y, simd::sub(y, r));
            min_z = simd::min(min_z, simd::sub(z, r));

            max_x = simd::max(max_x, simd::add(x, r));
            max_y = simd::max(max_y, simd::add(y, r));
            max_z = simd::max(max_z, simd::add(z, r));
        }

        aabb.min = {simd::horizontal_min(min_x), simd::horizontal_min(min_y), simd::horizontal_min(min_z)};
        aabb.max = {simd::horizontal_max(max_x), simd::horizontal_max(max_y), simd::horizontal_max(max_z)};
    }

    for (; i < count; i++) {
        const vec3 p = {in_x[i], in_y[i], in_z[i]};
        aabb.min = math::min(aabb.min, p);
        aabb.max = math::max(aabb.max, p);
    }

    return aabb;
}

vec3 compute_com(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, int64 count) {
    if (count == 0) return vec3(0);
    if (count == 1) return {in_x[0], in_y[0], in_z[0]};

    vec3 sum{0};
    int64 i = 0;

    const int64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
    if (simd_count > SIMD_WIDTH) {
        SIMD_TYPE_F x = SIMD_LOAD_F(in_x);
        SIMD_TYPE_F y = SIMD_LOAD_F(in_y);
        SIMD_TYPE_F z = SIMD_LOAD_F(in_z);
        i += SIMD_WIDTH;
        for (; i < simd_count; i += SIMD_WIDTH) {
            x = simd::add(x, SIMD_LOAD_F(in_x + i));
            y = simd::add(y, SIMD_LOAD_F(in_y + i));
            z = simd::add(z, SIMD_LOAD_F(in_z + i));
        }
        sum.x = simd::horizontal_add(x);
        sum.y = simd::horizontal_add(y);
        sum.z = simd::horizontal_add(z);
    }

    for (; i < count; i++) {
        sum.x += in_x[i];
        sum.y += in_y[i];
        sum.z += in_z[i];
    }
    return sum / (float)count;
}

vec3 compute_com(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, const float* RESTRICT in_m, int64 count) {
    if (count == 0) return vec3(0);
    if (count == 1) return {in_x[0], in_y[0], in_z[0]};

    vec3 vec_sum{0, 0, 0};
    float mass_sum = 0.0f;
    int64 i = 0;

    const int64 simd_count = (count / SIMD_WIDTH) * SIMD_WIDTH;
    if (simd_count > SIMD_WIDTH) {
        SIMD_TYPE_F m = SIMD_LOAD_F(in_m);
        SIMD_TYPE_F x = simd::mul(SIMD_LOAD_F(in_x), m);
        SIMD_TYPE_F y = simd::mul(SIMD_LOAD_F(in_y), m);
        SIMD_TYPE_F z = simd::mul(SIMD_LOAD_F(in_z), m);
        i += SIMD_WIDTH;
        for (; i < simd_count; i += SIMD_WIDTH) {
            const SIMD_TYPE_F mass = SIMD_LOAD_F(in_m + i);
            x = simd::add(x, simd::mul(SIMD_LOAD_F(in_x + i), mass));
            y = simd::add(y, simd::mul(SIMD_LOAD_F(in_y + i), mass));
            z = simd::add(z, simd::mul(SIMD_LOAD_F(in_z + i), mass));
            m = simd::add(m, mass);
        }
        vec_sum.x = simd::horizontal_add(x);
        vec_sum.y = simd::horizontal_add(y);
        vec_sum.z = simd::horizontal_add(z);
        mass_sum = simd::horizontal_add(m);
    }

    for (; i < count; i++) {
        const float mass = in_m[i];
        vec_sum.x += in_x[i] * mass;
        vec_sum.y += in_y[i] * mass;
        vec_sum.z += in_z[i] * mass;
        mass_sum += mass;
    }
    return vec_sum / mass_sum;
}

vec3 compute_com_periodic(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, const float* RESTRICT in_m, int64 count, const mat3& box) {
    if (count == 0) return vec3(0);
    if (count == 1) return {in_x[0], in_y[0], in_z[0]};

    const vec3 full_ext = box * vec3(1.0f);
    const vec3 half_ext = box * vec3(0.5f);

    float mass_sum = in_m[0];
    vec3 vec_sum = vec3(in_x[0], in_y[0], in_z[0]) * in_m[0];

    for (int64 i = 1; i < count; i++) {
        const vec3 p = {in_x[i], in_y[i], in_z[i]};
        const float mass = in_m[i];
        const vec3 com = vec_sum / mass_sum;
        const vec3 dp = de_periodize(com, p, full_ext, half_ext);
        vec_sum += dp * mass;
        mass_sum += mass;
    }
    return vec_sum / mass_sum;
}

// @TODO: Finalize implementation
/*
vec3 compute_com_periodic_vectorized(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, const float* RESTRICT in_m, int64 count, const mat3& box) {
    if (count == 0) return vec3(0);
    if (count == 1) return {in_x[0], in_y[0], in_z[0]};

    const vec3 full_ext = box * vec3(1.0f);
    const vec3 half_ext = box * vec3(0.5f);
    return {0, 0, 0};
}
*/

vec3 compute_com(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, const Element* RESTRICT in_element, int64 count) {
    if (count == 0) return {0, 0, 0};
    if (count == 1) return {in_x[0], in_y[0], in_z[0]};

    vec3 v_sum{0};
    float m_sum = 0.0f;
    for (int32 i = 0; i < count; i++) {
        const vec3 v = {in_x[i], in_y[i], in_z[i]};
        const float m = element::atomic_mass(in_element[i]);
        v_sum += v * m;
        m_sum += m;
    }

    return v_sum / m_sum;
}

mat3 compute_covariance_matrix(const float* RESTRICT x, const float* RESTRICT y, const float* RESTRICT z, const float* RESTRICT mass, int64 count, const vec3& com) {
    mat3 A{0};
    float mass_sum = 0.0f;
    for (int64 i = 0; i < count; i++) {
        // @TODO: Vectorize...
        const float qx = x[i] - com.x;
        const float qy = y[i] - com.y;
        const float qz = z[i] - com.z;
        const float m = mass[i];
        mass_sum += m;

        A[0][0] += m * qx * qx;
        A[0][1] += m * qy * qx;
        A[0][2] += m * qz * qx;
        A[1][0] += m * qx * qy;
        A[1][1] += m * qy * qy;
        A[1][2] += m * qz * qy;
        A[2][0] += m * qx * qz;
        A[2][1] += m * qy * qz;
        A[2][2] += m * qz * qz;
    }

    return A / mass_sum;
}

#define ARGS(M) M[0][0], M[1][0], M[2][0], M[0][1], M[1][1], M[2][1], M[0][2], M[1][2], M[2][2]
EigenFrame compute_eigen_frame(const float* RESTRICT in_x, const float* RESTRICT in_y, const float* RESTRICT in_z, const float* RESTRICT in_mass, int64 count, const vec3& com) {
    const mat3 M = compute_covariance_matrix(in_x, in_y, in_z, in_mass, count, com);
    mat3 U, S, V;
    svd(ARGS(M), ARGS(U), ARGS(S), ARGS(V));
    const mat3 Ut = glm::transpose(U);

    const float e_val[] = {math::sqrt(S[0][0]), math::sqrt(S[1][1]), math::sqrt(S[2][2])};
    const vec3 e_vec[] = {Ut[0], Ut[1], Ut[2]};
    int l[3] = {0, 1, 2};

    const auto swap = [](int& x, int& y) {
        int tmp = x;
        x = y;
        y = tmp;
    };

    if (e_val[l[0]] < e_val[l[1]]) swap(l[0], l[1]);
    if (e_val[l[1]] < e_val[l[2]]) swap(l[1], l[2]);
    if (e_val[l[0]] < e_val[l[1]]) swap(l[0], l[1]);

    EigenFrame ef;
    ef.value[0] = e_val[l[0]];
    ef.value[1] = e_val[l[1]];
    ef.value[2] = e_val[l[2]];

    ef.vector[0] = e_vec[l[0]];
    ef.vector[1] = e_vec[l[1]];
    ef.vector[2] = e_vec[l[2]];

    return ef;
}
#undef ARGS

void recenter_trajectory(MoleculeDynamic* dynamic, Bitfield atom_mask) {
    ASSERT(dynamic);
    if (!dynamic->operator bool()) {
        LOG_ERROR("Dynamic is not valid.");
        return;
    }

    const auto& mol = dynamic->molecule;
    auto& traj = dynamic->trajectory;

    int32 count = bitfield::number_of_bits_set(atom_mask);
    void* mem = TMP_MALLOC(count * sizeof(float) * 4);
    defer { TMP_FREE(mem); };
    float* x = (float*)mem + 0 * count;
    float* y = (float*)mem + 1 * count;
    float* z = (float*)mem + 2 * count;
    float* m = (float*)mem + 3 * count;

    bitfield::extract_data_from_mask(m, mol.atom.mass, atom_mask);

    for (auto& frame : traj.frame_buffer) {
        bitfield::extract_data_from_mask(x, frame.atom_position.x, atom_mask);
        bitfield::extract_data_from_mask(y, frame.atom_position.y, atom_mask);
        bitfield::extract_data_from_mask(z, frame.atom_position.z, atom_mask);
        const vec3 com = compute_com_periodic(x, y, z, m, count, frame.box);
        translate(frame.atom_position.x, frame.atom_position.y, frame.atom_position.z, mol.atom.count, -com);
        apply_pbc(frame.atom_position.x, frame.atom_position.y, frame.atom_position.z, mol.atom.mass, mol.sequences, frame.box);
    }
}

// clang-format off
void linear_interpolation_scalar(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
								 const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
                                 const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
								 int64 count, float t)
// clang-format on
{
    for (int64 i = 0; i < count; i++) {
        out_x[i] = in_x0[i] * (1.0f - t) + in_x1[i] * t;
        out_y[i] = in_y0[i] * (1.0f - t) + in_y1[i] * t;
        out_z[i] = in_z0[i] * (1.0f - t) + in_z1[i] * t;
    }
}

// clang-format off
void linear_interpolation(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
						  const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
                          const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
						  int64 count, float t)
// clang-format on
{
    for (int64 i = 0; i < count; i += SIMD_WIDTH) {
        const SIMD_TYPE_F x0 = SIMD_LOAD_F(in_x0 + i);
        const SIMD_TYPE_F y0 = SIMD_LOAD_F(in_y0 + i);
        const SIMD_TYPE_F z0 = SIMD_LOAD_F(in_z0 + i);

        const SIMD_TYPE_F x1 = SIMD_LOAD_F(in_x1 + i);
        const SIMD_TYPE_F y1 = SIMD_LOAD_F(in_y1 + i);
        const SIMD_TYPE_F z1 = SIMD_LOAD_F(in_z1 + i);

        const SIMD_TYPE_F x = simd::lerp(x0, x1, t);
        const SIMD_TYPE_F y = simd::lerp(y0, y1, t);
        const SIMD_TYPE_F z = simd::lerp(z0, z1, t);

        simd::store(out_x + i, x);
        simd::store(out_y + i, y);
        simd::store(out_z + i, z);
    }
}

// clang-format off
void linear_interpolation_128(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
							  const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
                              const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
							  int64 count, float t)
// clang-format on
{
    for (int64 i = 0; i < count; i += 4) {
        const __m128 x0 = simd::load_f128(in_x0 + i);
        const __m128 y0 = simd::load_f128(in_y0 + i);
        const __m128 z0 = simd::load_f128(in_z0 + i);

        const __m128 x1 = simd::load_f128(in_x1 + i);
        const __m128 y1 = simd::load_f128(in_y1 + i);
        const __m128 z1 = simd::load_f128(in_z1 + i);

        const __m128 x = simd::lerp(x0, x1, t);
        const __m128 y = simd::lerp(y0, y1, t);
        const __m128 z = simd::lerp(z0, z1, t);

        simd::store(out_x + i, x);
        simd::store(out_y + i, y);
        simd::store(out_z + i, z);
    }
}

#ifdef __AVX__
// clang-format off
void linear_interpolation_256(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
							  const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
                              const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
							  int64 count, float t)
// clang-format on
{
    for (int64 i = 0; i < count; i += 8) {
        const __m256 x0 = simd::load_f256(in_x0 + i);
        const __m256 y0 = simd::load_f256(in_y0 + i);
        const __m256 z0 = simd::load_f256(in_z0 + i);

        const __m256 x1 = simd::load_f256(in_x1 + i);
        const __m256 y1 = simd::load_f256(in_y1 + i);
        const __m256 z1 = simd::load_f256(in_z1 + i);

        const __m256 x = simd::lerp(x0, x1, t);
        const __m256 y = simd::lerp(y0, y1, t);
        const __m256 z = simd::lerp(z0, z1, t);

        simd::store(out_x + i, x);
        simd::store(out_y + i, y);
        simd::store(out_z + i, z);
    }
}
#endif

// clang-format off
void linear_interpolation_pbc_scalar(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
									 const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
                                     const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
									 int64 count, float t, const mat3& sim_box)
// clang-format on
{

    const float full_ext_x = sim_box[0][0];
    const float full_ext_y = sim_box[1][1];
    const float full_ext_z = sim_box[2][2];

    const float half_ext_x = full_ext_x * 0.5f;
    const float half_ext_y = full_ext_y * 0.5f;
    const float half_ext_z = full_ext_z * 0.5f;

    for (int64 i = 0; i < count; i++) {
        float x0 = in_x0[i];
        float y0 = in_y0[i];
        float z0 = in_z0[i];

        float x1 = in_x1[i];
        float y1 = in_y1[i];
        float z1 = in_z1[i];

        x1 = de_periodize(x0, x1, full_ext_x, half_ext_x);
        y1 = de_periodize(y0, y1, full_ext_y, half_ext_y);
        z1 = de_periodize(z0, z1, full_ext_z, half_ext_z);

        const float x = x0 * (1.0f - t) + x1 * t;
        const float y = y0 * (1.0f - t) + y1 * t;
        const float z = z0 * (1.0f - t) + z1 * t;

        out_x[i] = x;
        out_y[i] = y;
        out_z[i] = z;
    }
}

// clang-format off
void linear_interpolation_pbc(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
							  const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
                              const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
							  int64 count, float t, const mat3& sim_box)
// clang-format on
{
    const SIMD_TYPE_F full_box_ext_x = SIMD_SET_F(sim_box[0][0]);
    const SIMD_TYPE_F full_box_ext_y = SIMD_SET_F(sim_box[1][1]);
    const SIMD_TYPE_F full_box_ext_z = SIMD_SET_F(sim_box[2][2]);

    const SIMD_TYPE_F half_box_ext_x = simd::mul(full_box_ext_x, SIMD_SET_F(0.5f));
    const SIMD_TYPE_F half_box_ext_y = simd::mul(full_box_ext_y, SIMD_SET_F(0.5f));
    const SIMD_TYPE_F half_box_ext_z = simd::mul(full_box_ext_z, SIMD_SET_F(0.5f));

    for (int64 i = 0; i < count; i += SIMD_WIDTH) {
        SIMD_TYPE_F x0 = SIMD_LOAD_F(in_x0 + i);
        SIMD_TYPE_F y0 = SIMD_LOAD_F(in_y0 + i);
        SIMD_TYPE_F z0 = SIMD_LOAD_F(in_z0 + i);

        SIMD_TYPE_F x1 = SIMD_LOAD_F(in_x1 + i);
        SIMD_TYPE_F y1 = SIMD_LOAD_F(in_y1 + i);
        SIMD_TYPE_F z1 = SIMD_LOAD_F(in_z1 + i);

        x1 = de_periodize(x0, x1, full_box_ext_x, half_box_ext_x);
        y1 = de_periodize(y0, y1, full_box_ext_y, half_box_ext_y);
        z1 = de_periodize(z0, z1, full_box_ext_z, half_box_ext_z);

        const SIMD_TYPE_F x = simd::lerp(x0, x1, t);
        const SIMD_TYPE_F y = simd::lerp(y0, y1, t);
        const SIMD_TYPE_F z = simd::lerp(z0, z1, t);

        simd::store(out_x + i, x);
        simd::store(out_y + i, y);
        simd::store(out_z + i, z);
    }
}

// clang-format off
void linear_interpolation_pbc_128(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
								  const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
								  const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
								  int64 count, float t, const mat3& sim_box)
// clang-format on
{
    const __m128 full_box_ext_x = simd::set_f128(sim_box[0][0]);
    const __m128 full_box_ext_y = simd::set_f128(sim_box[1][1]);
    const __m128 full_box_ext_z = simd::set_f128(sim_box[2][2]);

    const __m128 half_box_ext_x = simd::mul(full_box_ext_x, simd::set_f128(0.5f));
    const __m128 half_box_ext_y = simd::mul(full_box_ext_y, simd::set_f128(0.5f));
    const __m128 half_box_ext_z = simd::mul(full_box_ext_z, simd::set_f128(0.5f));

    for (int64 i = 0; i < count; i += 4) {
        __m128 x0 = simd::load_f128(in_x0 + i);
        __m128 y0 = simd::load_f128(in_y0 + i);
        __m128 z0 = simd::load_f128(in_z0 + i);

        __m128 x1 = simd::load_f128(in_x1 + i);
        __m128 y1 = simd::load_f128(in_y1 + i);
        __m128 z1 = simd::load_f128(in_z1 + i);

        x1 = de_periodize(x0, x1, full_box_ext_x, half_box_ext_x);
        y1 = de_periodize(y0, y1, full_box_ext_y, half_box_ext_y);
        z1 = de_periodize(z0, z1, full_box_ext_z, half_box_ext_z);

        const __m128 x = simd::lerp(x0, x1, t);
        const __m128 y = simd::lerp(y0, y1, t);
        const __m128 z = simd::lerp(z0, z1, t);

        simd::store(out_x + i, x);
        simd::store(out_y + i, y);
        simd::store(out_z + i, z);
    }
}

#ifdef __AVX__
// clang-format off
void linear_interpolation_pbc_256(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
								  const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
								  const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
								  int64 count, float t, const mat3& sim_box)
// clang-format on
{
    const __m256 full_box_ext_x = simd::set_f256(sim_box[0][0]);
    const __m256 full_box_ext_y = simd::set_f256(sim_box[1][1]);
    const __m256 full_box_ext_z = simd::set_f256(sim_box[2][2]);

    const __m256 half_box_ext_x = simd::mul(full_box_ext_x, simd::set_f256(0.5f));
    const __m256 half_box_ext_y = simd::mul(full_box_ext_y, simd::set_f256(0.5f));
    const __m256 half_box_ext_z = simd::mul(full_box_ext_z, simd::set_f256(0.5f));

    for (int64 i = 0; i < count; i += 8) {
        __m256 x0 = simd::load_f256(in_x0 + i);
        __m256 y0 = simd::load_f256(in_y0 + i);
        __m256 z0 = simd::load_f256(in_z0 + i);

        __m256 x1 = simd::load_f256(in_x1 + i);
        __m256 y1 = simd::load_f256(in_y1 + i);
        __m256 z1 = simd::load_f256(in_z1 + i);

        x1 = de_periodize(x0, x1, full_box_ext_x, half_box_ext_x);
        y1 = de_periodize(y0, y1, full_box_ext_y, half_box_ext_y);
        z1 = de_periodize(z0, z1, full_box_ext_z, half_box_ext_z);

        const __m256 x = simd::lerp(x0, x1, t);
        const __m256 y = simd::lerp(y0, y1, t);
        const __m256 z = simd::lerp(z0, z1, t);

        simd::store(out_x + i, x);
        simd::store(out_y + i, y);
        simd::store(out_z + i, z);
    }
}
#endif

// clang-format off
void cubic_interpolation(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
						 const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
						 const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
                         const float* RESTRICT in_x2, const float* RESTRICT in_y2, const float* RESTRICT in_z2,
						 const float* RESTRICT in_x3, const float* RESTRICT in_y3, const float* RESTRICT in_z3,
						 int64 count, float t)
// clang-format on
{
    for (int i = 0; i < count; i += SIMD_WIDTH) {
        const SIMD_TYPE_F x0 = SIMD_LOAD_F(in_x0 + i);
        const SIMD_TYPE_F y0 = SIMD_LOAD_F(in_y0 + i);
        const SIMD_TYPE_F z0 = SIMD_LOAD_F(in_z0 + i);

        const SIMD_TYPE_F x1 = SIMD_LOAD_F(in_x1 + i);
        const SIMD_TYPE_F y1 = SIMD_LOAD_F(in_y1 + i);
        const SIMD_TYPE_F z1 = SIMD_LOAD_F(in_z1 + i);

        const SIMD_TYPE_F x2 = SIMD_LOAD_F(in_x2 + i);
        const SIMD_TYPE_F y2 = SIMD_LOAD_F(in_y2 + i);
        const SIMD_TYPE_F z2 = SIMD_LOAD_F(in_z2 + i);

        const SIMD_TYPE_F x3 = SIMD_LOAD_F(in_x3 + i);
        const SIMD_TYPE_F y3 = SIMD_LOAD_F(in_y3 + i);
        const SIMD_TYPE_F z3 = SIMD_LOAD_F(in_z3 + i);

        const SIMD_TYPE_F x = simd::cubic_spline(x0, x1, x2, x3, t);
        const SIMD_TYPE_F y = simd::cubic_spline(y0, y1, y2, y3, t);
        const SIMD_TYPE_F z = simd::cubic_spline(z0, z1, z2, z3, t);

        simd::store(out_x + i, x);
        simd::store(out_y + i, y);
        simd::store(out_z + i, z);
    }
}

// clang-format off
void cubic_interpolation_128(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
							 const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
							 const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
							 const float* RESTRICT in_x2, const float* RESTRICT in_y2, const float* RESTRICT in_z2,
							 const float* RESTRICT in_x3, const float* RESTRICT in_y3, const float* RESTRICT in_z3,
							 int64 count, float t)
// clang-format on
{
    for (int i = 0; i < count; i += 4) {
        const __m128 x0 = simd::load_f128(in_x0 + i);
        const __m128 y0 = simd::load_f128(in_y0 + i);
        const __m128 z0 = simd::load_f128(in_z0 + i);

        const __m128 x1 = simd::load_f128(in_x1 + i);
        const __m128 y1 = simd::load_f128(in_y1 + i);
        const __m128 z1 = simd::load_f128(in_z1 + i);

        const __m128 x2 = simd::load_f128(in_x2 + i);
        const __m128 y2 = simd::load_f128(in_y2 + i);
        const __m128 z2 = simd::load_f128(in_z2 + i);

        const __m128 x3 = simd::load_f128(in_x3 + i);
        const __m128 y3 = simd::load_f128(in_y3 + i);
        const __m128 z3 = simd::load_f128(in_z3 + i);

        const __m128 x = simd::cubic_spline(x0, x1, x2, x3, t);
        const __m128 y = simd::cubic_spline(y0, y1, y2, y3, t);
        const __m128 z = simd::cubic_spline(z0, z1, z2, z3, t);

        simd::store(out_x + i, x);
        simd::store(out_y + i, y);
        simd::store(out_z + i, z);
    }
}

#ifdef __AVX__
// clang-format off
void cubic_interpolation_256(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
							 const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
							 const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
							 const float* RESTRICT in_x2, const float* RESTRICT in_y2, const float* RESTRICT in_z2,
							 const float* RESTRICT in_x3, const float* RESTRICT in_y3, const float* RESTRICT in_z3,
							 int64 count, float t)
// clang-format on
{
    for (int i = 0; i < count; i += 8) {
        const __m256 x0 = simd::load_f256(in_x0 + i);
        const __m256 y0 = simd::load_f256(in_y0 + i);
        const __m256 z0 = simd::load_f256(in_z0 + i);

        const __m256 x1 = simd::load_f256(in_x1 + i);
        const __m256 y1 = simd::load_f256(in_y1 + i);
        const __m256 z1 = simd::load_f256(in_z1 + i);

        const __m256 x2 = simd::load_f256(in_x2 + i);
        const __m256 y2 = simd::load_f256(in_y2 + i);
        const __m256 z2 = simd::load_f256(in_z2 + i);

        const __m256 x3 = simd::load_f256(in_x3 + i);
        const __m256 y3 = simd::load_f256(in_y3 + i);
        const __m256 z3 = simd::load_f256(in_z3 + i);

        const __m256 x = simd::cubic_spline(x0, x1, x2, x3, t);
        const __m256 y = simd::cubic_spline(y0, y1, y2, y3, t);
        const __m256 z = simd::cubic_spline(z0, z1, z2, z3, t);

        simd::store(out_x + i, x);
        simd::store(out_y + i, y);
        simd::store(out_z + i, z);
    }
}
#endif

// clang-format off
void cubic_interpolation_pbc_scalar(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
								    const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
								    const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
								    const float* RESTRICT in_x2, const float* RESTRICT in_y2, const float* RESTRICT in_z2,
								    const float* RESTRICT in_x3, const float* RESTRICT in_y3, const float* RESTRICT in_z3,
								    int64 count, float t, const mat3& sim_box)
// clang-format on
{
    const vec3 full_box_ext = sim_box * vec3(1);
    const vec3 half_box_ext = full_box_ext * 0.5f;

    for (int64 i = 0; i < count; i++) {
        float x0 = in_x0[i];
        float y0 = in_y0[i];
        float z0 = in_z0[i];

        float x1 = in_x1[i];
        float y1 = in_y1[i];
        float z1 = in_z1[i];

        float x2 = in_x2[i];
        float y2 = in_y2[i];
        float z2 = in_z2[i];

        float x3 = in_x3[i];
        float y3 = in_y3[i];
        float z3 = in_z3[i];

        x0 = de_periodize(x1, x0, full_box_ext.x, half_box_ext.x);
        x2 = de_periodize(x1, x2, full_box_ext.x, half_box_ext.x);
        x3 = de_periodize(x1, x3, full_box_ext.x, half_box_ext.x);

        y0 = de_periodize(y1, y0, full_box_ext.y, half_box_ext.y);
        y2 = de_periodize(y1, y2, full_box_ext.y, half_box_ext.y);
        y3 = de_periodize(y1, y3, full_box_ext.y, half_box_ext.y);

        z0 = de_periodize(z1, z0, full_box_ext.z, half_box_ext.z);
        z2 = de_periodize(z1, z2, full_box_ext.z, half_box_ext.z);
        z3 = de_periodize(z1, z3, full_box_ext.z, half_box_ext.z);

        const float x = math::cubic_spline(x0, x1, x2, x3, t);
        const float y = math::cubic_spline(y0, y1, y2, y3, t);
        const float z = math::cubic_spline(z0, z1, z2, z3, t);

        out_x[i] = x;
        out_y[i] = y;
        out_z[i] = z;
    }
}

// clang-format off
void cubic_interpolation_pbc(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
							 const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
                             const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
							 const float* RESTRICT in_x2, const float* RESTRICT in_y2, const float* RESTRICT in_z2,
							 const float* RESTRICT in_x3, const float* RESTRICT in_y3, const float* RESTRICT in_z3,
							 int64 count, float t, const mat3& sim_box)
// clang-format on
{
    const SIMD_TYPE_F full_box_ext_x = SIMD_SET_F(sim_box[0][0]);
    const SIMD_TYPE_F full_box_ext_y = SIMD_SET_F(sim_box[1][1]);
    const SIMD_TYPE_F full_box_ext_z = SIMD_SET_F(sim_box[2][2]);

    const SIMD_TYPE_F half_box_ext_x = simd::mul(full_box_ext_x, SIMD_SET_F(0.5f));
    const SIMD_TYPE_F half_box_ext_y = simd::mul(full_box_ext_y, SIMD_SET_F(0.5f));
    const SIMD_TYPE_F half_box_ext_z = simd::mul(full_box_ext_z, SIMD_SET_F(0.5f));

    for (int64 i = 0; i < count; i += SIMD_WIDTH) {
        const SIMD_TYPE_F x0 = SIMD_LOAD_F(in_x0 + i);
        const SIMD_TYPE_F y0 = SIMD_LOAD_F(in_y0 + i);
        const SIMD_TYPE_F z0 = SIMD_LOAD_F(in_z0 + i);

        const SIMD_TYPE_F x1 = SIMD_LOAD_F(in_x1 + i);
        const SIMD_TYPE_F y1 = SIMD_LOAD_F(in_y1 + i);
        const SIMD_TYPE_F z1 = SIMD_LOAD_F(in_z1 + i);

        const SIMD_TYPE_F x2 = SIMD_LOAD_F(in_x2 + i);
        const SIMD_TYPE_F y2 = SIMD_LOAD_F(in_y2 + i);
        const SIMD_TYPE_F z2 = SIMD_LOAD_F(in_z2 + i);

        const SIMD_TYPE_F x3 = SIMD_LOAD_F(in_x3 + i);
        const SIMD_TYPE_F y3 = SIMD_LOAD_F(in_y3 + i);
        const SIMD_TYPE_F z3 = SIMD_LOAD_F(in_z3 + i);

        const SIMD_TYPE_F dp_x0 = de_periodize(x1, x0, full_box_ext_x, half_box_ext_x);
        const SIMD_TYPE_F dp_x2 = de_periodize(x1, x2, full_box_ext_x, half_box_ext_x);
        const SIMD_TYPE_F dp_x3 = de_periodize(x1, x3, full_box_ext_x, half_box_ext_x);

        const SIMD_TYPE_F dp_y0 = de_periodize(y1, y0, full_box_ext_y, half_box_ext_y);
        const SIMD_TYPE_F dp_y2 = de_periodize(y1, y2, full_box_ext_y, half_box_ext_y);
        const SIMD_TYPE_F dp_y3 = de_periodize(y1, y3, full_box_ext_y, half_box_ext_y);

        const SIMD_TYPE_F dp_z0 = de_periodize(z1, z0, full_box_ext_z, half_box_ext_z);
        const SIMD_TYPE_F dp_z2 = de_periodize(z1, z2, full_box_ext_z, half_box_ext_z);
        const SIMD_TYPE_F dp_z3 = de_periodize(z1, z3, full_box_ext_z, half_box_ext_z);

        const SIMD_TYPE_F x = simd::cubic_spline(dp_x0, x1, dp_x2, dp_x3, t);
        const SIMD_TYPE_F y = simd::cubic_spline(dp_y0, y1, dp_y2, dp_y3, t);
        const SIMD_TYPE_F z = simd::cubic_spline(dp_z0, z1, dp_z2, dp_z3, t);

        simd::store(out_x + i, x);
        simd::store(out_y + i, y);
        simd::store(out_z + i, z);
    }
}

// clang-format off
void cubic_interpolation_pbc_128(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
								 const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
								 const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
								 const float* RESTRICT in_x2, const float* RESTRICT in_y2, const float* RESTRICT in_z2,
								 const float* RESTRICT in_x3, const float* RESTRICT in_y3, const float* RESTRICT in_z3,
								 int64 count, float t, const mat3& sim_box)
// clang-format on
{
    const __m128 full_box_ext_x = simd::set_f128(sim_box[0][0]);
    const __m128 full_box_ext_y = simd::set_f128(sim_box[1][1]);
    const __m128 full_box_ext_z = simd::set_f128(sim_box[2][2]);

    const __m128 half_box_ext_x = simd::mul(full_box_ext_x, simd::set_f128(0.5f));
    const __m128 half_box_ext_y = simd::mul(full_box_ext_y, simd::set_f128(0.5f));
    const __m128 half_box_ext_z = simd::mul(full_box_ext_z, simd::set_f128(0.5f));

    for (int64 i = 0; i < count; i += 4) {
        const __m128 x0 = simd::load_f128(in_x0 + i);
        const __m128 y0 = simd::load_f128(in_y0 + i);
        const __m128 z0 = simd::load_f128(in_z0 + i);

        const __m128 x1 = simd::load_f128(in_x1 + i);
        const __m128 y1 = simd::load_f128(in_y1 + i);
        const __m128 z1 = simd::load_f128(in_z1 + i);

        const __m128 x2 = simd::load_f128(in_x2 + i);
        const __m128 y2 = simd::load_f128(in_y2 + i);
        const __m128 z2 = simd::load_f128(in_z2 + i);

        const __m128 x3 = simd::load_f128(in_x3 + i);
        const __m128 y3 = simd::load_f128(in_y3 + i);
        const __m128 z3 = simd::load_f128(in_z3 + i);

        const __m128 dp_x0 = de_periodize(x1, x0, full_box_ext_x, half_box_ext_x);
        const __m128 dp_x2 = de_periodize(x1, x2, full_box_ext_x, half_box_ext_x);
        const __m128 dp_x3 = de_periodize(x1, x3, full_box_ext_x, half_box_ext_x);

        const __m128 dp_y0 = de_periodize(y1, y0, full_box_ext_y, half_box_ext_y);
        const __m128 dp_y2 = de_periodize(y1, y2, full_box_ext_y, half_box_ext_y);
        const __m128 dp_y3 = de_periodize(y1, y3, full_box_ext_y, half_box_ext_y);

        const __m128 dp_z0 = de_periodize(z1, z0, full_box_ext_z, half_box_ext_z);
        const __m128 dp_z2 = de_periodize(z1, z2, full_box_ext_z, half_box_ext_z);
        const __m128 dp_z3 = de_periodize(z1, z3, full_box_ext_z, half_box_ext_z);

        const __m128 x = simd::cubic_spline(dp_x0, x1, dp_x2, dp_x3, t);
        const __m128 y = simd::cubic_spline(dp_y0, y1, dp_y2, dp_y3, t);
        const __m128 z = simd::cubic_spline(dp_z0, z1, dp_z2, dp_z3, t);

        simd::store(out_x + i, x);
        simd::store(out_y + i, y);
        simd::store(out_z + i, z);
    }
}

#ifdef __AVX__
// clang-format off
void cubic_interpolation_pbc_256(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
								 const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
								 const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
								 const float* RESTRICT in_x2, const float* RESTRICT in_y2, const float* RESTRICT in_z2,
								 const float* RESTRICT in_x3, const float* RESTRICT in_y3, const float* RESTRICT in_z3,
								 int64 count, float t, const mat3& sim_box)
// clang-format on
{
    const __m256 full_box_ext_x = simd::set_f256(sim_box[0][0]);
    const __m256 full_box_ext_y = simd::set_f256(sim_box[1][1]);
    const __m256 full_box_ext_z = simd::set_f256(sim_box[2][2]);

    const __m256 half_box_ext_x = simd::mul(full_box_ext_x, simd::set_f256(0.5f));
    const __m256 half_box_ext_y = simd::mul(full_box_ext_y, simd::set_f256(0.5f));
    const __m256 half_box_ext_z = simd::mul(full_box_ext_z, simd::set_f256(0.5f));

    for (int64 i = 0; i < count; i += 8) {
        __m256 x0 = simd::load_f256(in_x0 + i);
        __m256 y0 = simd::load_f256(in_y0 + i);
        __m256 z0 = simd::load_f256(in_z0 + i);

        __m256 x1 = simd::load_f256(in_x1 + i);
        __m256 y1 = simd::load_f256(in_y1 + i);
        __m256 z1 = simd::load_f256(in_z1 + i);

        __m256 x2 = simd::load_f256(in_x2 + i);
        __m256 y2 = simd::load_f256(in_y2 + i);
        __m256 z2 = simd::load_f256(in_z2 + i);

        __m256 x3 = simd::load_f256(in_x3 + i);
        __m256 y3 = simd::load_f256(in_y3 + i);
        __m256 z3 = simd::load_f256(in_z3 + i);

        const __m256 dp_x0 = de_periodize(x1, x0, full_box_ext_x, half_box_ext_x);
        const __m256 dp_x2 = de_periodize(x1, x2, full_box_ext_x, half_box_ext_x);
        const __m256 dp_x3 = de_periodize(x1, x3, full_box_ext_x, half_box_ext_x);

        const __m256 dp_y0 = de_periodize(y1, y0, full_box_ext_y, half_box_ext_y);
        const __m256 dp_y2 = de_periodize(y1, y2, full_box_ext_y, half_box_ext_y);
        const __m256 dp_y3 = de_periodize(y1, y3, full_box_ext_y, half_box_ext_y);

        const __m256 dp_z0 = de_periodize(z1, z0, full_box_ext_z, half_box_ext_z);
        const __m256 dp_z2 = de_periodize(z1, z2, full_box_ext_z, half_box_ext_z);
        const __m256 dp_z3 = de_periodize(z1, z3, full_box_ext_z, half_box_ext_z);

        const __m256 x = simd::cubic_spline(dp_x0, x1, dp_x2, dp_x3, t);
        const __m256 y = simd::cubic_spline(dp_y0, y1, dp_y2, dp_y3, t);
        const __m256 z = simd::cubic_spline(dp_z0, z1, dp_z2, dp_z3, t);

        simd::store(out_x + i, x);
        simd::store(out_y + i, y);
        simd::store(out_z + i, z);
    }
}

#endif

// clang-format off
void compute_velocities_scalar(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
							   const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
							   const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
							   int64 count, float dt)
// clang-format on
{
    for (int64 i = 0; i < count; i++) {
        out_x[i] = (in_x1[i] - in_x0[i]) * dt;
        out_y[i] = (in_y1[i] - in_y0[i]) * dt;
        out_z[i] = (in_z1[i] - in_z0[i]) * dt;
    }
}

// clang-format off
void compute_velocities(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
						const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
						const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
						int64 count, float dt)
// clang-format on
{
    const SIMD_TYPE_F dt_wide = SIMD_SET_F(dt);

    for (int64 i = 0; i < count; i += SIMD_WIDTH) {
        const SIMD_TYPE_F x0 = SIMD_LOAD_F(in_x0 + i);
        const SIMD_TYPE_F y0 = SIMD_LOAD_F(in_y0 + i);
        const SIMD_TYPE_F z0 = SIMD_LOAD_F(in_z0 + i);

        const SIMD_TYPE_F x1 = SIMD_LOAD_F(in_x1 + i);
        const SIMD_TYPE_F y1 = SIMD_LOAD_F(in_y1 + i);
        const SIMD_TYPE_F z1 = SIMD_LOAD_F(in_z1 + i);

        const SIMD_TYPE_F dx = simd::mul(simd::sub(x1, x0), dt_wide);
        const SIMD_TYPE_F dy = simd::mul(simd::sub(y1, y0), dt_wide);
        const SIMD_TYPE_F dz = simd::mul(simd::sub(z1, z0), dt_wide);

        SIMD_STORE(out_x + i, dx);
        SIMD_STORE(out_y + i, dy);
        SIMD_STORE(out_z + i, dz);
    }
}

// clang-format off
void compute_velocities_128(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
							const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
							const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
							int64 count, float dt)
// clang-format on
{
    const __m128 dt128 = simd::set_f128(dt);

    for (int i = 0; i < count; i += 4) {
        const __m128 x0 = simd::load_f128(in_x0 + i);
        const __m128 y0 = simd::load_f128(in_y0 + i);
        const __m128 z0 = simd::load_f128(in_z0 + i);

        const __m128 x1 = simd::load_f128(in_x1 + i);
        const __m128 y1 = simd::load_f128(in_y1 + i);
        const __m128 z1 = simd::load_f128(in_z1 + i);

        const __m128 dx = simd::mul(simd::sub(x1, x0), dt128);
        const __m128 dy = simd::mul(simd::sub(y1, y0), dt128);
        const __m128 dz = simd::mul(simd::sub(z1, z0), dt128);

        simd::store(out_x + i, dx);
        simd::store(out_y + i, dy);
        simd::store(out_z + i, dz);
    }
}

#ifdef __AVX__
// clang-format off
void compute_velocities_256(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
							const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
							const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
							int64 count, float dt)
// clang-format on
{
    const __m256 dt256 = simd::set_f256(dt);

    for (int i = 0; i < count; i += 8) {
        const __m256 x0 = simd::load_f256(in_x0 + i);
        const __m256 y0 = simd::load_f256(in_y0 + i);
        const __m256 z0 = simd::load_f256(in_z0 + i);

        const __m256 x1 = simd::load_f256(in_x1 + i);
        const __m256 y1 = simd::load_f256(in_y1 + i);
        const __m256 z1 = simd::load_f256(in_z1 + i);

        const __m256 dx = simd::mul(simd::sub(x1, x0), dt256);
        const __m256 dy = simd::mul(simd::sub(y1, y0), dt256);
        const __m256 dz = simd::mul(simd::sub(z1, z0), dt256);

        simd::store(out_x + i, dx);
        simd::store(out_y + i, dy);
        simd::store(out_z + i, dz);
    }
}
#endif

// clang-format off
void compute_velocities_pbc(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
							const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
                            const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
							int64 count, float dt, const mat3& sim_box)
// clang-format on
{
    const SIMD_TYPE_F full_box_ext_x = SIMD_SET_F(sim_box[0][0]);
    const SIMD_TYPE_F full_box_ext_y = SIMD_SET_F(sim_box[1][1]);
    const SIMD_TYPE_F full_box_ext_z = SIMD_SET_F(sim_box[2][2]);

    const SIMD_TYPE_F half_box_ext_x = simd::mul(full_box_ext_x, SIMD_SET_F(0.5f));
    const SIMD_TYPE_F half_box_ext_y = simd::mul(full_box_ext_y, SIMD_SET_F(0.5f));
    const SIMD_TYPE_F half_box_ext_z = simd::mul(full_box_ext_z, SIMD_SET_F(0.5f));

    const SIMD_TYPE_F dt_wide = SIMD_SET_F(dt);

    for (int i = 0; i < count; i += SIMD_WIDTH) {
        const SIMD_TYPE_F x0 = SIMD_LOAD_F(in_x0 + i);
        const SIMD_TYPE_F y0 = SIMD_LOAD_F(in_y0 + i);
        const SIMD_TYPE_F z0 = SIMD_LOAD_F(in_z0 + i);

        const SIMD_TYPE_F x1 = SIMD_LOAD_F(in_x1 + i);
        const SIMD_TYPE_F y1 = SIMD_LOAD_F(in_y1 + i);
        const SIMD_TYPE_F z1 = SIMD_LOAD_F(in_z1 + i);

        const SIMD_TYPE_F dp_x0 = de_periodize(x1, x0, full_box_ext_x, half_box_ext_x);
        const SIMD_TYPE_F dp_y0 = de_periodize(y1, y0, full_box_ext_y, half_box_ext_y);
        const SIMD_TYPE_F dp_z0 = de_periodize(z1, z0, full_box_ext_z, half_box_ext_z);

        const SIMD_TYPE_F dx = simd::mul(simd::sub(x1, dp_x0), dt_wide);
        const SIMD_TYPE_F dy = simd::mul(simd::sub(y1, dp_y0), dt_wide);
        const SIMD_TYPE_F dz = simd::mul(simd::sub(z1, dp_z0), dt_wide);

        simd::store(out_x + i, dx);
        simd::store(out_y + i, dy);
        simd::store(out_z + i, dz);
    }
}

// clang-format off
void compute_velocities_pbc_128(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
								const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
                                const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
								int64 count, float dt, const mat3& sim_box)
// clang-format on
{
    const __m128 full_box_ext_x = simd::set_f128(sim_box[0][0]);
    const __m128 full_box_ext_y = simd::set_f128(sim_box[1][1]);
    const __m128 full_box_ext_z = simd::set_f128(sim_box[2][2]);

    const __m128 half_box_ext_x = simd::mul(full_box_ext_x, simd::set_f128(0.5f));
    const __m128 half_box_ext_y = simd::mul(full_box_ext_y, simd::set_f128(0.5f));
    const __m128 half_box_ext_z = simd::mul(full_box_ext_z, simd::set_f128(0.5f));

    const __m128 dt128 = simd::set_f128(dt);

    for (int i = 0; i < count; i += 4) {
        const __m128 x0 = simd::load_f128(in_x0 + i);
        const __m128 y0 = simd::load_f128(in_y0 + i);
        const __m128 z0 = simd::load_f128(in_z0 + i);

        const __m128 x1 = simd::load_f128(in_x1 + i);
        const __m128 y1 = simd::load_f128(in_y1 + i);
        const __m128 z1 = simd::load_f128(in_z1 + i);

        const __m128 dp_x0 = de_periodize(x1, x0, full_box_ext_x, half_box_ext_x);
        const __m128 dp_y0 = de_periodize(y1, y0, full_box_ext_y, half_box_ext_y);
        const __m128 dp_z0 = de_periodize(z1, z0, full_box_ext_z, half_box_ext_z);

        const __m128 dx = simd::mul(simd::sub(x1, dp_x0), dt128);
        const __m128 dy = simd::mul(simd::sub(y1, dp_y0), dt128);
        const __m128 dz = simd::mul(simd::sub(z1, dp_z0), dt128);

        simd::store(out_x + i, dx);
        simd::store(out_y + i, dy);
        simd::store(out_z + i, dz);
    }
}

#ifdef __AVX__
// clang-format off
void compute_velocities_pbc_256(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
								const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
                                const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
								int64 count, float dt, const mat3& sim_box)
// clang-format on
{
    const __m256 full_box_ext_x = simd::set_f256(sim_box[0][0]);
    const __m256 full_box_ext_y = simd::set_f256(sim_box[1][1]);
    const __m256 full_box_ext_z = simd::set_f256(sim_box[2][2]);

    const __m256 half_box_ext_x = simd::mul(full_box_ext_x, simd::set_f256(0.5f));
    const __m256 half_box_ext_y = simd::mul(full_box_ext_y, simd::set_f256(0.5f));
    const __m256 half_box_ext_z = simd::mul(full_box_ext_z, simd::set_f256(0.5f));

    const __m256 dt256 = simd::set_f256(dt);

    for (int i = 0; i < count; i += 8) {
        const __m256 x0 = simd::load_f256(in_x0 + i);
        const __m256 y0 = simd::load_f256(in_y0 + i);
        const __m256 z0 = simd::load_f256(in_z0 + i);
        const __m256 x1 = simd::load_f256(in_x1 + i);
        const __m256 y1 = simd::load_f256(in_y1 + i);
        const __m256 z1 = simd::load_f256(in_z1 + i);

        const __m256 dp_x0 = de_periodize(x1, x0, full_box_ext_x, half_box_ext_x);
        const __m256 dp_y0 = de_periodize(y1, y0, full_box_ext_y, half_box_ext_y);
        const __m256 dp_z0 = de_periodize(z1, z0, full_box_ext_z, half_box_ext_z);

        const __m256 dx = simd::mul(simd::sub(x1, dp_x0), dt256);
        const __m256 dy = simd::mul(simd::sub(y1, dp_y0), dt256);
        const __m256 dz = simd::mul(simd::sub(z1, dp_z0), dt256);

        simd::store(out_x + i, dx);
        simd::store(out_y + i, dy);
        simd::store(out_z + i, dz);
    }
}
#endif

/*
void apply_pbc_residues(Array<vec3> positions, Array<const Residue> residues, const mat3& sim_box) {
    const glm_vec4 full_box_ext = _mm_set_ps(0.f, sim_box[2][2], sim_box[1][1], sim_box[0][0]);
    for (const auto& res : residues) {
        const auto count = res.atom_idx.end - res.atom_idx.beg;
        const glm_vec4 scl = _mm_set_ps1(1.f / (float)count);
        glm_vec4 center = _mm_setzero_ps();
        for (const auto& pos : positions.subarray(res.atom_idx.beg, count)) {
            const glm_vec4 pos_vec = *reinterpret_cast<const glm_vec4*>(&pos);
            center = glm_vec4_add(center, pos_vec);
        }
        center = glm_vec4_mul(center, scl);
        const glm_vec4 pbc_cent = apply_pbc(center, full_box_ext);
        const glm_vec4 delta = glm_vec4_sub(pbc_cent, center);

        // @TODO: if delta is zero, skip this
        // if (!all_zero(delta)) {
        for (auto& pos : positions.subarray(res.atom_idx.beg, count)) {
            glm_vec4& pos_vec = *reinterpret_cast<glm_vec4*>(&pos);
            pos_vec = glm_vec4_add(pos_vec, delta);
        }
        //}
    }
}

void apply_pbc_chains(Array<vec3> positions, Array<const Chain> chains, Array<const Residue> residues, const mat3& sim_box) {
    const glm_vec4 full_box_ext = _mm_set_ps(0.f, sim_box[2][2], sim_box[1][1], sim_box[0][0]);
    for (const auto& chain : chains) {
        const auto beg_idx = residues[chain.res_idx.beg].atom_idx.beg;
        const auto end_idx = residues[chain.res_idx.end - 1].atom_idx.end;
        const auto count = end_idx - beg_idx;
        const glm_vec4 scl = _mm_set_ps1(1.f / (float)count);
        glm_vec4 center = _mm_setzero_ps();
        for (const auto& pos : positions.subarray(beg_idx, count)) {
            const glm_vec4 pos_vec = *reinterpret_cast<const glm_vec4*>(&pos);
            center = glm_vec4_add(center, pos_vec);
        }
        center = glm_vec4_mul(center, scl);
        const glm_vec4 pbc_cent = apply_pbc(center, full_box_ext);
        const glm_vec4 delta = glm_vec4_sub(pbc_cent, center);

        // if (!all_zero(delta)) {
        for (auto& pos : positions.subarray(beg_idx, count)) {
            const glm_vec4 pos_vec = *reinterpret_cast<glm_vec4*>(&pos);
            const glm_vec4 res = glm_vec4_add(pos_vec, delta);
            pos = *reinterpret_cast<const vec3*>(&res);
        }
        //}
    }
}
*/

void apply_pbc(float* RESTRICT x, float* RESTRICT y, float* RESTRICT z, const float* RESTRICT mass, ArrayView<const Sequence> sequences, const mat3& sim_box) {
    const vec3 ext = sim_box * vec3(1.0f);
    const vec3 one_over_ext = 1.0f / ext;

    for (const auto& seq : sequences) {
        const auto& range = seq.atom_range;
        float* seq_x = x + range.beg;
        float* seq_y = y + range.beg;
        float* seq_z = z + range.beg;
        const float* seq_mass = mass + range.beg;
        const vec3 com = compute_com_periodic(seq_x, seq_y, seq_z, seq_mass, range.size(), sim_box);
        const vec3 com_dp = math::fract(com * one_over_ext) * ext;

        for (int i = 0; i < range.size(); i++) {
            seq_x[i] = de_periodize(com_dp.x, seq_x[i], ext.x, ext.x * 0.5f);
            seq_y[i] = de_periodize(com_dp.y, seq_y[i], ext.y, ext.y * 0.5f);
            seq_z[i] = de_periodize(com_dp.z, seq_z[i], ext.z, ext.z * 0.5f);
        }

        /*
        const vec3 delta = com_dp - com;
        if (math::dot(delta, delta) > 0.0001f) {
            translate_ref(seq_x, seq_y, seq_z, range.size(), delta);
        }
        */
    }
}

inline bool covelent_bond_heuristic(float x0, float y0, float z0, Element e0, float x1, float y1, float z1, Element e1) {
    const float d = element::covalent_radius(e0) + element::covalent_radius(e1);
    const float d_max = d + 0.3f;
    const float d_min = d - 0.5f;
    const float dx = x1 - x0;
    const float dy = y1 - y0;
    const float dz = z1 - z0;
    const float d2 = dx * dx + dy * dy + dz * dz;
    return (d_min * d_min) < d2 && d2 < (d_max * d_max);
}

bool has_covalent_bond(const Residue& res_a, const Residue& res_b) { return (res_a.bond_idx.beg < res_b.bond_idx.end && res_b.bond_idx.beg < res_a.bond_idx.end); }

bool valid_segment(const BackboneSegment& segment) { return segment.ca_idx != -1 && segment.c_idx != -1 && segment.n_idx != -1 && segment.o_idx != -1; }

// Computes covalent bonds between a set of atoms with given positions and elements.
// The approach is inspired by the technique used in NGL (https://github.com/arose/ngl)
DynamicArray<Bond> compute_covalent_bonds(ArrayView<Residue> residues, const float* pos_x, const float* pos_y, const float* pos_z, const Element* element, int64 count) {
    UNUSED(count);

    if (residues.count == 0) {
        LOG_WARNING("Cannot compute covalent bonds, no residues were given.");
        return {};
    }

    constexpr float max_covelent_bond_length = 4.0f;
    DynamicArray<Bond> bonds;
    spatialhash::Frame frame;

    // @NOTE: The assumtion is that a bond is either within a single residue or between concecutive residues.
    for (ResIdx ri = 0; ri < (ResIdx)residues.size(); ri++) {
        auto& res = residues[ri];
        const float* res_pos_x = pos_x + res.atom_range.beg;
        const float* res_pos_y = pos_y + res.atom_range.beg;
        const float* res_pos_z = pos_z + res.atom_range.beg;

        spatialhash::compute_frame(&frame, res_pos_x, res_pos_y, res_pos_z, res.atom_range.size(), vec3(max_covelent_bond_length));

        if (ri > 0) {
            // Include potential shared bonds from previous residue
            res.bond_idx.beg = residues[ri - 1].bond_idx.end_internal;
            res.bond_idx.beg_internal = res.bond_idx.end_internal = res.bond_idx.end = residues[ri - 1].bond_idx.end;
        } else {
            res.bond_idx.beg = res.bond_idx.end = (BondIdx)bonds.size();
            res.bond_idx.beg_internal = res.bond_idx.end_internal = res.bond_idx.end = (BondIdx)bonds.size();
        }

        // Internal bonds
        for (AtomIdx i = res.atom_range.beg; i < res.atom_range.end; i++) {
            const vec3& pos_xyz = {pos_x[i], pos_y[i], pos_z[i]};
            spatialhash::for_each_within(frame, pos_xyz, max_covelent_bond_length, [&bonds, &res, offset = res.atom_range.beg, pos_x, pos_y, pos_z, element, i](int j, const vec3& atom_j_pos) {
                (void)atom_j_pos;
                j += offset;  // @NOTE: Map residue idx j (given by spatial hash) back to full atomic idx
                const bool has_bond = covelent_bond_heuristic(pos_x[i], pos_y[i], pos_z[i], element[i], pos_x[j], pos_y[j], pos_z[j], element[j]);
                if (i < j && has_bond) {
                    bonds.push_back({{i, j}});
                    res.bond_idx.end++;
                }
            });
        }
        res.bond_idx.end_internal = res.bond_idx.end;

        // Locate potential external bonds to next residue
        if (ri < (ResIdx)residues.size() - 1) {
            auto& next_res = residues[ri + 1];
            const float* next_res_pos_x = pos_x + next_res.atom_range.beg;
            const float* next_res_pos_y = pos_y + next_res.atom_range.beg;
            const float* next_res_pos_z = pos_z + next_res.atom_range.beg;

            spatialhash::compute_frame(&frame, next_res_pos_x, next_res_pos_y, next_res_pos_z, next_res.atom_range.size(), vec3(max_covelent_bond_length));

            for (AtomIdx i = res.atom_range.beg; i < res.atom_range.end; i++) {
                const vec3& pos_xyz = {pos_x[i], pos_y[i], pos_z[i]};
                spatialhash::for_each_within(frame, pos_xyz, max_covelent_bond_length,
                                             [&bonds, &res, offset = next_res.atom_range.beg, pos_x, pos_y, pos_z, element, i](int j, const vec3& atom_j_pos) {
                                                 (void)atom_j_pos;
                                                 j += offset;  // @NOTE: Map residue idx j (given by spatial hash) back to full atomic idx
                                                 const bool has_bond = covelent_bond_heuristic(pos_x[i], pos_y[i], pos_z[i], element[i], pos_x[j], pos_y[j], pos_z[j], element[j]);
                                                 if (has_bond) {
                                                     bonds.push_back({{i, j}});
                                                     res.bond_idx.end++;
                                                 }
                                             });
            }
        }
    }

    return bonds;
}

DynamicArray<Bond> compute_covalent_bonds(const float* pos_x, const float* pos_y, const float* pos_z, const Element* element, int64 count) {
    constexpr float max_covelent_bond_length = 4.0f;
    spatialhash::Frame frame = spatialhash::compute_frame(pos_x, pos_y, pos_z, count, vec3(max_covelent_bond_length));
    DynamicArray<Bond> bonds;

    for (int i = 0; i < count; i++) {
        const vec3& pos_xyz = {pos_x[i], pos_y[i], pos_z[i]};
        spatialhash::for_each_within(frame, pos_xyz, max_covelent_bond_length, [&bonds, pos_x, pos_y, pos_z, element, i](int j, const vec3& atom_j_pos) {
            (void)atom_j_pos;
            bool has_bond = covelent_bond_heuristic(pos_x[i], pos_y[i], pos_z[i], element[i], pos_x[j], pos_y[j], pos_z[j], element[j]);
            if (i < j && has_bond) {
                bonds.push_back({{i, j}});
            }
        });
    }

    return bonds;
}

// @NOTE this method is sub-optimal and can surely be improved...
// Residues should have no more than 2 potential connections to other residues.
DynamicArray<Sequence> compute_sequences(ArrayView<const Residue> residues) {

    /*
DynamicArray<Bond> residue_bonds;
for (ResIdx i = 0; i < (ResIdx)residues.size() - 1; i++) {
    if (has_covalent_bond(residues[i], residues[i + 1])) {
        residue_bonds.push_back({{i, i + 1}});
    }
}

if (residue_bonds.size() == 0) {
    // No residue bonds, return residues as individual sequences
    DynamicArray<Sequence> seq(residues.size());
    for (int64 i = 0; i < residues.size(); i++) {
        seq[i].atom_range = residues[i].atom_range;
        seq[i].res_range = {(ResIdx)i, (ResIdx)i + 1};
    }
    return seq;
}

DynamicArray<int> residue_sequences(residues.size(), -1);
if (residue_bonds.size() > 0) {
    int curr_seq_idx = 0;
    int res_bond_idx = 0;
    for (int i = 0; i < residues.count; i++) {
        if (residue_sequences[i] == -1) residue_sequences[i] = curr_seq_idx++;
        for (; res_bond_idx < residue_bonds.size(); res_bond_idx++) {
            const auto& res_bond = residue_bonds[res_bond_idx];
            if (i == res_bond.idx[0]) {
                residue_sequences[res_bond.idx[1]] = residue_sequences[res_bond.idx[0]];
            } else if (res_bond.idx[0] > i)
                break;
        }
    }
}
    */

    DynamicArray<Sequence> seq;
    seq.push_back({"", {0, 1}, residues[0].atom_range});
    for (ResIdx i = 0; i < (ResIdx)residues.size() - 1; i++) {
        if (has_covalent_bond(residues[i], residues[i + 1])) {
            seq.back().res_range.end++;
            seq.back().atom_range.end = residues[i + 1].atom_range.end;
        } else {
            seq.push_back({"", {i + 1, i + 2}, residues[i + 1].atom_range});
        }
    }

    return seq;

    /*
DynamicArray<Sequence> sequences;
int curr_seq_idx = -1;
for (int i = 0; i < residue_sequences.size(); i++) {
    if (residue_sequences[i] != curr_seq_idx) {
        sequences.push_back({});
        sequences.back().res_range = {(ResIdx)i, (ResIdx)i};
    }
    if (sequences.size() > 0) {
        sequences.back().res_range.end++;
    }
}

for (auto& s : sequences) {
    s.atom_range.beg = residues[s.res_range.beg].atom_range.beg;
    s.atom_range.end = residues[s.res_range.end - 1].atom_range.end;
}

return sequences;
    */
}

DynamicArray<BackboneSequence> compute_backbone_sequences(ArrayView<const BackboneSegment> segments, ArrayView<const Residue> residues) {
    if (segments.count == 0) return {};
    ASSERT(segments.count == residues.count);

    DynamicArray<BackboneSequence> bb_sequences;
    for (ResIdx i = 0; i < (ResIdx)residues.size(); i++) {
        if (valid_segment(segments[i])) {
            bb_sequences.push_back({i, i + 1});
            while (i < (ResIdx)residues.size() - 1 && valid_segment(segments[i + 1]) && has_covalent_bond(residues[i], residues[i + 1])) {
                bb_sequences.back().end++;
                i++;
            }
        }
    }

    return bb_sequences;
}

template <int64 N>
bool match(const Label& lbl, const char (&cstr)[N]) {
    for (int64 i = 0; i < N; i++) {
        if (tolower(lbl[i]) != tolower(cstr[i])) return false;
    }
    return true;
}

DynamicArray<BackboneSegment> compute_backbone_segments(ArrayView<const Residue> residues, ArrayView<const Label> atom_labels) {
    DynamicArray<BackboneSegment> segments;
    int64 invalid_segments = 0;
    constexpr int32 min_atom_count = 4;  // Must contain at least 4 atoms to be considered as an amino acid.
    for (auto& res : residues) {
        const int32 atom_count = res.atom_range.end - res.atom_range.beg;
        if (atom_count < min_atom_count) {
            segments.push_back({-1, -1, -1, -1});
            invalid_segments++;
            continue;
        }

        BackboneSegment seg{};

        // find atoms
        for (int32 i = res.atom_range.beg; i < res.atom_range.end; i++) {
            const auto& lbl = atom_labels[i];
            if (seg.ca_idx == -1 && match(lbl, "CA")) seg.ca_idx = i;
            if (seg.n_idx == -1 && match(lbl, "N")) seg.n_idx = i;
            if (seg.c_idx == -1 && match(lbl, "C")) seg.c_idx = i;
            if (seg.o_idx == -1 && match(lbl, "O")) seg.o_idx = i;
        }

        // Could not match "O"
        if (seg.o_idx == -1) {
            // Pick first atom containing O after C atom
            for (int32 i = seg.c_idx; i < res.atom_range.end; i++) {
                const auto& lbl = atom_labels[i];
                if (lbl[0] == 'o' || lbl[0] == 'O') {
                    seg.o_idx = i;
                    break;
                }
            }
        }

        if (!valid_segment(seg)) {
            // LOG_ERROR("Could not identify all backbone indices for residue %s.", res.name.beg());
            invalid_segments++;
        }
        //} else {
        //    invalid_segments++;
        //}
        segments.push_back(seg);
    }

    if (invalid_segments == segments.size()) return {};

    return segments;
}

/*
DynamicArray<SplineSegment> compute_spline(Array<const vec3> atom_pos, Array<const uint32> colors, Array<const BackboneSegment> backbone, int32 num_subdivisions, float tension) {
    if (backbone.count < 4) return {};

    DynamicArray<vec3> p_tmp;
    DynamicArray<vec3> o_tmp;
    DynamicArray<vec3> c_tmp;
    DynamicArray<int> ca_idx;

    // Pad front with duplicated vectors
    auto d_p0 = atom_pos[backbone[1].ca_idx] - atom_pos[backbone[0].ca_idx];
    p_tmp.push_back(atom_pos[backbone[0].ca_idx] - d_p0);

    auto d_o0 = atom_pos[backbone[1].o_idx] - atom_pos[backbone[0].o_idx];
    o_tmp.push_back(atom_pos[backbone[0].o_idx] - d_o0);

    auto d_c0 = atom_pos[backbone[1].c_idx] - atom_pos[backbone[0].c_idx];
    c_tmp.push_back(atom_pos[backbone[0].c_idx] - d_c0);

    ca_idx.push_back(backbone[0].ca_idx);

    // Fetch vectors
    const int size = (int)(backbone.count);
    for (auto i = 0; i < size; i++) {
        p_tmp.push_back(atom_pos[backbone[i].ca_idx]);
        o_tmp.push_back(atom_pos[backbone[i].o_idx]);
        c_tmp.push_back(atom_pos[backbone[i].c_idx]);
        ca_idx.push_back(backbone[i].ca_idx);
    }

    // Pad back with duplicated vectors
    auto d_pn = atom_pos[backbone[size - 1].ca_idx] - atom_pos[backbone[size - 2].ca_idx];
    p_tmp.push_back(atom_pos[backbone[size - 1].ca_idx] + d_pn);
    p_tmp.push_back(p_tmp.back() + d_pn);

    auto d_on = atom_pos[backbone[size - 1].o_idx] - atom_pos[backbone[size - 2].o_idx];
    o_tmp.push_back(atom_pos[backbone[size - 1].o_idx] + d_on);
    o_tmp.push_back(o_tmp.back() + d_on);

    auto d_cn = atom_pos[backbone[size - 1].c_idx] - atom_pos[backbone[size - 2].c_idx];
    c_tmp.push_back(atom_pos[backbone[size - 1].c_idx] + d_cn);
    c_tmp.push_back(c_tmp.back() + d_cn);

    ca_idx.push_back(backbone[size - 1].ca_idx);
    ca_idx.push_back(backbone[size - 1].ca_idx);

    // @NOTE: Flip direction of support vector (O <- C) if pointing the 'wrong way' as seen from previous segment
    for (int64 i = 1; i < o_tmp.size(); i++) {
        vec3 v0 = o_tmp[i - 1] - c_tmp[i - 1];
        vec3 v1 = o_tmp[i] - c_tmp[i];

        if (glm::dot(v0, v1) < 0) {
            o_tmp[i] = c_tmp[i] - v1;
        }
    }

    DynamicArray<SplineSegment> segments;

    for (int64 i = 1; i < p_tmp.size() - 2; i++) {
        auto p0 = p_tmp[i - 1];
        auto p1 = p_tmp[i];
        auto p2 = p_tmp[i + 1];
        auto p3 = p_tmp[i + 2];

        auto o0 = o_tmp[i - 1];
        auto o1 = o_tmp[i];
        auto o2 = o_tmp[i + 1];
        auto o3 = o_tmp[i + 2];

        auto c0 = c_tmp[i - 1];
        auto c1 = c_tmp[i];
        auto c2 = c_tmp[i + 1];
        auto c3 = c_tmp[i + 2];

        uint32 idx = ca_idx[i];
        uint32 color = colors[idx];

        auto count = (i < (p_tmp.size() - 3)) ? num_subdivisions : num_subdivisions + 1;
        for (int n = 0; n < count; n++) {
            auto t = n / (float)(num_subdivisions);

            vec3 p = math::spline(p0, p1, p2, p3, t, tension);
            vec3 o = math::spline(o0, o1, o2, o3, t, tension);
            vec3 c = math::spline(c0, c1, c2, c3, t, tension);

            vec3 v_dir = math::normalize(o - c);

            vec3 tangent = math::spline_tangent(p0, p1, p2, p3, t, tension);
            vec3 normal = math::normalize(math::cross(v_dir, tangent));
            vec3 binormal = math::normalize(math::cross(tangent, normal));

            segments.push_back({p, tangent, normal, binormal, idx, color});
        }
    }

    return segments;
}
*/

DynamicArray<BackboneAngle> compute_backbone_angles(ArrayView<const BackboneSegment> backbone, const float* pos_x, const float* pos_y, const float* pos_z) {
    if (backbone.count == 0) return {};
    DynamicArray<BackboneAngle> angles(backbone.count);
    compute_backbone_angles(angles, backbone, pos_x, pos_y, pos_z);
    return angles;
}

void compute_backbone_angles(ArrayView<BackboneAngle> dst, ArrayView<const BackboneSegment> backbone_segments, const float* pos_x, const float* pos_y, const float* pos_z) {
    ASSERT(dst.count >= backbone_segments.count);
    float phi = 0, psi = 0;

    if (backbone_segments.size() < 2) {
        return;
    }

    ASSERT(valid_segment(backbone_segments[0]));
    vec3 n = {pos_x[backbone_segments[0].n_idx], pos_y[backbone_segments[0].n_idx], pos_z[backbone_segments[0].n_idx]};
    vec3 ca = {pos_x[backbone_segments[0].ca_idx], pos_y[backbone_segments[0].ca_idx], pos_z[backbone_segments[0].ca_idx]};
    vec3 c = {pos_x[backbone_segments[0].c_idx], pos_y[backbone_segments[0].c_idx], pos_z[backbone_segments[0].c_idx]};

    vec3 c_prev = c;
    vec3 n_next = {pos_x[backbone_segments[1].n_idx], pos_y[backbone_segments[1].n_idx], pos_z[backbone_segments[1].n_idx]};
    phi = 0.0f;
    psi = math::dihedral_angle(n, ca, c, n_next);
    dst[0] = {phi, psi};

    for (int64 i = 1; i < backbone_segments.count - 1; i++) {
        ASSERT(valid_segment(backbone_segments[i]));

        c_prev = c;
        n = n_next;
        ca = {pos_x[backbone_segments[i].ca_idx], pos_y[backbone_segments[i].ca_idx], pos_z[backbone_segments[i].ca_idx]};
        c = {pos_x[backbone_segments[i].c_idx], pos_y[backbone_segments[i].c_idx], pos_z[backbone_segments[i].c_idx]};
        n_next = {pos_x[backbone_segments[i + 1].n_idx], pos_y[backbone_segments[i + 1].n_idx], pos_z[backbone_segments[i + 1].n_idx]};

        phi = math::dihedral_angle(c_prev, n, ca, c);
        psi = math::dihedral_angle(n, ca, c, n_next);
        dst[i] = {phi, psi};
    }

    auto N = backbone_segments.count - 1;
    ASSERT(valid_segment(backbone_segments[N]));

    c_prev = c;
    n = n_next;
    ca = {pos_x[backbone_segments[N].ca_idx], pos_y[backbone_segments[N].ca_idx], pos_z[backbone_segments[N].ca_idx]};
    c = {pos_x[backbone_segments[N].c_idx], pos_y[backbone_segments[N].c_idx], pos_z[backbone_segments[N].c_idx]};

    phi = math::dihedral_angle(c_prev, n, ca, c);
    psi = 0.0f;
    dst[N] = {phi, psi};
}

DynamicArray<BackboneAngle> compute_backbone_angles(ArrayView<const BackboneSegment> segments, ArrayView<const BackboneSequence> sequences, const float* pos_x, const float* pos_y,
                                                    const float* pos_z) {
    if (segments.size() == 0) return {};
    DynamicArray<BackboneAngle> angles(segments.count);
    compute_backbone_angles(angles, segments, sequences, pos_x, pos_y, pos_z);
    return angles;
}

void compute_backbone_angles(ArrayView<BackboneAngle> dst, ArrayView<const BackboneSegment> segments, ArrayView<const BackboneSequence> sequences, const float* pos_x, const float* pos_y,
                             const float* pos_z) {
    for (const auto& seq : sequences) {
        compute_backbone_angles(dst.subarray(seq.beg, seq.end - seq.beg), segments.subarray(seq.beg, seq.end - seq.beg), pos_x, pos_y, pos_z);
    }
}

void init_backbone_angles_trajectory(BackboneAnglesTrajectory* data, const MoleculeDynamic& dynamic) {
    ASSERT(data);
    if (!dynamic.molecule || !dynamic.trajectory) return;

    if (data->angle_data) {
        FREE(data->angle_data.ptr);
    }

    int32 alloc_count = (int32)dynamic.molecule.backbone.segments.count * (int32)dynamic.trajectory.frame_buffer.count;
    data->num_segments = (int32)dynamic.molecule.backbone.segments.count;
    data->num_frames = 0;
    data->angle_data = {(BackboneAngle*)CALLOC(alloc_count, sizeof(BackboneAngle)), alloc_count};
}

void free_backbone_angles_trajectory(BackboneAnglesTrajectory* data) {
    ASSERT(data);
    if (data->angle_data) {
        FREE(data->angle_data.ptr);
        *data = {};
    }
}

void compute_backbone_angles_trajectory(BackboneAnglesTrajectory* data, const MoleculeDynamic& dynamic) {
    ASSERT(dynamic);
    if (dynamic.trajectory.num_frames == 0 || dynamic.molecule.backbone.segments.count == 0) return;

    //@NOTE: Trajectory may be loading while this is taking place, therefore read num_frames once and stick to that
    const int32 traj_num_frames = dynamic.trajectory.num_frames;

    // @NOTE: If we are up to date, no need to compute anything
    if (traj_num_frames == data->num_frames) {
        return;
    }

    // @TODO: parallelize?
    // @NOTE: Only compute data for indices which are new
    for (int32 f_idx = data->num_frames; f_idx < traj_num_frames; f_idx++) {
        auto pos_x = get_trajectory_position_x(dynamic.trajectory, f_idx);
        auto pos_y = get_trajectory_position_y(dynamic.trajectory, f_idx);
        auto pos_z = get_trajectory_position_z(dynamic.trajectory, f_idx);

        ArrayView<BackboneAngle> frame_angles = get_backbone_angles(*data, f_idx);
        for (const auto& bb_seq : dynamic.molecule.backbone.sequences) {
            auto bb_segments = get_backbone(dynamic.molecule, bb_seq);
            auto bb_angles = frame_angles.subarray(bb_seq);

            if (bb_segments.size() < 2) {
                memset(bb_angles.ptr, 0, bb_angles.size_in_bytes());
            } else {
                compute_backbone_angles(bb_angles, bb_segments, pos_x.data(), pos_y.data(), pos_z.data());
            }
        }
    }
    data->num_frames = traj_num_frames;  // update current count
}

DynamicArray<float> compute_atom_radii(ArrayView<const Element> elements) {
    DynamicArray<float> radii(elements.size(), 0);
    compute_atom_radii(radii.data(), elements.data(), radii.size());
    return radii;
}

void compute_atom_radii(float* out_radius, const Element* element, int64 count) {
    for (int64 i = 0; i < count; i++) {
        out_radius[i] = element::vdw_radius(element[i]);
    }
}

DynamicArray<float> compute_atom_masses(ArrayView<const Element> elements) {
    DynamicArray<float> mass(elements.size(), 0);
    compute_atom_radii(mass.data(), elements.data(), mass.size());
    return mass;
}

void compute_atom_masses(float* out_mass, const Element* element, int64 count) {
    for (int64 i = 0; i < count; i++) {
        out_mass[i] = element::vdw_radius(element[i]);
    }
}

bool is_amino_acid(const Residue& res) { return aminoacid::get_from_string(res.name) != AminoAcid::Unknown; }

static constexpr CStringView dna_residues[12] = {"DA", "DA3", "DA5", "DC", "DC3", "DC5", "DG", "DG3", "DG5", "DT", "DT3", "DT5"};
bool is_dna(const Residue& res) {
    for (auto dna_res : dna_residues) {
        if (compare(res.name, dna_res)) return true;
    }
    return false;
}

DynamicArray<Label> get_unique_residue_types(const MoleculeStructure& mol) {
    DynamicArray<Label> types = {};
    Label cur_lbl = {};
    for (const auto& res : mol.residues) {
        if (res.name != cur_lbl) {
            types.push_back(res.name);
        }
    }
    return types;
}

DynamicArray<ResIdx> get_residues_by_name(const MoleculeStructure& mol, CStringView name) {
    DynamicArray<ResIdx> residues;
    for (ResIdx i = 0; i < (ResIdx)mol.residues.size(); i++) {
        if (mol.residues[i].name == name) {
            residues.push_back(i);
        }
    }
    return residues;
}

bool atom_ranges_match(const MoleculeStructure& mol, AtomRange range_a, AtomRange range_b) {
    if (range_a.size() == 0 || range_b.size() == 0) return false;
    if (range_a.size() != range_b.size()) return false;

    const auto ele_a = get_elements(mol).subarray(range_a);
    const auto ele_b = get_elements(mol).subarray(range_b);
    const auto lbl_a = get_labels(mol).subarray(range_a);
    const auto lbl_b = get_labels(mol).subarray(range_b);
    int a_i = range_a.beg;
    int b_i = range_b.beg;
    while (a_i != range_a.end && b_i != range_b.end) {
        if (ele_a[a_i] != ele_b[b_i]) return false;
        if (lbl_a[a_i] != lbl_b[b_i]) return false;
    }
    return true;
}

DynamicArray<AtomRange> find_equivalent_structures(const MoleculeStructure& mol, AtomRange ref) {
    DynamicArray<AtomRange> matches = {};

    const auto ele = get_elements(mol);
    const auto lbl = get_labels(mol);

    const auto ele_ref = ele.subarray(ref);
    const auto lbl_ref = lbl.subarray(ref);

    int j = 0;
    for (int i = 0; i < mol.atom.count; i++) {
        if (ref.beg <= i && i < ref.end) {
            j = 0;
            i = ref.end;
        }
        if (ele[i] == ele_ref[j] && (compare(lbl[i], lbl_ref[j]))) {
            j++;
        }
        if (j == ref.size()) {
            matches.push_back({i - j + 1, i + 1});
            j = 0;
        }
    }

    return matches;
}

// @NOTE: THIS IS STUPID
bool structure_match(const MoleculeStructure& mol, Bitfield mask, int mask_offset, int mask_count, int structure_offset) {

    const auto ele = get_elements(mol);
    //const auto lbl = get_labels(mol);
    // const auto res_idx = get_residue_indices(mol);

    for (int i = 0; i < mask_count; ++i) {
        if (mask[mask_offset + i]) {
            const auto ref_element = ele[mask_offset + i];
            const auto element = ele[structure_offset + i];
            if (element != ref_element) return false;

            //const auto& ref_label = lbl[mask_offset + i];
            //const auto& label = lbl[structure_offset + i];
            //if (compare(label, ref_label) == false) return false;

            // if (compare(mol.residues[res_idx[mask_offset + i]].name, mol.residues[res_idx[mask_offset + i]].name) == false) return false;
            // if (mol.residues[res_idx[mask_offset + i]].id != mol.residues[res_idx[mask_offset + i]].id) return false;
        }
    }

    return true;
};

DynamicArray<Bitfield> find_equivalent_structures(const MoleculeStructure& mol, Bitfield ref_mask) {
    DynamicArray<Bitfield> matches = {};

    const auto ele = get_elements(mol);
    const auto lbl = get_labels(mol);

    const int32 first_bit = bitfield::find_first_bit_set(ref_mask);
    const int32 last_bit = bitfield::find_last_bit_set(ref_mask);

    if (first_bit == -1 || last_bit == -1) {
        LOG_WARNING("find_equivalent_structures: Reference bitfield was empty");
        return matches;
    }

    int32 count = bitfield::number_of_bits_set(ref_mask);
    void* mem = TMP_MALLOC(count * sizeof(float) * 4);
    defer { TMP_FREE(mem); };
    float* x = (float*)mem + 0 * count;
    float* y = (float*)mem + 1 * count;
    float* z = (float*)mem + 2 * count;
    float* m = (float*)mem + 3 * count;

    bitfield::extract_data_from_mask(x, mol.atom.position.x, ref_mask);
    bitfield::extract_data_from_mask(y, mol.atom.position.y, ref_mask);
    bitfield::extract_data_from_mask(z, mol.atom.position.z, ref_mask);
    bitfield::extract_data_from_mask(m, mol.atom.mass, ref_mask);

    const EigenFrame ref_eigen = compute_eigen_frame(x, y, z, m, count, compute_com(x, y, z, m, count));

    const int32 mask_offset = first_bit;
    const int32 mask_count = last_bit - first_bit + 1;
    for (int i = 0; i < mol.atom.count - mask_count; ++i) {
        if (structure_match(mol, ref_mask, mask_offset, mask_count, i)) {
            Bitfield mask;
            bitfield::init(&mask, mol.atom.count);
            bitfield::clear_all(mask);
            for (int j = 0; j < mask_count; ++j) {
                if (ref_mask[mask_offset + j]) bitfield::set_bit(mask, i + j);
            }

            bitfield::extract_data_from_mask(x, mol.atom.position.x, mask);
            bitfield::extract_data_from_mask(y, mol.atom.position.y, mask);
            bitfield::extract_data_from_mask(z, mol.atom.position.z, mask);
            bitfield::extract_data_from_mask(m, mol.atom.mass, mask);

            const EigenFrame eigen = compute_eigen_frame(x, y, z, m, count, compute_com(x, y, z, m, count));

            const float ratio_x = ref_eigen.value[0] / eigen.value[0];
            const float ratio_y = ref_eigen.value[1] / eigen.value[1];
            const float ratio_z = ref_eigen.value[2] / eigen.value[2];

            matches.push_back(mask);
            i += count;
        }
    }

#ifdef DEBUG
    for (const auto& match : matches) {
        ASSERT(bitfield::number_of_bits_set(match) == bitfield::number_of_bits_set(ref_mask));
    }
#endif

    return matches;
}